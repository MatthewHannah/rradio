use num_complex::Complex32;

/// Power-averaging AGC targeting unit RMS amplitude.
/// Initial acquisition phase accumulates samples to estimate gain,
/// then switches to exponential power averaging.
pub struct PowerAgc {
    rate: f32,
    target_power: f32,
    avg_power: f32,
    pub gain: f32,
    max_gain: f32,
    // Acquisition phase
    acquired: bool,
    acq_accum: f32,
    acq_count: usize,
    acq_threshold: usize,
}

impl PowerAgc {
    pub fn new(rate: f32, target_rms: f32, max_gain: f32, acq_samples: usize) -> Self {
        PowerAgc {
            rate,
            target_power: target_rms * target_rms,
            avg_power: 1.0,
            gain: 1.0,
            max_gain,
            acquired: false,
            acq_accum: 0.0,
            acq_count: 0,
            acq_threshold: acq_samples,
        }
    }

    pub fn process(&mut self, sample: Complex32) -> Complex32 {
        if !self.acquired {
            let power = sample.re * sample.re + sample.im * sample.im;
            self.acq_accum += power;
            self.acq_count += 1;
            if self.acq_count >= self.acq_threshold {
                self.avg_power = self.acq_accum / self.acq_count as f32;
                if self.avg_power > 1e-20 {
                    self.gain = (self.target_power / self.avg_power).sqrt().min(self.max_gain);
                }
                self.acquired = true;
            }
            sample
        } else {
            let power = sample.re * sample.re + sample.im * sample.im;
            self.avg_power = (1.0 - self.rate) * self.avg_power + self.rate * power;
            if self.avg_power > 1e-20 {
                self.gain = (self.target_power / self.avg_power).sqrt().min(self.max_gain);
            }
            sample * self.gain
        }
    }
}

/// Create polyphase filter bank from prototype filter.
/// Returns Vec of `nfilters` arms, each with `ceil(len/nfilters)` taps.
/// Arms are reversed for forward dot-product convolution.
fn polyphase_bank(prototype_filt: &[f32], k: usize) -> Vec<Vec<f32>> {
    let arm_len = (prototype_filt.len() + k - 1) / k;
    let mut bank = vec![vec![0.0; arm_len]; k];
    for (i, &coef) in prototype_filt.iter().enumerate() {
        bank[i % k][i / k] = coef;
    }
    for arm in &mut bank {
        arm.reverse();
    }
    bank
}

/// Compute derivative filter from prototype via central difference.
fn differentiate(h: &[f32]) -> Vec<f32> {
    let n = h.len();
    let mut dh = vec![0.0_f32; n];
    for i in 0..n {
        let prev = if i > 0 { h[i - 1] } else { 0.0 };
        let next = if i < n - 1 { h[i + 1] } else { 0.0 };
        dh[i] = 0.5 * (next - prev);
    }
    dh
}

fn filter_eval(circ: &[Complex32], oldest: usize, coeffs: &[f32]) -> Complex32 {
    circ[oldest..oldest + coeffs.len()]
        .iter()
        .zip(coeffs)
        .map(|(&s, &c)| s * c)
        .sum()
}

/// Polyphase filterbank symbol synchronizer.
///
/// Generalized clock tracking approach following GNU Radio's symbol_sync:
/// phase accumulates at the INPUT sample rate. A chip boundary fires when
/// the phase wraps past `nfilters`. The fractional phase at the wrap point
/// selects the polyphase arm for sub-sample interpolation.
///
/// This handles any SPS (integer or fractional) and naturally absorbs
/// clock drift into the instantaneous period estimate.
///
/// The PI loop filter tracks the symbol period in input samples
/// (clock_tracking_loop formulation). K_TED is in units of TED output
/// per input sample of timing offset.
pub struct SymbolSync {
    // Polyphase filter bank
    circbuf: Vec<Complex32>,
    head: usize,
    circ_len: usize,
    nfilters: usize,
    filt: Vec<Vec<f32>>,
    dfilt: Vec<Vec<f32>>,

    // Phase accumulator (runs at input sample rate)
    phase: f32,         // [0, nfilters), fractional arm position
    inst_rate: f32,     // phase advance per input sample (= nfilters / inst_period)

    // PI loop filter (period-based, clock_tracking_loop style)
    avg_period: f32,    // integrator: average period in input samples
    inst_period: f32,   // instantaneous period (avg + proportional correction)
    alpha: f32,         // proportional gain
    beta: f32,          // integral gain
    min_period: f32,
    max_period: f32,

    pub agc: PowerAgc,
}

impl SymbolSync {
    /// Create a new SymbolSync.
    ///
    /// * `prototype_filt` — RRC prototype at `nfilters × sps` rate
    /// * `sps` — nominal samples per chip at input rate
    /// * `max_period_deviation` — max period deviation in input samples
    /// * `nfilters` — number of polyphase arms
    /// * `loop_bw` — ωn_norm (normalized natural frequency at chip rate)
    /// * `damping` — ζ
    /// * `k_ted` — TED gain per input sample of timing offset (= K_TED_acc × nfilters)
    pub fn new(prototype_filt: Vec<f32>, sps: usize, max_period_deviation: f32,
               nfilters: usize, loop_bw: f32, damping: f32, k_ted_per_sample: f32) -> Self {
        let mut derivative_filt = differentiate(&prototype_filt);
        derivative_filt.iter_mut().for_each(|x| *x *= nfilters as f32);

        let filt = polyphase_bank(&prototype_filt, nfilters);
        let dfilt = polyphase_bank(&derivative_filt, nfilters);

        let nominal_period = sps as f32;
        let (alpha, beta) = crate::pi_loop::calculate_gains(loop_bw, damping, k_ted_per_sample);

        SymbolSync {
            circbuf: vec![Complex32::new(0.0, 0.0); filt[0].len() * 2],
            head: 0,
            circ_len: filt[0].len(),
            nfilters,
            filt,
            dfilt,

            phase: nfilters as f32 / 2.0,
            inst_rate: nfilters as f32 / nominal_period,
            avg_period: nominal_period,
            inst_period: nominal_period,
            alpha,
            beta,
            min_period: nominal_period - max_period_deviation,
            max_period: nominal_period + max_period_deviation,

            agc: PowerAgc::new(0.01, 1.0, 1e5, 100),
        }
    }

    /// Pull samples from the iterator until a chip boundary fires.
    /// Returns the AGC-normalized matched filter output, or None on end-of-stream.
    pub fn next(&mut self, iter: &mut impl Iterator<Item = Complex32>) -> Option<Complex32> {
        loop {
            let sample = iter.next()?;

            // Push into circular delay line
            self.circbuf[self.head] = sample;
            self.circbuf[self.head + self.circ_len] = sample;
            self.head = (self.head + 1) % self.circ_len;

            // Advance phase accumulator
            self.phase += self.inst_rate;

            // Check for chip boundary
            if self.phase >= self.nfilters as f32 {
                self.phase -= self.nfilters as f32;

                // Select arm from fractional phase at the wrap point
                let arm = (self.phase as usize).min(self.nfilters - 1);

                // Evaluate matched filter and derivative at selected arm
                let mf_out = filter_eval(&self.circbuf, self.head, &self.filt[arm]);
                let dmf_out = filter_eval(&self.circbuf, self.head, &self.dfilt[arm]);

                // AGC: normalize to unit power, apply same gain to dMF
                let mf_out = self.agc.process(mf_out);
                let dmf_out = dmf_out * self.agc.gain;

                // ML TED: e = Re{conj(MF) · dMF} / 2
                let error = (mf_out.conj() * dmf_out).re / 2.0;

                // PI loop filter (period-based, clock_tracking_loop style)
                self.avg_period += self.beta * error;
                self.avg_period = self.avg_period.clamp(self.min_period, self.max_period);
                self.inst_period = self.avg_period + self.alpha * error;
                if self.inst_period <= 0.0 {
                    self.inst_period = self.avg_period;
                }

                // Update phase advance rate for next chip interval
                self.inst_rate = self.nfilters as f32 / self.inst_period;

                return Some(mf_out);
            }
        }
    }
}

pub struct SymbolSyncIter<I> where I: Iterator<Item = Complex32> {
    symsync: SymbolSync,
    iter: I,
}


#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f32) -> Complex32 {
        Complex32::new(re, 0.0)
    }

    #[test]
    fn test_circ_buf_filt_arb() {
        let coeff = vec![1.0, 2.0, 3.0];
        let circ = vec![c(4.0), c(5.0), c(6.0), c(4.0), c(5.0), c(6.0)];

        let result = filter_eval(&circ, 0, &coeff);
        assert_eq!(result, c(32.0));
    }
}