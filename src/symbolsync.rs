use num_complex::Complex32;

use crate::pi_loop::PiLoopFilter;

// Timing loop filter mechanics derived from GNU Radio clock_tracking_loop
struct TimingLoop {
    pi_loop: PiLoopFilter,

    // accumulated resampler period (in input sample periods)
    resampler_period: f32,
    interp_factor: f32,

    // coarse and fine adjustments for resampling
    coarse: isize,
    fine: usize,
}

impl TimingLoop {
    fn new(loop_bw: f32, damping: f32, k_ted: f32, expected_period: f32, max_deviation: f32, interp_factor: f32) -> Self {
        TimingLoop {
            pi_loop: PiLoopFilter::new(loop_bw, damping, k_ted, expected_period - max_deviation, expected_period + max_deviation),
            resampler_period: 0.0,
            interp_factor,
            coarse: 0,
            fine: 0,
        }
    }

    fn advance(&mut self, error: f32) {
        // PI loop adjustments
        let mut instantaneous_period = self.pi_loop.advance(error);

        if instantaneous_period < 0.0 {
            instantaneous_period = self.pi_loop.integrator; // bad condition, reset ourselves to expected period
        }

        // coarse+interpolator adjustments
        self.coarse = 0;
        self.resampler_period += instantaneous_period;
        while self.resampler_period >= self.interp_factor {
            self.resampler_period -= self.interp_factor;
            self.coarse += 1;
        }
        while self.resampler_period <= 0.0 {
            self.resampler_period += self.interp_factor;
            self.coarse -= 1;
        }
        self.fine = (self.resampler_period.round() as usize).clamp(0, self.interp_factor as usize - 1);
    }
}

pub struct SymbolSync {
    circbuf: Vec<Complex32>,
    head: usize,
    circ_len: usize,
    samples_per_symbol: usize,
    filt: Vec<Vec<f32>>,
    dfilt: Vec<Vec<f32>>,
    timing_loop: TimingLoop,
}

/// Create polyphase filter bank from prototype filter
fn polyphase_bank(prototype_filt: &[f32], k: usize) -> Vec<Vec<f32>> {
    let arm_len = (prototype_filt.len() + k - 1) / k; // ceiling division
    let mut bank = vec![vec![0.0; arm_len]; k];
    for (i, &coef) in prototype_filt.iter().enumerate() {
        bank[i % k][i / k] = coef;
    }
    // Reverse each arm so that a forward dot product with oldest→newest
    // buffer order implements standard convolution (h[0] × newest sample).
    for arm in &mut bank {
        arm.reverse();
    }
    bank
}

/// Compute derivative filter from prototype via central difference
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

impl SymbolSync {
    pub fn new(prototype_filt: Vec<f32>, samples_per_symbol: usize, max_deviation: f32, interp_factor: usize, loop_bw: f32, damping: f32, k_ted: f32) -> Self {
        let mut derivative_filt = differentiate(&prototype_filt);
        derivative_filt.iter_mut().for_each(|x| *x *= interp_factor as f32); // scale prototype for interpolation

        let filt = polyphase_bank(&prototype_filt, interp_factor);
        let dfilt = polyphase_bank(&derivative_filt, interp_factor);

        SymbolSync {
            circbuf: vec![Complex32::new(0.0, 0.0); filt[0].len() * 2], // provides contiguous access to samples
            head: 0,
            circ_len: filt[0].len(),
            samples_per_symbol,
            filt,
            dfilt,
            timing_loop: TimingLoop::new(loop_bw, damping, k_ted, samples_per_symbol as f32, max_deviation, interp_factor as f32),
        }
    }

    pub fn next(&mut self, iter: &mut impl Iterator<Item = Complex32>) -> Option<Complex32> {
        let consume = self.samples_per_symbol as isize + self.timing_loop.coarse;
        for _ in 0..consume {
            let sample = iter.next()?;
            self.circbuf[self.head] = sample;
            self.circbuf[self.head + self.circ_len] = sample; // duplicate for wraparound
            self.head = (self.head + 1) % self.circ_len;
        }

        let arm_idx = self.timing_loop.fine;
        let sym_estimate = filter_eval(&self.circbuf, self.head, &self.filt[arm_idx]);
        let sym_derivative = filter_eval(&self.circbuf, self.head, &self.dfilt[arm_idx]);
        let error = SymbolSync::compute_error(sym_estimate, sym_derivative);

        self.timing_loop.advance(error);

        Some(sym_estimate)
    }

    fn compute_error(sym: Complex32, sym_derivative: Complex32) -> f32 {
        (sym.conj() * sym_derivative).re / 2.0
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