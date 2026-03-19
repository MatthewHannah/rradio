/// Polyphase filterbank symbol synchronizer — port of liquid-dsp's symsync.
///
/// Uses a polyphase matched filter bank with timing error detection from
/// the derivative matched filter (Mengali method). This is the same
/// algorithm used in liquid-dsp's symsync_crcf, which Redsea depends on.

use num_complex::Complex32;

/// Generate Root Raised Cosine filter taps at a given sample rate.
fn generate_rrc_prototype(k: usize, m: usize, npfb: usize, beta: f32) -> Vec<f32> {
    let h_len = 2 * m * k * npfb + 1;
    let t_sym = (k * npfb) as f32; // samples per symbol in prototype rate
    let mut h = vec![0.0_f32; h_len];
    let center = (h_len - 1) / 2;

    for i in 0..h_len {
        let t = (i as f32 - center as f32) / t_sym;
        h[i] = rrc_impulse(t, beta);
    }

    // Normalize energy
    let energy: f32 = h.iter().map(|x| x * x).sum::<f32>().sqrt();
    if energy > 0.0 {
        for x in h.iter_mut() {
            *x /= energy;
        }
    }

    h
}

/// RRC impulse response at normalized time t (in symbol periods).
fn rrc_impulse(t: f32, beta: f32) -> f32 {
    let pi = std::f32::consts::PI;

    if t.abs() < 1e-7 {
        // t = 0
        1.0 + beta * (4.0 / pi - 1.0)
    } else if ((4.0 * beta * t).abs() - 1.0).abs() < 1e-4 {
        // t = ±1/(4β) singularity
        (beta / std::f32::consts::SQRT_2)
            * ((1.0 + 2.0 / pi) * (pi / (4.0 * beta)).sin()
                + (1.0 - 2.0 / pi) * (pi / (4.0 * beta)).cos())
    } else {
        let pi_t = pi * t;
        let num = (pi_t * (1.0 - beta)).sin() + 4.0 * beta * t * (pi_t * (1.0 + beta)).cos();
        let den = pi_t * (1.0 - (4.0 * beta * t).powi(2));
        num / den
    }
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

/// Decompose prototype into polyphase arms.
/// Arm k gets taps at indices k, k+npfb, k+2*npfb, ...
fn polyphase_decompose(proto: &[f32], npfb: usize) -> Vec<Vec<f32>> {
    let arm_len = (proto.len() + npfb - 1) / npfb;
    let mut banks = Vec::with_capacity(npfb);
    for k in 0..npfb {
        let arm: Vec<f32> = (0..arm_len)
            .map(|i| {
                let idx = k + i * npfb;
                if idx < proto.len() { proto[idx] } else { 0.0 }
            })
            .collect();
        banks.push(arm);
    }
    banks
}

pub struct SymSync {
    // Polyphase filter banks
    mf_banks: Vec<Vec<f32>>,      // matched filter arms
    dmf_banks: Vec<Vec<f32>>,     // derivative matched filter arms
    delay_line: Vec<Complex32>,   // circular buffer
    write_pos: usize,
    arm_len: usize,

    // Timing state
    npfb: usize,
    k: usize,
    k_out: usize,
    tau: f32,
    bf: f32,
    b: i32,
    rate: f32,
    del: f32,
    decim_counter: usize,

    // IIR loop filter state (1st-order, transposed direct form II)
    iir_b0: f32,
    iir_a1: f32,
    iir_v: f32,        // filter state
    rate_adjustment: f32,

    // Output buffer
    output_buf: Vec<Complex32>,
}

impl SymSync {
    /// Create a new polyphase filterbank symbol synchronizer.
    ///
    /// * `k` — input samples per symbol (e.g., 3)
    /// * `m` — filter semi-length in symbols (e.g., 3)
    /// * `beta` — RRC roll-off factor (e.g., 0.8)
    /// * `npfb` — number of polyphase filter banks (e.g., 32)
    pub fn new(k: usize, m: usize, beta: f32, npfb: usize) -> Self {
        // Generate prototype filters at npfb × input rate
        let proto_mf = generate_rrc_prototype(k, m, npfb, beta);
        let proto_dmf = differentiate(&proto_mf);

        // Decompose into polyphase arms
        let mf_banks = polyphase_decompose(&proto_mf, npfb);
        let dmf_banks = polyphase_decompose(&proto_dmf, npfb);
        let arm_len = mf_banks[0].len();

        let k_out = 1_usize;
        let rate = k as f32 / k_out as f32;

        SymSync {
            mf_banks,
            dmf_banks,
            delay_line: vec![Complex32::new(0.0, 0.0); arm_len],
            write_pos: 0,
            arm_len,
            npfb,
            k,
            k_out,
            tau: 0.0,
            bf: 0.0,
            b: 0,
            rate,
            del: rate,
            decim_counter: 0,
            iir_b0: 0.0,
            iir_a1: 0.0,
            iir_v: 0.0,
            rate_adjustment: 0.0,
            output_buf: Vec::with_capacity(4),
        }
    }

    /// Set the timing loop bandwidth (normalized, 0..1).
    /// Matches liquid-dsp's symsync_crcf_set_bandwidth().
    pub fn set_bandwidth(&mut self, bt: f32) {
        let alpha = 1.0 - bt;
        let beta_coeff = 0.22 * bt;
        let a = 0.5_f32;
        let b_coeff = 0.495_f32;

        // IIR SOS coefficients (normalized by a0)
        let a0 = 1.0 - a * alpha;
        self.iir_b0 = beta_coeff / a0;
        self.iir_a1 = -(b_coeff * alpha) / a0;  // negative of denominator coeff
        self.iir_v = 0.0;

        self.rate_adjustment = 0.5 * bt;
    }

    /// Set output rate (symbols per input symbol period).
    pub fn set_output_rate(&mut self, k_out: usize) {
        self.k_out = k_out;
        self.rate = self.k as f32 / k_out as f32;
        self.del = self.rate;
    }

    /// Push one input sample into the delay line.
    #[inline]
    fn push(&mut self, sample: Complex32) {
        self.delay_line[self.write_pos] = sample;
        self.write_pos = (self.write_pos + 1) % self.arm_len;
    }

    /// Evaluate polyphase arm `arm` against the delay line.
    #[inline]
    fn evaluate_arm(&self, bank: &[Vec<f32>], arm: usize) -> Complex32 {
        let taps = &bank[arm];
        let mut sum = Complex32::new(0.0, 0.0);
        for i in 0..self.arm_len {
            let idx = (self.write_pos + self.arm_len - 1 - i) % self.arm_len;
            sum += self.delay_line[idx] * taps[i];
        }
        sum
    }

    /// Process the internal timing loop (called at decimation points).
    fn advance_loop(&mut self, mf: Complex32, dmf: Complex32) {
        // Timing error: Re(conj(mf) * dmf) — Mengali method
        let q = (mf.conj() * dmf).re;
        let q = q.clamp(-1.0, 1.0);

        // IIR loop filter (1st-order, transposed direct form II)
        let q_hat = self.iir_b0 * q + self.iir_v;
        self.iir_v = -self.iir_a1 * q_hat; // a1 already negated in set_bandwidth

        // Update rate and step size
        self.rate += self.rate_adjustment * q_hat;
        self.del = self.rate + q_hat;
    }

    /// Execute the synchronizer on one input sample.
    /// Returns a slice of output symbols (0 or 1 symbols typically).
    pub fn execute(&mut self, input: Complex32) -> &[Complex32] {
        self.output_buf.clear();

        // Push into shared delay line
        self.push(input);

        // Process polyphase filterbank
        while self.b < self.npfb as i32 {
            // Compute matched filter output for arm b
            let mf = self.evaluate_arm(&self.mf_banks, self.b as usize);

            // Scale by 1/k
            let y = mf / self.k as f32;
            self.output_buf.push(y);

            // At decimation point, update timing loop
            self.decim_counter += 1;
            if self.decim_counter >= self.k_out {
                let dmf = self.evaluate_arm(&self.dmf_banks, self.b as usize);
                self.advance_loop(mf, dmf);
                self.decim_counter = 0;
            }

            // Update timing phase
            self.tau += self.del;
            self.bf = self.tau * self.npfb as f32;
            self.b = self.bf.round() as i32;
        }

        // Wrap
        self.tau -= 1.0;
        self.bf -= self.npfb as f32;
        self.b -= self.npfb as i32;

        &self.output_buf
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rrc_energy_normalized() {
        let proto = generate_rrc_prototype(3, 3, 32, 0.8);
        let energy: f32 = proto.iter().map(|x| x * x).sum::<f32>().sqrt();
        assert!((energy - 1.0).abs() < 0.01, "Energy should be ~1.0, got {}", energy);
    }

    #[test]
    fn test_polyphase_decomposition() {
        let proto = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let banks = polyphase_decompose(&proto, 3);
        assert_eq!(banks.len(), 3);
        assert_eq!(banks[0], vec![1.0, 4.0]);
        assert_eq!(banks[1], vec![2.0, 5.0]);
        assert_eq!(banks[2], vec![3.0, 6.0]);
    }

    #[test]
    fn test_output_rate() {
        let mut ss = SymSync::new(3, 3, 0.8, 32);
        ss.set_bandwidth(0.01);
        ss.set_output_rate(1);

        // Feed 300 input samples (100 symbols at 3 samp/sym)
        let mut total_outputs = 0;
        for i in 0..300 {
            let phase = 2.0 * std::f32::consts::PI * (i as f32 / 3.0);
            let input = Complex32::new(phase.cos(), 0.0);
            total_outputs += ss.execute(input).len();
        }

        // Should get approximately 100 outputs (1 per symbol)
        assert!(
            total_outputs > 80 && total_outputs < 120,
            "Expected ~100 outputs, got {}",
            total_outputs
        );
    }
}
