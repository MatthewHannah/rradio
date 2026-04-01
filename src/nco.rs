/// Numerically Controlled Oscillator with optional PLL.
///
/// Provides a complex sinusoid at a configurable frequency. Can mix a real
/// input signal down to complex baseband. Supports external phase error
/// injection for decision-directed carrier tracking (à la Redsea).

use num_complex::Complex32;

pub struct Nco {
    phase: f32,
    base_freq: f32,    // radians per sample
    pll_freq: f32,     // PLL frequency correction (radians per sample)
    pll_alpha: f32,    // proportional gain
    pll_beta: f32,     // integral gain
}

impl Nco {
    /// Create a new NCO at the given frequency.
    ///
    /// * `freq_hz` — center frequency in Hz
    /// * `sample_rate` — NCO step rate in Hz
    /// * `pll_bw_hz` — PLL bandwidth in Hz (0 to disable PLL)
    /// * `feedback_rate` — rate at which step_pll() will be called (Hz)
    pub fn new(freq_hz: f32, sample_rate: f32, pll_bw_hz: f32, feedback_rate: f32) -> Self {
        let base_freq = 2.0 * std::f32::consts::PI * freq_hz / sample_rate;

        // liquid-dsp PLL formula: alpha = bw, beta = sqrt(bw)
        // bw is normalized to the feedback rate
        let bw_norm = pll_bw_hz / feedback_rate;
        let alpha = bw_norm;
        let beta = bw_norm.sqrt();

        // Scale to NCO step rate: frequency correction per step = dphi * alpha / (steps_per_feedback)
        let steps_per_feedback = sample_rate / feedback_rate;

        Nco {
            phase: 0.0,
            base_freq,
            pll_freq: 0.0,
            pll_alpha: alpha / steps_per_feedback,
            pll_beta: beta,
        }
    }

    /// Mix a real sample down to complex baseband: sample × exp(-jφ).
    #[inline]
    pub fn mix_down(&self, sample: f32) -> Complex32 {
        Complex32::new(sample * self.phase.cos(), -sample * self.phase.sin())
    }

    /// Advance the oscillator phase by one sample.
    #[inline]
    pub fn step(&mut self) {
        self.phase += self.base_freq + self.pll_freq;
        // Keep phase in [0, 2π) — avoid unbounded growth
        if self.phase >= 2.0 * std::f32::consts::PI {
            self.phase -= 2.0 * std::f32::consts::PI;
        } else if self.phase < 0.0 {
            self.phase += 2.0 * std::f32::consts::PI;
        }
    }

    /// Apply external phase error to the PLL.
    /// Matches liquid-dsp: adjust_frequency(dphi * alpha), adjust_phase(dphi * beta)
    #[inline]
    pub fn step_pll(&mut self, dphi: f32) {
        self.pll_freq += dphi * self.pll_alpha;
        self.phase += dphi * self.pll_beta;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mix_down_frequency() {
        let fs = 240000.0;
        let freq = 57000.0;
        let mut nco = Nco::new(freq, fs, 0.0, 1.0);

        // Mix a cosine at 57 kHz — should produce DC
        let mut sum = Complex32::new(0.0, 0.0);
        let n = 1000;
        for i in 0..n {
            let t = i as f32 / fs;
            let input = (2.0 * std::f32::consts::PI * freq * t).cos();
            sum += nco.mix_down(input);
            nco.step();
        }
        // DC component should be significant (positive real, near zero imaginary)
        let avg = sum / n as f32;
        assert!(avg.re > 0.3, "Expected positive DC, got {}", avg.re);
        assert!(avg.im.abs() < 0.1, "Expected near-zero imaginary, got {}", avg.im);
    }

    #[test]
    fn test_pll_tracks_offset() {
        let fs = 240000.0;
        let mut nco = Nco::new(57000.0, fs, 1.0, 1000.0);

        // Feed a signal at 57001 Hz — PLL should track the 1 Hz offset
        let actual_freq = 57001.0;
        for i in 0..24000 {
            let t = i as f32 / fs;
            let input = (2.0 * std::f32::consts::PI * actual_freq * t).cos();
            let baseband = nco.mix_down(input);
            // Simple phase error: angle of baseband (should be near 0 when locked)
            let error = baseband.im.atan2(baseband.re);
            nco.step_pll(error * 0.1);
            nco.step();
        }

        // After settling, the PLL frequency correction should be close to
        // the offset: 2π × 1 / 240000 ≈ 2.6e-5 rad/sample
        let expected = 2.0 * std::f32::consts::PI * 1.0 / fs;
        assert!(
            (nco.pll_freq - expected).abs() < expected * 0.5,
            "PLL freq {:.6} should be near {:.6}",
            nco.pll_freq, expected
        );
    }
}
