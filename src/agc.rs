/// Automatic Gain Control — liquid-dsp compatible.
///
/// Uses the liquid-dsp AGC algorithm:
///   y = x * g
///   power_estimate += alpha * (|y|² - power_estimate)
///   g *= exp(-0.5 * alpha * log(power_estimate))
///
/// Supports both real (f32) and complex (Complex32) operation.

use crate::filterable::Filter;
use num_complex::Complex32;

pub struct Agc {
    gain: f32,
    power_estimate: f32,
    alpha: f32,
    scale: f32,
}

impl Agc {
    /// Create AGC with liquid-dsp compatible parameters.
    ///
    /// * `alpha` — normalized bandwidth (0..1), controls tracking speed
    /// * `initial_gain` — starting gain value
    pub fn new(alpha: f32, initial_gain: f32) -> Self {
        Agc {
            gain: initial_gain,
            power_estimate: 1.0,
            alpha: alpha.clamp(0.0, 1.0),
            scale: 1.0,
        }
    }

    #[inline]
    fn update_gain(&mut self, output_power: f32) {
        // Smooth power estimate
        self.power_estimate += self.alpha * (output_power - self.power_estimate);

        // liquid-dsp gain update: pull output power toward 1.0
        if self.power_estimate > 1e-16 {
            self.gain *= (-0.5 * self.alpha * self.power_estimate.ln()).exp();
        }

        // Clamp gain
        self.gain = self.gain.min(1e6);
    }
}

impl Filter<f32> for Agc {
    fn process(&mut self, x: f32) -> f32 {
        let y = x * self.gain;
        self.update_gain(y * y);
        y * self.scale
    }
}

impl Filter<Complex32> for Agc {
    fn process(&mut self, x: Complex32) -> Complex32 {
        let y = x * self.gain;
        self.update_gain(y.norm_sqr());
        y * self.scale
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_converges_to_target() {
        let mut agc = Agc::new(0.03, 1.0);
        let input_amp = 0.1;
        let mut last_output = 0.0;
        for i in 0..5000 {
            let x = input_amp * if i % 4 < 2 { 1.0 } else { -1.0 };
            last_output = agc.process(x).abs();
        }
        assert!(
            (last_output - 1.0).abs() < 0.15,
            "Expected output near 1.0, got {}", last_output
        );
    }

    #[test]
    fn test_attenuates_strong_signal() {
        let mut agc = Agc::new(0.03, 1.0);
        let input_amp = 10.0;
        let mut last_output = 0.0;
        for i in 0..5000 {
            let x = input_amp * if i % 4 < 2 { 1.0 } else { -1.0 };
            last_output = agc.process(x).abs();
        }
        assert!(
            (last_output - 1.0).abs() < 0.15,
            "Expected output near 1.0, got {}", last_output
        );
    }

    #[test]
    fn test_handles_zero_input() {
        let mut agc = Agc::new(0.03, 1.0);
        for _ in 0..1000 {
            let out = agc.process(0.0);
            assert!(out.is_finite(), "Output should be finite, got {}", out);
            assert_eq!(out, 0.0);
        }
    }

    #[test]
    fn test_tracks_amplitude_step() {
        let mut agc = Agc::new(0.03, 1.0);
        for i in 0..5000 {
            let x = 0.5 * if i % 4 < 2 { 1.0 } else { -1.0 };
            agc.process(x);
        }
        let mut last_output = 0.0;
        for i in 0..5000 {
            let x = 5.0 * if i % 4 < 2 { 1.0 } else { -1.0 };
            last_output = agc.process(x).abs();
        }
        assert!(
            (last_output - 1.0).abs() < 0.3,
            "After step, expected output near 1.0, got {}", last_output
        );
    }
}
