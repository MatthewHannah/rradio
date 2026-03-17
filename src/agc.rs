/// Developed using Claude Opus 4.6
/// Automatic Gain Control (AGC).
///
/// Normalizes signal amplitude to a target RMS level using an exponential
/// moving average of signal power. Implements the `Filter<f32>` trait for
/// composability with the DSP chain.

use crate::filterable::Filter;

pub struct Agc {
    power_estimate: f32,
    alpha: f32,
    target_power: f32,
    max_gain: f32,
    initialized: bool,
}

impl Agc {
    /// Create a new AGC.
    ///
    /// - `sample_rate`: input sample rate in Hz
    /// - `bandwidth_hz`: AGC tracking bandwidth (how fast it adapts)
    /// - `target_rms`: desired output RMS amplitude
    pub fn new(sample_rate: f32, bandwidth_hz: f32, target_rms: f32) -> Self {
        let alpha = 2.0 * std::f32::consts::PI * bandwidth_hz / sample_rate;
        Agc {
            power_estimate: 0.0,
            alpha: alpha.clamp(0.0, 1.0),
            target_power: target_rms * target_rms,
            max_gain: 100.0,
            initialized: false,
        }
    }
}

impl Filter<f32> for Agc {
    fn process(&mut self, x: f32) -> f32 {
        // Seed the power estimate from the first nonzero sample
        if !self.initialized && x != 0.0 {
            self.power_estimate = x * x;
            self.initialized = true;
        }

        self.power_estimate += self.alpha * (x * x - self.power_estimate);

        let gain = if self.power_estimate > 0.0 {
            (self.target_power / self.power_estimate).sqrt().min(self.max_gain)
        } else {
            1.0 // pass through until we have signal
        };

        x * gain
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_converges_to_target() {
        let mut agc = Agc::new(9600.0, 50.0, 1.0);
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
        let mut agc = Agc::new(9600.0, 50.0, 1.0);
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
        let mut agc = Agc::new(9600.0, 50.0, 1.0);
        for _ in 0..1000 {
            let out = agc.process(0.0);
            assert!(out.is_finite(), "Output should be finite, got {}", out);
            assert_eq!(out, 0.0);
        }
    }

    #[test]
    fn test_tracks_amplitude_step() {
        let mut agc = Agc::new(9600.0, 50.0, 1.0);
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
            (last_output - 1.0).abs() < 0.15,
            "After step, expected output near 1.0, got {}", last_output
        );
    }
}
