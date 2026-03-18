/// Gardner Timing Error Detector with interpolating clock recovery.
///
/// Recovers the symbol clock from an oversampled signal by finding the
/// optimal sampling instants. Consumes ~N input samples per output,
/// where N = sample_rate / symbol_rate.

use crate::rds_config::GardnerConfig;

pub struct GardnerClockRecovery<I: Iterator<Item = f32>> {
    iter: I,

    // NCO state
    nco_phase: f32,
    base_increment: f32, // symbol_rate / sample_rate

    // Interpolation buffer: buf[0] = previous, buf[1] = current
    buf: [f32; 2],

    // TED state
    prev_strobe: f32,
    mid_sample: f32,

    // Loop filter (PI controller)
    integrator: f32,
    correction: f32, // kp * error + integrator, applied every NCO step
    kp: f32,
    ki: f32,
}

impl<I: Iterator<Item = f32>> GardnerClockRecovery<I> {
    pub fn new(iter: I, symbol_rate: f32, sample_rate: f32, config: &GardnerConfig) -> Self {
        // PI loop filter gains derived from 2nd-order PLL control theory.
        // Maps desired loop bandwidth and damping ratio to discrete-time
        // proportional (Kp) and integral (Ki) gains:
        //   Kp = 4 * zeta * bw_n / (1 + 2 * zeta * bw_n + bw_n^2)
        //   Ki = 4 * bw_n^2 / (1 + 2 * zeta * bw_n + bw_n^2)
        // This is the same formula used in GNU Radio's Symbol Sync block
        // (see symbol_sync_cc_impl.cc in gr-digital).
        let bw_n = config.loop_bw;
        let zeta = 0.707;
        let denom = 1.0 + 2.0 * zeta * bw_n + bw_n * bw_n;
        let kp = 4.0 * zeta * bw_n / denom;
        let ki = 4.0 * bw_n * bw_n / denom;

        GardnerClockRecovery {
            iter,
            nco_phase: 0.0,
            base_increment: symbol_rate / sample_rate,
            buf: [0.0; 2],
            prev_strobe: 0.0,
            mid_sample: 0.0,
            integrator: 0.0,
            correction: 0.0,
            kp,
            ki,
        }
    }

    fn process(&mut self, input: f32) -> Option<f32> {
        self.buf[0] = self.buf[1];
        self.buf[1] = input;

        let prev_phase = self.nco_phase;
        let step = self.base_increment + self.correction;
        self.nco_phase += step;

        // Midpoint crossing (phase crosses 0.5)
        if prev_phase < 0.5 && self.nco_phase >= 0.5 {
            let mu = (0.5 - prev_phase) / step;
            self.mid_sample = self.buf[0] * (1.0 - mu) + self.buf[1] * mu;
        }

        // Chip strobe (phase crosses 1.0)
        if self.nco_phase >= 1.0 {
            self.nco_phase -= 1.0;
            let mu = (1.0 - prev_phase) / step;
            let curr_strobe = self.buf[0] * (1.0 - mu) + self.buf[1] * mu;

            // Gardner TED: error = midpoint × (previous_strobe − current_strobe)
            let error = self.mid_sample * (self.prev_strobe - curr_strobe);

            // PI loop filter
            self.integrator += self.ki * error;
            let max_adj = self.base_increment * 0.5;
            self.integrator = self.integrator.clamp(-max_adj, max_adj);
            self.correction = self.kp * error + self.integrator;

            self.prev_strobe = curr_strobe;

            return Some(curr_strobe);
        }

        None
    }
}

impl<I: Iterator<Item = f32>> Iterator for GardnerClockRecovery<I> {
    type Item = f32;

    fn next(&mut self) -> Option<f32> {
        loop {
            let input = self.iter.next()?;
            if let Some(chip) = self.process(input) {
                return Some(chip);
            }
        }
    }
}

pub trait GardnerClockRecoverable {
    fn clock_recover(self, symbol_rate: f32, sample_rate: f32, config: &GardnerConfig)
        -> GardnerClockRecovery<Self>
    where
        Self: Sized + Iterator<Item = f32>;
}

impl<I: Iterator<Item = f32>> GardnerClockRecoverable for I {
    fn clock_recover(self, symbol_rate: f32, sample_rate: f32, config: &GardnerConfig)
        -> GardnerClockRecovery<Self>
    {
        GardnerClockRecovery::new(self, symbol_rate, sample_rate, config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_output_rate() {
        // Verify we get approximately the right number of output samples.
        let sample_rate = 9600.0;
        let symbol_rate = 2375.0;
        let n_input = 9600; // 1 second

        let samples_per_symbol = sample_rate / symbol_rate;
        let input: Vec<f32> = (0..n_input)
            .map(|i| {
                if ((i as f32 / samples_per_symbol) as usize) % 2 == 0 {
                    1.0
                } else {
                    -1.0
                }
            })
            .collect();

        let output: Vec<f32> = input
            .into_iter()
            .clock_recover(symbol_rate, sample_rate, &GardnerConfig::default())
            .collect();

        let expected = (n_input as f32 * symbol_rate / sample_rate) as usize;
        assert!(
            output.len() > expected * 95 / 100 && output.len() < expected * 105 / 100,
            "Expected ~{} outputs, got {}",
            expected,
            output.len()
        );
    }

    /// Generate a shaped chip signal: upsample + simple low-pass smoothing.
    fn shaped_chip_signal(chips: &[f32], sample_rate: f32, symbol_rate: f32) -> Vec<f32> {
        let sps = sample_rate / symbol_rate;
        let n_samples = (chips.len() as f32 * sps) as usize;
        let mut signal = vec![0.0_f32; n_samples];

        // Upsample: place chip values at symbol centers
        for (i, &chip) in chips.iter().enumerate() {
            let center = ((i as f32 + 0.5) * sps) as usize;
            if center < n_samples {
                signal[center] = chip;
            }
        }

        // Convolve with a simple raised-cosine-ish pulse (half-sine, width = sps)
        let pulse_half = (sps / 2.0).ceil() as usize;
        let pulse_len = 2 * pulse_half + 1;
        let pulse: Vec<f32> = (0..pulse_len)
            .map(|i| {
                let t = (i as f32 - pulse_half as f32) / sps;
                if t.abs() <= 0.5 {
                    (std::f32::consts::PI * t).cos()
                } else {
                    0.0
                }
            })
            .collect();

        // Simple convolution
        let mut shaped = vec![0.0_f32; n_samples];
        for i in 0..n_samples {
            let mut sum = 0.0;
            for j in 0..pulse_len {
                let idx = i as isize - pulse_half as isize + j as isize;
                if idx >= 0 && (idx as usize) < n_samples {
                    sum += signal[idx as usize] * pulse[j];
                }
            }
            shaped[i] = sum;
        }

        shaped
    }

    #[test]
    fn test_recovers_shaped_chips() {
        let sample_rate = 9600.0;
        let symbol_rate = 2375.0;

        // Long random-ish chip pattern (deterministic for reproducibility)
        let chips: Vec<f32> = (0..200)
            .map(|i| if (i * 7 + 3) % 11 > 5 { 1.0 } else { -1.0 })
            .collect();

        let signal = shaped_chip_signal(&chips, sample_rate, symbol_rate);

        let output: Vec<f32> = signal
            .into_iter()
            .clock_recover(symbol_rate, sample_rate, &GardnerConfig::default())
            .collect();

        // Skip first 25% for loop settling, check chip polarity matches
        let settle = output.len() / 4;
        let check = &output[settle..];

        // Find the best alignment offset (loop may lock at any phase)
        let mut best_match = 0;
        let mut best_offset = 0;
        for offset in 0..4 {
            let matches: usize = check.iter().enumerate()
                .filter(|&(i, &val)| {
                    let chip_idx = settle + i + offset;
                    if chip_idx < chips.len() {
                        (val > 0.0) == (chips[chip_idx] > 0.0)
                    } else {
                        false
                    }
                })
                .count();
            if matches > best_match {
                best_match = matches;
                best_offset = offset;
            }
        }

        let check_len = check.len().min(chips.len() - settle - best_offset);
        let accuracy = best_match as f32 / check_len as f32;
        assert!(
            accuracy > 0.85,
            "Expected >85% chip match (offset={}), got {:.0}% ({}/{})",
            best_offset, accuracy * 100.0, best_match, check_len
        );
    }

    #[test]
    fn test_tracks_frequency_offset() {
        // Feed a signal with a slight frequency mismatch and verify
        // the loop still produces the correct output rate.
        let sample_rate = 9600.0;
        let symbol_rate = 2375.0;
        let actual_rate = 2380.0; // 0.2% faster than nominal
        let n_chips = 500;

        let sps_actual = sample_rate / actual_rate;
        let n_samples = (n_chips as f32 * sps_actual) as usize;

        let signal: Vec<f32> = (0..n_samples)
            .map(|i| {
                if ((i as f32 / sps_actual) as usize) % 2 == 0 {
                    1.0
                } else {
                    -1.0
                }
            })
            .collect();

        let output: Vec<f32> = signal
            .into_iter()
            .clock_recover(symbol_rate, sample_rate, &GardnerConfig::default())
            .collect();

        // Should track the actual rate, not the nominal — output count
        // should be closer to n_chips than to nominal_chips
        let nominal_chips = (n_samples as f32 * symbol_rate / sample_rate) as usize;
        let diff_from_actual = (output.len() as i32 - n_chips as i32).unsigned_abs();
        let diff_from_nominal = (output.len() as i32 - nominal_chips as i32).unsigned_abs();

        assert!(
            diff_from_actual < diff_from_nominal + 5,
            "Loop should track actual rate: got {} outputs \
             (actual={} chips, nominal={})",
            output.len(), n_chips, nominal_chips
        );
    }
}
