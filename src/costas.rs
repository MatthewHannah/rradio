/// Costas loop for BPSK phase synchronization.
///
/// Accepts Complex32 input, removes residual frequency offset and phase
/// error, and outputs the full phase-corrected complex sample.

use num_complex::Complex32;
use crate::rds_config::CostasConfig;

pub struct CostasLoop {
    phase: f32,
    freq: f32,
    alpha: f32,
    beta: f32,
    // Lock detection: exponential average of |error|
    avg_error: f32,
    lock_alpha: f32,
    sample_count: u64,
}

impl CostasLoop {
    pub fn new(config: &CostasConfig) -> Self {
        let damping = 0.707_f32;
        let denom = 1.0 + 2.0 * damping * config.loop_bw + config.loop_bw * config.loop_bw;
        let alpha = 4.0 * damping * config.loop_bw / denom;
        let beta = 4.0 * config.loop_bw * config.loop_bw / denom;

        CostasLoop {
            phase: 0.0,
            freq: 0.0,
            alpha,
            beta,
            avg_error: 1.0,
            lock_alpha: 0.01,  // smoothing for lock indicator
            sample_count: 0,
        }
    }

    /// Process one complex sample. Returns the phase-corrected complex output.
    pub fn process(&mut self, input: Complex32) -> Complex32 {
        let correction = Complex32::new(self.phase.cos(), -self.phase.sin());
        let rotated = input * correction;

        // BPSK Costas error: Re × Im (zero when aligned to real axis)
        let error = rotated.re * rotated.im;

        // Update lock indicator
        self.avg_error += self.lock_alpha * (error.abs() - self.avg_error);
        self.sample_count += 1;

        // Update loop
        self.freq += self.beta * error;
        self.phase += self.freq + self.alpha * error;

        while self.phase >= 2.0 * std::f32::consts::PI {
            self.phase -= 2.0 * std::f32::consts::PI;
        }
        while self.phase < 0.0 {
            self.phase += 2.0 * std::f32::consts::PI;
        }

        rotated
    }

    /// Average absolute error — small when locked.
    pub fn avg_error(&self) -> f32 { self.avg_error }

    /// Current frequency estimate (radians/sample).
    pub fn freq(&self) -> f32 { self.freq }

    /// Current phase estimate (radians).
    pub fn phase(&self) -> f32 { self.phase }

    /// Samples processed so far.
    pub fn sample_count(&self) -> u64 { self.sample_count }
}

/// Iterator adapter: consumes Complex32, yields Complex32.
/// Logs lock diagnostics periodically to stderr.
pub struct CostasIter<I: Iterator<Item = Complex32>> {
    iter: I,
    costas: CostasLoop,
    log_interval: u64,
}

impl<I: Iterator<Item = Complex32>> Iterator for CostasIter<I> {
    type Item = Complex32;

    fn next(&mut self) -> Option<Complex32> {
        let sample = self.iter.next()?;
        let out = self.costas.process(sample);

        if self.costas.sample_count % self.log_interval == 0 && self.costas.sample_count > 0 {
            eprintln!("COSTAS n={} avg_err={:.6} freq={:.6} rad/samp phase={:.3} rad",
                self.costas.sample_count,
                self.costas.avg_error(),
                self.costas.freq(),
                self.costas.phase());
        }

        Some(out)
    }
}

pub trait CostasDemodulable {
    fn costas_demod(self, config: &CostasConfig) -> CostasIter<Self>
    where
        Self: Sized + Iterator<Item = Complex32>;
}

impl<I: Iterator<Item = Complex32>> CostasDemodulable for I {
    fn costas_demod(self, config: &CostasConfig) -> CostasIter<Self> {
        CostasIter {
            iter: self,
            costas: CostasLoop::new(config),
            log_interval: 500,  // every 500 symbols ≈ 0.42s
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_removes_phase_offset() {
        let mut costas = CostasLoop::new(&CostasConfig::default());
        let offset = 30.0_f32.to_radians();
        let rot = Complex32::new(offset.cos(), offset.sin());

        let mut output = Vec::new();
        for i in 0..500 {
            let bit = if (i * 7 + 3) % 11 > 5 { 1.0 } else { -1.0 };
            let sample = Complex32::new(bit, 0.0) * rot;
            output.push(costas.process(sample).re);
        }

        let settled = &output[200..];
        let near_peak = settled.iter().filter(|&&x| x.abs() > 0.7).count();
        let ratio = near_peak as f32 / settled.len() as f32;
        assert!(ratio > 0.9, "Expected >90% near ±1, got {:.0}%", ratio * 100.0);
    }

    #[test]
    fn test_tracks_frequency_offset() {
        let mut costas = CostasLoop::new(&CostasConfig::default());
        let freq_offset = 0.01_f32; // radians per sample

        let mut output = Vec::new();
        for i in 0..2000 {
            let bit = if (i * 7 + 3) % 11 > 5 { 1.0 } else { -1.0 };
            let phase = freq_offset * i as f32;
            let sample = Complex32::new(bit, 0.0) * Complex32::new(phase.cos(), phase.sin());
            output.push(costas.process(sample).re);
        }

        let settled = &output[500..];
        let near_peak = settled.iter().filter(|&&x| x.abs() > 0.5).count();
        let ratio = near_peak as f32 / settled.len() as f32;
        assert!(ratio > 0.7, "Expected >70% near peaks, got {:.0}%", ratio * 100.0);
    }

    #[test]
    fn test_no_offset_passthrough() {
        let mut costas = CostasLoop::new(&CostasConfig::default());

        let symbols: Vec<f32> = (0..500)
            .map(|i| if (i * 3 + 1) % 5 > 2 { 1.0 } else { -1.0 })
            .collect();

        let output: Vec<f32> = symbols.iter()
            .map(|&s| costas.process(Complex32::new(s, 0.0)).re)
            .collect();

        let settled_in = &symbols[100..];
        let settled_out = &output[100..];
        let matches = settled_in.iter().zip(settled_out.iter())
            .filter(|&(&a, &b)| (a > 0.0) == (b > 0.0))
            .count();

        let accuracy = matches as f32 / settled_in.len() as f32;
        assert!(accuracy > 0.95, "Passthrough: {:.0}% match", accuracy * 100.0);
    }

    #[test]
    fn test_imaginary_near_zero_after_lock() {
        let mut costas = CostasLoop::new(&CostasConfig::default());
        let offset = 45.0_f32.to_radians();
        let rot = Complex32::new(offset.cos(), offset.sin());

        let mut output = Vec::new();
        for i in 0..500 {
            let bit = if (i * 7 + 3) % 11 > 5 { 1.0 } else { -1.0 };
            let sample = Complex32::new(bit, 0.0) * rot;
            output.push(costas.process(sample));
        }

        let settled = &output[200..];
        let avg_im = settled.iter().map(|s| s.im.abs()).sum::<f32>() / settled.len() as f32;
        assert!(avg_im < 0.15, "Imaginary should be near zero after lock, got {:.3}", avg_im);
    }
}
