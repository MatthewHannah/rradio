/// Costas loop for BPSK phase synchronization.
///
/// Accepts Complex32 input, removes residual frequency offset and phase
/// error, and outputs the phase-corrected real component (f32).
/// This collapses the complex BPSK signal to a real ±amplitude stream.

use num_complex::Complex32;

pub struct CostasLoop {
    phase: f32,
    freq: f32,
    alpha: f32,
    beta: f32,
}

impl CostasLoop {
    pub fn new(loop_bw: f32) -> Self {
        // PI loop filter gains from 2nd-order PLL control theory.
        // Same formula as Gardner — see gardner_clock_recovery.rs for details.
        let damping = 0.707_f32;
        let denom = 1.0 + 2.0 * damping * loop_bw + loop_bw * loop_bw;
        let alpha = 4.0 * damping * loop_bw / denom;
        let beta = 4.0 * loop_bw * loop_bw / denom;

        CostasLoop {
            phase: 0.0,
            freq: 0.0,
            alpha,
            beta,
        }
    }

    /// Process one complex sample. Returns the phase-corrected real component.
    pub fn process(&mut self, input: Complex32) -> f32 {
        // Rotate input by negative estimated phase
        let correction = Complex32::new(self.phase.cos(), -self.phase.sin());
        let rotated = input * correction;

        // BPSK Costas error: Re × Im
        let error = rotated.re * rotated.im;

        // Update loop
        self.freq += self.beta * error;
        self.phase += self.freq + self.alpha * error;

        // Keep phase in [0, 2π)
        while self.phase >= 2.0 * std::f32::consts::PI {
            self.phase -= 2.0 * std::f32::consts::PI;
        }
        while self.phase < 0.0 {
            self.phase += 2.0 * std::f32::consts::PI;
        }

        rotated.re
    }
}

/// Iterator adapter: consumes Complex32, yields f32.
pub struct CostasIter<I: Iterator<Item = Complex32>> {
    iter: I,
    costas: CostasLoop,
}

impl<I: Iterator<Item = Complex32>> Iterator for CostasIter<I> {
    type Item = f32;

    fn next(&mut self) -> Option<f32> {
        let sample = self.iter.next()?;
        Some(self.costas.process(sample))
    }
}

pub trait CostasDemodulable {
    fn costas_demod(self, loop_bw: f32) -> CostasIter<Self>
    where
        Self: Sized + Iterator<Item = Complex32>;
}

impl<I: Iterator<Item = Complex32>> CostasDemodulable for I {
    fn costas_demod(self, loop_bw: f32) -> CostasIter<Self> {
        CostasIter {
            iter: self,
            costas: CostasLoop::new(loop_bw),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn test_removes_phase_offset() {
        let mut costas = CostasLoop::new(0.05);
        let offset = 30.0_f32.to_radians();
        let rot = Complex32::new(offset.cos(), offset.sin());

        let mut output = Vec::new();
        for i in 0..500 {
            let bit = if (i * 7 + 3) % 11 > 5 { 1.0 } else { -1.0 };
            let sample = Complex32::new(bit, 0.0) * rot;
            output.push(costas.process(sample));
        }

        let settled = &output[200..];
        let near_peak = settled.iter().filter(|&&x| x.abs() > 0.7).count();
        let ratio = near_peak as f32 / settled.len() as f32;
        assert!(ratio > 0.9, "Expected >90% near ±1, got {:.0}%", ratio * 100.0);
    }

    #[test]
    fn test_tracks_frequency_offset() {
        let mut costas = CostasLoop::new(0.05);
        let freq_offset = 0.01_f32; // radians per sample

        let mut output = Vec::new();
        for i in 0..2000 {
            let bit = if (i * 7 + 3) % 11 > 5 { 1.0 } else { -1.0 };
            let phase = freq_offset * i as f32;
            let sample = Complex32::new(bit, 0.0) * Complex32::new(phase.cos(), phase.sin());
            output.push(costas.process(sample));
        }

        let settled = &output[500..];
        let near_peak = settled.iter().filter(|&&x| x.abs() > 0.5).count();
        let ratio = near_peak as f32 / settled.len() as f32;
        assert!(ratio > 0.7, "Expected >70% near peaks, got {:.0}%", ratio * 100.0);
    }

    #[test]
    fn test_no_offset_passthrough() {
        let mut costas = CostasLoop::new(0.05);

        let symbols: Vec<f32> = (0..500)
            .map(|i| if (i * 3 + 1) % 5 > 2 { 1.0 } else { -1.0 })
            .collect();

        let output: Vec<f32> = symbols.iter()
            .map(|&s| costas.process(Complex32::new(s, 0.0)))
            .collect();

        let settled_in = &symbols[100..];
        let settled_out = &output[100..];
        let matches = settled_in.iter().zip(settled_out.iter())
            .filter(|&(&a, &b)| (a > 0.0) == (b > 0.0))
            .count();

        let accuracy = matches as f32 / settled_in.len() as f32;
        assert!(accuracy > 0.95, "Passthrough: {:.0}% match", accuracy * 100.0);
    }
}
