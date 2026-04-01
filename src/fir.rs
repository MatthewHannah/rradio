/// Developed using Claude Opus 4.6

use crate::filterable::{Filterable, Filter};
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct Fir<Num> {
    coeffs: Vec<f32>,
    delay_line: Vec<Num>,
    head: usize,
    len: usize,
}

impl<Num: Filterable<Num>> Fir<Num> {
    pub fn new(coeffs: Vec<f32>) -> Self {
        let len = coeffs.len();
        Fir {
            coeffs,
            delay_line: vec![Num::zero(); len * 2], // doubled for contiguous access
            head: 0,
            len,
        }
    }

    /// Push a sample into the delay line without computing output.
    #[inline]
    pub fn push(&mut self, x: Num) {
        self.delay_line[self.head] = x;
        self.delay_line[self.head + self.len] = x; // duplicate for wraparound
        self.head = (self.head + 1) % self.len;
    }

    /// Compute the filter output from the current delay line state.
    /// head points to the oldest sample; head+len-1 is the newest.
    #[inline]
    pub fn execute(&self) -> Num {
        self.delay_line[self.head..self.head + self.len]
            .iter()
            .zip(self.coeffs.iter().rev())
            .fold(Num::zero(), |acc, (&s, &c)| acc + s * c)
    }

    pub fn process(&mut self, x: Num) -> Num {
        self.push(x);
        self.execute()
    }
}

impl<Num> Filter<Num> for Fir<Num> where Num: Filterable<Num> {
    fn process(&mut self, x: Num) -> Num {
        self.process(x)
    }
}

// ── FIR filter design ──────────────────────────────────────────────

/// Window function types for FIR design.
#[derive(Debug, Clone)]
pub enum WindowType {
    Blackman,
    Hamming,
    Hann,
    Rectangular,
}

impl WindowType {
    fn apply(&self, n_taps: usize) -> Vec<f64> {
        match self {
            WindowType::Rectangular => vec![1.0; n_taps],
            WindowType::Blackman => (0..n_taps)
                .map(|i| {
                    let x = i as f64 / (n_taps - 1) as f64;
                    0.42 - 0.5 * (2.0 * PI * x).cos() + 0.08 * (4.0 * PI * x).cos()
                })
                .collect(),
            WindowType::Hamming => (0..n_taps)
                .map(|i| {
                    let x = i as f64 / (n_taps - 1) as f64;
                    0.54 - 0.46 * (2.0 * PI * x).cos()
                })
                .collect(),
            WindowType::Hann => (0..n_taps)
                .map(|i| {
                    let x = i as f64 / (n_taps - 1) as f64;
                    0.5 * (1.0 - (2.0 * PI * x).cos())
                })
                .collect(),
        }
    }
}

/// Normalized sinc function: sin(πx) / (πx), with sinc(0) = 1.
fn sinc(x: f64) -> f64 {
    if x.abs() < 1e-12 {
        1.0
    } else {
        (PI * x).sin() / (PI * x)
    }
}

/// Generate a windowed-sinc FIR lowpass filter, normalized to unit DC gain.
pub fn generate_lowpass_taps(fs: f64, cutoff: f64, num_taps: usize, window: &WindowType) -> Vec<f32> {
    let num_taps = if num_taps % 2 == 0 { num_taps + 1 } else { num_taps };
    let n_half = (num_taps / 2) as i64;
    let fc = cutoff / fs;

    let mut h: Vec<f64> = (0..num_taps)
        .map(|i| {
            let n = i as i64 - n_half;
            2.0 * fc * sinc(2.0 * fc * n as f64)
        })
        .collect();

    let w = window.apply(num_taps);
    for (h, w) in h.iter_mut().zip(w.iter()) {
        *h *= w;
    }

    let dc_gain: f64 = h.iter().sum();
    if dc_gain.abs() > 0.0 {
        for x in h.iter_mut() {
            *x /= dc_gain;
        }
    }

    h.iter().map(|&x| x as f32).collect()
}

/// Generate Root Raised Cosine (RRC) filter taps, normalized to unit energy.
pub fn generate_rrc_taps(fs: f64, symbol_rate: f64, beta: f64, num_symbols: usize) -> Vec<f32> {
    let sps = fs / symbol_rate;
    let half_len = (num_symbols as f64 * sps).round() as usize;
    let len = 2 * half_len + 1;
    let t_sym = 1.0 / symbol_rate;
    let mut taps = vec![0.0_f64; len];

    for i in 0..len {
        let t = (i as f64 - half_len as f64) / fs;
        let t_norm = t / t_sym;

        taps[i] = if t_norm.abs() < 1e-10 {
            (1.0 + beta * (4.0 / PI - 1.0)) / t_sym
        } else if ((4.0 * beta * t_norm).abs() - 1.0).abs() < 1e-10 {
            beta / (t_sym * std::f64::consts::SQRT_2) *
                ((1.0 + 2.0 / PI) * (PI / (4.0 * beta)).sin() +
                 (1.0 - 2.0 / PI) * (PI / (4.0 * beta)).cos())
        } else {
            let pi_t = PI * t_norm;
            let num = (pi_t * (1.0 - beta)).sin() +
                      4.0 * beta * t_norm * (pi_t * (1.0 + beta)).cos();
            let den = pi_t * (1.0 - (4.0 * beta * t_norm).powi(2));
            num / den / t_sym
        };
    }

    let energy: f64 = taps.iter().map(|x| x * x).sum::<f64>().sqrt();
    if energy > 0.0 {
        for t in taps.iter_mut() {
            *t /= energy;
        }
    }

    taps.iter().map(|&x| x as f32).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fir() {
        let coeffs = vec![0.25, 0.25, 0.25, 0.25];
        let mut fir = Fir::new(coeffs);

        let input = vec![1.0, 2.0, 3.0, 4.0];
        let expected_output = vec![0.25, 0.75, 1.5, 2.5];

        for (x, expected) in input.into_iter().zip(expected_output.into_iter()) {
            let y = fir.process(x);
            assert!((y - expected).abs() < 1e-6);
        }
    }

    #[test]
    fn test_fir_impulse_response() {
        let coeffs = vec![0.25, 0.25, 0.25, 0.25];
        let mut fir = Fir::new(coeffs);

        let input = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let expected_output = vec![0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.0];

        for (x, expected) in input.into_iter().zip(expected_output.into_iter()) {
            let y = fir.process(x);
            assert!((y - expected).abs() < 1e-6);
        }
    }

    #[test]
    fn test_fir_no_op() {
        let coeffs = vec![1.0];
        let mut fir = Fir::new(coeffs);

        let input = vec![1.0, 2.0, 3.0, 4.0];
        let expected_output = vec![1.0, 2.0, 3.0, 4.0];

        for (x, expected) in input.into_iter().zip(expected_output.into_iter()) {
            let y = fir.process(x);
            assert!((y - expected).abs() < 1e-6);
        }
    }
}
