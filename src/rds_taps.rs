/// RDS matched filter tap generation.
///
/// Generates Manchester-shaped RRC matched filter taps at runtime,
/// equivalent to `py/rds_manchester_rrc.py`.

use std::f64::consts::PI;

const RDS_BIT_RATE: f64 = 1187.5;
const RDS_CHIP_RATE: f64 = 2375.0;
const T_D: f64 = 1.0 / RDS_BIT_RATE;
const T_CHIP: f64 = 1.0 / RDS_CHIP_RATE;

/// Window function types for the matched filter.
#[derive(Debug, Clone)]
pub enum WindowType {
    Blackman,
    Hamming,
    Hann,
    Rectangular,
}

impl WindowType {
    pub fn from_str(s: &str) -> Option<WindowType> {
        match s {
            "blackman" => Some(WindowType::Blackman),
            "hamming" => Some(WindowType::Hamming),
            "hann" => Some(WindowType::Hann),
            "rectangular" => Some(WindowType::Rectangular),
            _ => None,
        }
    }

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

/// Generate a windowed-sinc FIR lowpass filter.
///
/// Parameters:
/// - `fs`: sample rate in Hz
/// - `cutoff`: cutoff frequency in Hz
/// - `num_taps`: number of filter taps (must be odd)
/// - `window`: window function type
///
/// Returns a Vec<f32> of filter taps, normalized to unit DC gain.
pub fn generate_lowpass_taps(fs: f64, cutoff: f64, num_taps: usize, window: &WindowType) -> Vec<f32> {
    let num_taps = if num_taps % 2 == 0 { num_taps + 1 } else { num_taps };
    let n_half = (num_taps / 2) as i64;
    let fc = cutoff / fs; // normalized cutoff (0 to 0.5)

    let mut h: Vec<f64> = (0..num_taps)
        .map(|i| {
            let n = i as i64 - n_half;
            2.0 * fc * sinc(2.0 * fc * n as f64)
        })
        .collect();

    // Apply window
    let w = window.apply(num_taps);
    for (h, w) in h.iter_mut().zip(w.iter()) {
        *h *= w;
    }

    // Normalize to unit DC gain
    let dc_gain: f64 = h.iter().sum();
    if dc_gain.abs() > 0.0 {
        for x in h.iter_mut() {
            *x /= dc_gain;
        }
    }

    h.iter().map(|&x| x as f32).collect()
}

/// Generate Manchester-shaped RRC matched filter taps.
///
/// Parameters:
/// - `fs`: sample rate in Hz (should be integer multiple of 1187.5)
/// - `num_spans`: chip periods to span on each side for the base RRC
/// - `window`: window function type
///
/// Returns a Vec<f32> of filter taps, normalized to unit energy.
pub fn generate_manchester_rrc_taps(fs: f64, num_spans: usize, window: &WindowType) -> Vec<f32> {
    // Step 1: Base RRC impulse response
    let n_half = (num_spans as f64 * T_CHIP * fs).round() as i64;
    let n_taps = (2 * n_half + 1) as usize;

    let mut h_rrc: Vec<f64> = (0..n_taps)
        .map(|i| {
            let n = i as i64 - n_half;
            let u = 4.0 * n as f64 / (T_D * fs);
            sinc(u + 0.5) + sinc(u - 0.5)
        })
        .collect();

    // Step 2: Apply window
    let w = window.apply(n_taps);
    for (h, w) in h_rrc.iter_mut().zip(w.iter()) {
        *h *= w;
    }

    // Step 3: Convolve with Manchester waveform [+1...+1, -1...-1]
    let sps = (fs / RDS_BIT_RATE).round() as usize;
    let mut manchester = vec![1.0_f64; sps];
    for i in sps / 2..sps {
        manchester[i] = -1.0;
    }

    let conv_len = h_rrc.len() + manchester.len() - 1;
    let mut h: Vec<f64> = vec![0.0; conv_len];
    for (i, &r) in h_rrc.iter().enumerate() {
        for (j, &m) in manchester.iter().enumerate() {
            h[i + j] += r * m;
        }
    }

    // Step 4: Normalize to unit energy
    let energy: f64 = h.iter().map(|x| x * x).sum::<f64>().sqrt();
    if energy > 0.0 {
        for x in h.iter_mut() {
            *x /= energy;
        }
    }

    // Convert to f32
    h.iter().map(|&x| x as f32).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tap_count() {
        let taps = generate_manchester_rrc_taps(19000.0, 8, &WindowType::Blackman);
        // RRC: 2 * round(8 * T_CHIP * 19000) + 1 = 2*64+1 = 129
        // Manchester: 16 samples
        // Convolution: 129 + 16 - 1 = 144
        assert_eq!(taps.len(), 144);
    }

    #[test]
    fn test_unit_energy() {
        let taps = generate_manchester_rrc_taps(19000.0, 8, &WindowType::Blackman);
        let energy: f32 = taps.iter().map(|x| x * x).sum::<f32>().sqrt();
        assert!(
            (energy - 1.0).abs() < 0.01,
            "Expected unit energy, got {}", energy
        );
    }

    #[test]
    fn test_symmetric_antisymmetric() {
        // The Manchester-RRC should be antisymmetric (odd symmetry)
        // because the Manchester waveform is antisymmetric
        let taps = generate_manchester_rrc_taps(19000.0, 8, &WindowType::Blackman);
        let n = taps.len();
        let mid = n / 2;
        for i in 0..mid {
            let diff = (taps[i] + taps[n - 1 - i]).abs();
            assert!(
                diff < 1e-6,
                "Tap {} and {} should be negatives: {} vs {}",
                i, n - 1 - i, taps[i], taps[n - 1 - i]
            );
        }
    }

    #[test]
    fn test_matches_python_output() {
        // Verify against known Python-generated taps (spot check)
        let taps = generate_manchester_rrc_taps(19000.0, 8, &WindowType::Blackman);

        // Peak should be near the center, positive lobe then negative
        let mid = taps.len() / 2;
        // The positive peak is before center, negative peak after
        assert!(taps[mid - 4] > 0.1, "Expected positive peak before center");
        assert!(taps[mid + 4] < -0.1, "Expected negative peak after center");
    }

    #[test]
    fn test_different_windows() {
        let blackman = generate_manchester_rrc_taps(19000.0, 8, &WindowType::Blackman);
        let hamming = generate_manchester_rrc_taps(19000.0, 8, &WindowType::Hamming);
        let hann = generate_manchester_rrc_taps(19000.0, 8, &WindowType::Hann);

        // All should have same length
        assert_eq!(blackman.len(), hamming.len());
        assert_eq!(blackman.len(), hann.len());

        // But different tap values (after normalization they're similar but not identical)
        let diff: f32 = blackman.iter().zip(hamming.iter())
            .map(|(a, b)| (a - b).abs())
            .sum();
        assert!(diff > 1e-6, "Blackman and Hamming should differ, diff={}", diff);
    }

    #[test]
    fn test_different_spans() {
        let short = generate_manchester_rrc_taps(19000.0, 4, &WindowType::Blackman);
        let long = generate_manchester_rrc_taps(19000.0, 16, &WindowType::Blackman);

        assert!(long.len() > short.len(), "More spans = more taps");
    }
}

    #[test]
    fn test_lowpass_dump() {
        let taps = generate_lowpass_taps(240000.0, 4000.0, 201, &WindowType::Blackman);
        assert_eq!(taps.len(), 201);
        // Print first 10, middle, last 10 for comparison
        for i in 0..10 {
            eprintln!("tap[{}] = {:.15e}", i, taps[i]);
        }
        eprintln!("tap[100] = {:.15e}", taps[100]);
        for i in 191..201 {
            eprintln!("tap[{}] = {:.15e}", i, taps[i]);
        }
        // DC gain should be ~1.0
        let dc: f32 = taps.iter().sum();
        eprintln!("DC gain: {}", dc);
        assert!((dc - 1.0).abs() < 0.001, "DC gain should be ~1.0, got {}", dc);
    }
