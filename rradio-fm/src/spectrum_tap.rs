use num_complex::Complex32;
use rustfft::FftPlanner;

use crate::telemetry::LatestWriter;

/// Snapshot of a magnitude spectrum for display.
#[derive(Clone)]
pub struct SpectrumSnapshot {
    /// Magnitude spectrum in dB.
    /// For complex input: fft_size bins, FFT-shifted so DC is centered.
    /// For real input: fft_size/2 bins, 0 Hz to fs/2.
    pub magnitudes_db: Vec<f32>,
    /// Sample rate of the input signal (for computing frequency axis).
    pub sample_rate: f32,
    /// True if this came from complex (IQ) input — spectrum spans -fs/2 to +fs/2.
    /// False if from real input — spectrum spans 0 to fs/2.
    pub is_complex: bool,
}

/// Trait abstracting over f32 and Complex32 for the spectrum tap.
pub trait FftSample: Copy + Send + 'static {
    const IS_COMPLEX: bool;
    fn to_complex(self) -> Complex32;
    fn apply_window(self, w: f32) -> Self;
}

impl FftSample for Complex32 {
    const IS_COMPLEX: bool = true;

    #[inline]
    fn to_complex(self) -> Complex32 {
        self
    }

    #[inline]
    fn apply_window(self, w: f32) -> Self {
        self * w
    }
}

impl FftSample for f32 {
    const IS_COMPLEX: bool = false;

    #[inline]
    fn to_complex(self) -> Complex32 {
        Complex32::new(self, 0.0)
    }

    #[inline]
    fn apply_window(self, w: f32) -> Self {
        self * w
    }
}

/// Iterator adapter that accumulates samples, computes an FFT,
/// and publishes the magnitude spectrum to a `LatestWriter`.
/// All samples pass through unchanged.
pub struct SpectrumTap<I: Iterator, S: FftSample> {
    inner: I,
    buffer: Vec<Complex32>,
    fft_size: usize,
    window: Vec<f32>,
    writer: LatestWriter<SpectrumSnapshot>,
    planner: FftPlanner<f32>,
    pos: usize,
    sample_rate: f32,
    _phantom: std::marker::PhantomData<S>,
}

impl<I, S> Iterator for SpectrumTap<I, S>
where
    I: Iterator<Item = S>,
    S: FftSample,
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        let sample = self.inner.next()?;

        self.buffer[self.pos] = sample.apply_window(self.window[self.pos]).to_complex();
        self.pos += 1;

        if self.pos >= self.fft_size {
            self.pos = 0;
            self.compute_and_publish();
        }

        Some(sample)
    }
}

impl<I: Iterator, S: FftSample> SpectrumTap<I, S> {
    fn compute_and_publish(&mut self) {
        let fft = self.planner.plan_fft_forward(self.fft_size);
        fft.process(&mut self.buffer);

        let scale = 1.0 / (self.fft_size as f32);

        let magnitudes_db: Vec<f32> = if S::IS_COMPLEX {
            // Complex input: all N bins are meaningful.
            // FFT output is [0, fs/N, ..., fs/2, -fs/2+fs/N, ..., -fs/N].
            // FFT-shift: put DC in the center → [−fs/2 ... 0 ... +fs/2).
            let n = self.fft_size;
            let half = n / 2;
            (0..n)
                .map(|i| {
                    // fftshift index: second half first, then first half
                    let src = (i + half) % n;
                    let mag = self.buffer[src].norm() * scale;
                    20.0 * mag.max(1e-10).log10()
                })
                .collect()
        } else {
            // Real input: only first N/2 bins are unique (0 to fs/2).
            let half = self.fft_size / 2;
            self.buffer[..half]
                .iter()
                .map(|c| {
                    let mag = c.norm() * scale;
                    20.0 * mag.max(1e-10).log10()
                })
                .collect()
        };

        self.writer.publish(SpectrumSnapshot {
            magnitudes_db,
            sample_rate: self.sample_rate,
            is_complex: S::IS_COMPLEX,
        });
    }
}

fn hann_window(size: usize) -> Vec<f32> {
    use std::f32::consts::PI;
    (0..size)
        .map(|i| {
            let x = (2.0 * PI * i as f32) / (size as f32 - 1.0);
            0.5 * (1.0 - x.cos())
        })
        .collect()
}

pub trait SpectrumTappable: Iterator + Sized {
    fn spectrum_tap<S: FftSample>(
        self,
        fft_size: usize,
        sample_rate: f32,
        writer: LatestWriter<SpectrumSnapshot>,
    ) -> SpectrumTap<Self, S>
    where
        Self: Iterator<Item = S>,
    {
        SpectrumTap {
            inner: self,
            buffer: vec![Complex32::new(0.0, 0.0); fft_size],
            fft_size,
            window: hann_window(fft_size),
            writer,
            planner: FftPlanner::new(),
            pos: 0,
            sample_rate,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<I: Iterator + Sized> SpectrumTappable for I {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::telemetry::latest_value;

    #[test]
    fn test_spectrum_tap_passthrough() {
        let input: Vec<f32> = (0..128).map(|i| i as f32).collect();
        let (tx, _rx) = latest_value::<SpectrumSnapshot>();

        let output: Vec<f32> = input
            .clone()
            .into_iter()
            .spectrum_tap(64, 48000.0, tx)
            .collect();

        assert_eq!(input, output, "tap must not alter the sample stream");
    }

    #[test]
    fn test_spectrum_tap_publishes_after_fft_size() {
        let fft_size = 64;
        let (tx, rx) = latest_value::<SpectrumSnapshot>();

        // Feed exactly one FFT block (real input)
        let input: Vec<f32> = vec![0.0; fft_size];
        let _: Vec<f32> = input
            .into_iter()
            .spectrum_tap(fft_size, 48000.0, tx)
            .collect();

        let snap = rx.take_latest().expect("should have published a snapshot");
        assert_eq!(snap.magnitudes_db.len(), fft_size / 2);
        assert_eq!(snap.sample_rate, 48000.0);
        assert!(!snap.is_complex);
    }

    #[test]
    fn test_spectrum_tap_no_publish_before_full() {
        let fft_size = 64;
        let (tx, rx) = latest_value::<SpectrumSnapshot>();

        // Feed less than one FFT block
        let input: Vec<f32> = vec![0.0; fft_size - 1];
        let _: Vec<f32> = input
            .into_iter()
            .spectrum_tap(fft_size, 48000.0, tx)
            .collect();

        assert!(
            rx.take_latest().is_none(),
            "should not publish until fft_size samples accumulated"
        );
    }

    #[test]
    fn test_spectrum_tap_complex_input() {
        let fft_size = 64;
        let (tx, rx) = latest_value::<SpectrumSnapshot>();

        let input: Vec<Complex32> = (0..fft_size * 2)
            .map(|i| Complex32::new((i as f32 * 0.1).sin(), (i as f32 * 0.1).cos()))
            .collect();

        let output: Vec<Complex32> = input
            .clone()
            .into_iter()
            .spectrum_tap(fft_size, 480000.0, tx)
            .collect();

        assert_eq!(input, output);
        let snap = rx.take_latest().expect("should have a snapshot");
        // Complex: all N bins (fftshift'd), not N/2
        assert_eq!(snap.magnitudes_db.len(), fft_size);
        assert!(snap.is_complex);
    }
}
