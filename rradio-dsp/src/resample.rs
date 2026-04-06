use crate::filterable::Filterable;

use num_traits::Zero;

pub struct Upsample<I, T> where I: Iterator<Item = T>, T: Copy + Zero {
    iter: I,
    factor: usize,
    idx: usize,
}

impl<I, T> Upsample<I, T> where I: Iterator<Item = T>, T: Copy + Zero {
    pub fn new(iter: I, factor: usize) -> Self {
        Upsample {
            iter: iter,
            factor: factor,
            idx: 0,
        }
    }
}

impl<I, T> Iterator for Upsample<I, T> where I: Iterator<Item = T>, T: Copy + Zero {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        let out = if self.idx == 0 {
            self.iter.next()
        } else {
            Some(T::zero())
        };
        self.idx = (self.idx + 1) % self.factor;
        out
    }
}

#[allow(dead_code)]
pub trait Upsampleable<T> {
    fn upsample(self, factor: usize) -> Upsample<Self, T> where Self: Sized + Iterator<Item = T>, T: Copy + Zero;
}

impl<I, T> Upsampleable<T> for I where I: Iterator<Item = T>, T: Copy + Zero {
    fn upsample(self, factor: usize) -> Upsample<I,T> {
        Upsample::new(self, factor)
    }
}

pub struct Downsample<I> {
    iter: I,
    factor: usize,
}

impl<I> Downsample<I> {
    pub fn new(iter: I, factor: usize) -> Self {
        Downsample {
            iter: iter,
            factor: factor,
        }
    }
}

impl<I> Iterator for Downsample<I> where I: Iterator {
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        let val = self.iter.next();
        for _ in 0..self.factor-1 {
            self.iter.next();
        }
        val
    }
}

pub trait Downsampleable {
    fn downsample(self, factor: usize) -> Downsample<Self> where Self: Sized;
}

impl<I> Downsampleable for I where I: Iterator {
    fn downsample(self, factor: usize) -> Downsample<I> {
        Downsample::new(self, factor)
    }
}

/// Polyphase rational resampler: combines FIR filtering with L/M resampling.
pub struct RationalResampler<Num: Filterable<Num>> {
    phases: Vec<Vec<f32>>,
    arm_len: usize,
    buffer: Vec<Num>,
    write_pos: usize,
    phase: usize,
    next_advance: usize, // how many input samples to consume before next output
    l: usize,
    m: usize,
}

impl<Num: Filterable<Num>> RationalResampler<Num> {
    /// Create a new resampler.
    ///
    /// * `prototype` - FIR filter taps designed at L × fs_out (= fs_in × L / M × L... no, at fs_in)
    /// * `l` - interpolation factor
    /// * `m` - decimation factor
    ///
    /// Output rate = input_rate × L / M
    pub fn new(prototype: Vec<f32>, l: usize, m: usize) -> Self {
        let arm_len = (prototype.len() + l - 1) / l;

        // Phase k gets taps[k], taps[k+L], taps[k+2L], ...
        // Stored reversed so compute() is a forward dot product over
        // the oldest-to-newest doubled buffer.
        let mut phases = Vec::with_capacity(l);
        for k in 0..l {
            let mut arm: Vec<f32> = (0..arm_len)
                .map(|i| {
                    let idx = k + i * l;
                    if idx < prototype.len() {
                        prototype[idx]
                    } else {
                        0.0
                    }
                })
                .collect();
            arm.reverse();
            phases.push(arm);
        }

        RationalResampler {
            phases,
            arm_len,
            buffer: vec![Num::zero(); arm_len * 2], // doubled for contiguous access
            write_pos: 0,
            phase: 0,
            next_advance: arm_len,
            l,
            m,
        }
    }

    #[inline]
    fn push(&mut self, sample: Num) {
        self.buffer[self.write_pos] = sample;
        self.buffer[self.write_pos + self.arm_len] = sample;
        self.write_pos = (self.write_pos + 1) % self.arm_len;
    }

    #[inline]
    fn compute(&self) -> Num {
        let arm = &self.phases[self.phase];
        let start = self.write_pos; // oldest sample in doubled layout
        let mut sum = Num::zero();
        for i in 0..self.arm_len {
            sum = sum + self.buffer[start + i] * arm[i];
        }
        sum
    }

    /// Process one output sample, pulling inputs as needed.
    /// Returns None when the input is exhausted.
    pub fn process<I: Iterator<Item = Num>>(&mut self, iter: &mut I) -> Option<Num> {
        for _ in 0..self.next_advance {
            self.push(iter.next()?);
        }

        let output = self.compute();

        self.phase += self.m;
        self.next_advance = self.phase / self.l;
        self.phase %= self.l;

        Some(output)
    }
}

// --- Iterator adapter ---

pub struct RationalResampleIter<I, Num>
where
    I: Iterator<Item = Num>,
    Num: Filterable<Num>,
{
    iter: I,
    resampler: RationalResampler<Num>,
}

impl<I, Num> Iterator for RationalResampleIter<I, Num>
where
    I: Iterator<Item = Num>,
    Num: Filterable<Num>,
{
    type Item = Num;

    fn next(&mut self) -> Option<Num> {
        self.resampler.process(&mut self.iter)
    }
}

pub trait RationalResampleable<Num: Filterable<Num>> {
    fn resample(
        self,
        prototype: Vec<f32>, // FIR filter taps designed at fs_in * l
        l: usize, // interpolation factor
        m: usize, // decimation factor
    ) -> RationalResampleIter<Self, Num>
    where
        Self: Sized + Iterator<Item = Num>;
}

impl<I, Num> RationalResampleable<Num> for I
where
    I: Iterator<Item = Num>,
    Num: Filterable<Num>,
{
    fn resample(
        self,
        prototype: Vec<f32>,
        l: usize,
        m: usize,
    ) -> RationalResampleIter<Self, Num> {
        RationalResampleIter {
            iter: self,
            resampler: RationalResampler::new(prototype, l, m),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex32;

    #[test]
    fn test_upsample() {
        let input = vec![Complex32::new(1.0, 0.0), Complex32::new(2.0, 0.0), Complex32::new(3.0, 0.0)];
        let upsampled: Vec<Complex32> = input.into_iter().upsample(3).collect();
        let expected = vec![
            Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0), Complex32::new(0.0, 0.0),
            Complex32::new(2.0, 0.0), Complex32::new(0.0, 0.0), Complex32::new(0.0, 0.0),
            Complex32::new(3.0, 0.0), Complex32::new(0.0, 0.0), Complex32::new(0.0, 0.0),
        ];
        assert_eq!(upsampled, expected);
    }

    #[test]
    fn test_downsample() {
        let input = vec![Complex32::new(1.0, 0.0), Complex32::new(2.0, 0.0), Complex32::new(3.0, 0.0), Complex32::new(4.0, 0.0), Complex32::new(5.0, 0.0)];
        let downsampled: Vec<Complex32> = input.into_iter().downsample(2).collect();
        let expected = vec![Complex32::new(1.0, 0.0), Complex32::new(3.0, 0.0), Complex32::new(5.0, 0.0)];
        assert_eq!(downsampled, expected);
    }

    #[test]
    fn test_rate_ratio() {
        // 4 input samples at rate 4, resample to rate 3 (L=3, M=4)
        // Prototype: simple 6-tap filter
        let taps = vec![1.0, 2.0, 3.0, 3.0, 2.0, 1.0];
        let input: Vec<f32> = (0..100).map(|i| i as f32).collect();
        let output: Vec<f32> = input
            .into_iter()
            .resample(taps, 3, 4)
            .collect();
        // Should produce ~75 output samples (100 * 3/4)
        assert_eq!(output.len(), 75);
    }

    #[test]
    fn test_identity_resample() {
        // L=1, M=1 should pass through (with filter applied)
        let taps = vec![1.0]; // identity filter
        let input: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let output: Vec<f32> = input
            .into_iter()
            .resample(taps, 1, 1)
            .collect();
        assert_eq!(output, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    }

    #[test]
    fn test_downsample_by_2() {
        // L=1, M=2 with identity filter = take every other sample
        let taps = vec![1.0];
        let input: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let output: Vec<f32> = input
            .into_iter()
            .resample(taps, 1, 2)
            .collect();
        assert_eq!(output, vec![1.0, 3.0, 5.0]);
    }
}