use crate::filterable::Filterable;

/// Polyphase rational resampler: combines FIR filtering with L/M resampling.
///
/// Decomposes a prototype FIR filter into L polyphase arms. For each output
/// sample, evaluates exactly one arm (ceil(N/L) taps), giving the same result
/// as upsample-by-L → filter → downsample-by-M but at a fraction of the cost.
pub struct PolyphaseResampler<Num: Filterable<Num>> {
    phases: Vec<Vec<f32>>,
    arm_len: usize,
    buffer: Vec<Num>,
    write_pos: usize,
    phase: usize,
    next_advance: usize,
    l: usize,
    m: usize,
}

impl<Num: Filterable<Num>> PolyphaseResampler<Num> {
    /// Create a new polyphase resampler.
    ///
    /// * `prototype` - FIR filter taps designed at L × fs_out (= fs_in × L / M × L... no, at fs_in)
    /// * `l` - interpolation factor
    /// * `m` - decimation factor
    ///
    /// Output rate = input_rate × L / M
    pub fn new(prototype: Vec<f32>, l: usize, m: usize) -> Self {
        let arm_len = (prototype.len() + l - 1) / l;

        // Phase k gets taps[k], taps[k+L], taps[k+2L], ...
        let mut phases = Vec::with_capacity(l);
        for k in 0..l {
            let arm: Vec<f32> = (0..arm_len)
                .map(|i| {
                    let idx = k + i * l;
                    if idx < prototype.len() {
                        prototype[idx]
                    } else {
                        0.0
                    }
                })
                .collect();
            phases.push(arm);
        }

        PolyphaseResampler {
            phases,
            arm_len,
            buffer: vec![Num::zero(); arm_len],
            write_pos: 0,
            phase: 0,
            next_advance: arm_len, // prime the buffer on first output
            l,
            m,
        }
    }

    #[inline]
    fn push(&mut self, sample: Num) {
        self.buffer[self.write_pos] = sample;
        self.write_pos = (self.write_pos + 1) % self.arm_len;
    }

    #[inline]
    fn compute(&self) -> Num {
        let arm = &self.phases[self.phase];
        let mut sum = Num::zero();
        for i in 0..self.arm_len {
            let idx = (self.write_pos + self.arm_len - 1 - i) % self.arm_len;
            sum = sum + self.buffer[idx] * arm[i];
        }
        sum
    }

    /// Process one output sample, pulling inputs as needed.
    /// Returns None when the input is exhausted.
    #[inline]
    pub fn next_from<I: Iterator<Item = Num>>(&mut self, iter: &mut I) -> Option<Num> {
        for _ in 0..self.next_advance {
            self.push(iter.next()?);
        }

        let output = self.compute();

        self.phase += self.m;
        self.next_advance = self.phase / self.l;
        self.phase %= self.l;

        Some(output)
    }

    pub fn arm_len(&self) -> usize {
        self.arm_len
    }

    pub fn num_phases(&self) -> usize {
        self.l
    }
}

// --- Iterator adapter ---

pub struct PolyphaseResampleIter<I, Num>
where
    I: Iterator<Item = Num>,
    Num: Filterable<Num>,
{
    iter: I,
    resampler: PolyphaseResampler<Num>,
}

impl<I, Num> Iterator for PolyphaseResampleIter<I, Num>
where
    I: Iterator<Item = Num>,
    Num: Filterable<Num>,
{
    type Item = Num;

    #[inline]
    fn next(&mut self) -> Option<Num> {
        self.resampler.next_from(&mut self.iter)
    }
}

pub trait PolyphaseResampleable<Num: Filterable<Num>> {
    fn polyphase_resample(
        self,
        prototype: Vec<f32>,
        l: usize,
        m: usize,
    ) -> PolyphaseResampleIter<Self, Num>
    where
        Self: Sized + Iterator<Item = Num>;
}

impl<I, Num> PolyphaseResampleable<Num> for I
where
    I: Iterator<Item = Num>,
    Num: Filterable<Num>,
{
    fn polyphase_resample(
        self,
        prototype: Vec<f32>,
        l: usize,
        m: usize,
    ) -> PolyphaseResampleIter<Self, Num> {
        PolyphaseResampleIter {
            iter: self,
            resampler: PolyphaseResampler::new(prototype, l, m),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rate_ratio() {
        // 4 input samples at rate 4, resample to rate 3 (L=3, M=4)
        // Prototype: simple 6-tap filter
        let taps = vec![1.0, 2.0, 3.0, 3.0, 2.0, 1.0];
        let input: Vec<f32> = (0..100).map(|i| i as f32).collect();
        let output: Vec<f32> = input
            .into_iter()
            .polyphase_resample(taps, 3, 4)
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
            .polyphase_resample(taps, 1, 1)
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
            .polyphase_resample(taps, 1, 2)
            .collect();
        assert_eq!(output, vec![1.0, 3.0, 5.0]);
    }
}
