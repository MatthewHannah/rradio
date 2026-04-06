use std::ops::Mul;

use num_complex::Complex32;

#[derive(Clone)]
pub struct Osc {
    phase: f32,
    increment: f32,
}

impl Osc {
    pub fn new(freq: f32, sample_rate: f32) -> Self {
        Osc {
            phase: 0.0,
            increment: freq / sample_rate,
        }
    }

    pub fn next(&mut self) -> Complex32 {
        let phase = self.phase * 2.0 * std::f32::consts::PI;
        self.phase = (self.phase + self.increment) % 1.0;
        Complex32::new(phase.cos(), phase.sin())
    }
}

impl Iterator for Osc {
    type Item = Complex32;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.next())
    }
}

pub struct MixableIter<T, I: Iterator<Item = T>> where T: Mul<Complex32, Output = Complex32> {
    iter: I,
    osc: Osc,
}

impl<T, I: Iterator<Item = T>> MixableIter<T, I> where T: Mul<Complex32, Output = Complex32> {
    pub fn new(iter: I, osc: Osc) -> Self {
        MixableIter { iter, osc }
    }

    pub fn next(&mut self) -> Option<Complex32> {
        let sample = self.iter.next()?;
        let osc_sample = self.osc.next();
        Some(sample * osc_sample)
    }
}

impl<T, I: Iterator<Item = T>> Iterator for MixableIter<T, I> where T: Mul<Complex32, Output = Complex32> {
    type Item = Complex32;

    fn next(&mut self) -> Option<Self::Item> {
        self.next()
    }
}

pub trait Mixable<T>: Iterator<Item = T> + Sized where T: Mul<Complex32, Output = Complex32> {
    fn mix(self, freq: f32, sample_rate: f32) -> MixableIter<T, Self>;
}

impl<T, I: Iterator<Item = T>> Mixable<T> for I where T: Mul<Complex32, Output = Complex32> {
    fn mix(self, freq: f32, sample_rate: f32) -> MixableIter<T, Self> {
        MixableIter::new(self, Osc::new(freq, sample_rate))
    }
}


#[cfg(test)]
mod tests {
    use num_complex::ComplexFloat;

    use super::*;

    #[test]
    fn test_osc() {
        let o = Osc::new(1e3, 1e6);
        let samples: Vec<Complex32> = o.take(2000).collect();
        assert!(samples[0].re() - 1.0 < 0.01);
        assert!(samples[1000].re() - 1.0 < 0.01);
        assert!(samples[250].im() - 1.0 < 0.01);
    }
}