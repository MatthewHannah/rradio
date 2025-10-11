use std::ops::{Mul, Add, Sub};
use num_traits::Zero;

pub trait Filterable<Num> :
    Add<Num, Output=Num> + Sub<Num, Output=Num> + Mul<f32,Output = Num> + Zero + Copy {}

impl<T> Filterable<T> for T where
    T: Add<T, Output=T> + Sub<T, Output=T> + Mul<f32,Output = T> + Zero + Copy {}

#[derive(Debug, Clone)]
pub struct Biquad<Num> {
    b0: f32,
    b1: f32,
    b2: f32,
    a1: f32,
    a2: f32,
    x1: Num,
    x2: Num,
    y1: Num,
    y2: Num,
}

impl<Num> Biquad<Num> where Num: Filterable<Num> {
    pub fn new(b0: f32, b1: f32, b2: f32, a1: f32, a2: f32) -> Self {
        Biquad {
            b0,
            b1,
            b2,
            a1,
            a2,
            x1: Num::zero(),
            x2: Num::zero(),
            y1: Num::zero(),
            y2: Num::zero(),
        }
    }

    pub fn process(&mut self, x: Num) -> Num {
        let y = x * self.b0 + self.x1 * self.b1 + self.x2 * self.b2 - self.y1 * self.a1 - self.y2 * self.a2;
        self.x2 = self.x1;
        self.x1 = x;
        self.y2 = self.y1;
        self.y1 = y;
        y
    }
}

pub struct BiquadIter<I, Num> where I: Iterator<Item = Num>, Num: Filterable<Num> {
    iter: I,
    biquad: Biquad<Num>,
}

impl<I, Num> BiquadIter<I, Num> where I: Iterator<Item = Num>, Num: Filterable<Num> {
    pub fn new(iter: I, b: Biquad<Num>) -> Self {
        BiquadIter {
            iter,
            biquad: b,
        }
    }
}

impl<I, Num> Iterator for BiquadIter<I, Num> where I: Iterator<Item = Num>, Num: Filterable<Num> {
    type Item = Num;

    fn next(&mut self) -> Option<Self::Item> {
        let sample = self.iter.next()?;
        Some(self.biquad.process(sample))
    }
}

pub trait Biquadable<Num> where Num: Filterable<Num> {
    fn biquad(self, b: Biquad<Num>) -> BiquadIter<Self, Self::Item> where Self: Sized, Self: Iterator<Item = Num>;
}

impl<I, Num> Biquadable<Num> for I where I: Iterator<Item = Num>, Num: Filterable<Num> {
    fn biquad(self, b: Biquad<Num>) -> BiquadIter<Self, Num> {
        BiquadIter::new(self, b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex32;

    #[test]
    fn test_biquad_complex() {
        let mut bq = Biquad::new(0.02008282, 0.04016564, 0.02008282, -1.56097580, 0.64130708);
        let input = vec![1.0, 0.0, 0.0, 0.0, 0.0].iter().map(|&x| Complex32::new(x, 0.0)).collect::<Vec<Complex32>>();
        let output: Vec<Complex32> = input.iter().map(|&x| bq.process(x)).collect();
        let expected =  vec![0.02008282, 0.07151443602, 0.1188358693, 0.139637202, 0.1417600088].iter().map(|&x| Complex32::new(x, 0.0)).collect::<Vec<Complex32>>();
        for (o, e) in output.iter().zip(expected.iter()) {
            assert!((o - e).norm() < 1e-4);
        }
    }

    #[test]
    fn test_biquad_f32() {
        let mut bq = Biquad::new(0.02008282, 0.04016564, 0.02008282, -1.56097580, 0.64130708);
        let input = vec![1.0, 0.0, 0.0, 0.0, 0.0];
        let output: Vec<f32> = input.iter().map(|&x| bq.process(x)).collect();
        let expected = vec![0.02008282, 0.07151443602, 0.1188358693, 0.139637202, 0.1417600088];
        for (o, e) in output.iter().zip(expected.iter()) {
            assert!((o - e).abs() < 1e-4);
        }
    }
}