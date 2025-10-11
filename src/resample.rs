use num_traits::Zero;

pub struct Upsample<I, T> where I: Iterator<Item = T>, T: Copy + Zero {
    iter: I,
    factor: usize,
    cur: T,
    idx: usize,
}

impl<I, T> Upsample<I, T> where I: Iterator<Item = T>, T: Copy + Zero {
    pub fn new(iter: I, factor: usize) -> Self {
        Upsample {
            iter: iter,
            factor: factor,
            cur: T::zero(),
            idx: 0,
        }
    }
}

impl<I, T> Iterator for Upsample<I, T> where I: Iterator<Item = T>, T: Copy + Zero {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx == 0 {
            self.cur = self.iter.next()?;
        }
        self.idx = (self.idx + 1) % self.factor;
        Some(self.cur)
    }
}

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


#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex32;

    #[test]
    fn test_upsample() {
        let input = vec![Complex32::new(1.0, 0.0), Complex32::new(2.0, 0.0), Complex32::new(3.0, 0.0)];
        let upsampled: Vec<Complex32> = input.into_iter().upsample(3).collect();
        let expected = vec![
            Complex32::new(1.0, 0.0), Complex32::new(1.0, 0.0), Complex32::new(1.0, 0.0),
            Complex32::new(2.0, 0.0), Complex32::new(2.0, 0.0), Complex32::new(2.0, 0.0),
            Complex32::new(3.0, 0.0), Complex32::new(3.0, 0.0), Complex32::new(3.0, 0.0),
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
}