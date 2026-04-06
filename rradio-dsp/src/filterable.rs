use std::ops::{Mul, Add, Sub};
use num_traits::Zero;

pub trait Filterable<Num> :
    Add<Num, Output=Num> + Sub<Num, Output=Num> + Mul<f32,Output = Num> + Zero + Copy {}

impl<T> Filterable<T> for T where
    T: Add<T, Output=T> + Sub<T, Output=T> + Mul<f32,Output = T> + Zero + Copy {}

pub trait Filter<Num> where Num: Filterable<Num> {
    fn process(&mut self, x: Num) -> Num;
}

pub struct FilterIter<I, F, Num> where I: Iterator<Item = Num> + Sized, F: Filter<Num> + Sized, Num: Filterable<Num> {
    iter: I,
    filter: F,
}

impl<I, F, Num> FilterIter<I, F, Num> where I: Iterator<Item = Num>, F: Filter<Num>, Num: Filterable<Num> {
    pub fn new(iter: I, f: F) -> Self {
        FilterIter {
            iter,
            filter: f,
        }
    }
}

impl<I, F, Num> Iterator for FilterIter<I, F, Num> where I: Iterator<Item = Num>, F: Filter<Num>, Num: Filterable<Num> {
    type Item = Num;

    fn next(&mut self) -> Option<Self::Item> {
        let sample = self.iter.next()?;
        Some(self.filter.process(sample))
    }
}

pub trait FilterableIter<Num> where Num: Filterable<Num> {
    fn dsp_filter<F>(self, f: F) -> FilterIter<Self, F, Self::Item> where Self: Sized, Self: Iterator<Item = Num>, F: Filter<Num>;
}

impl<I, Num> FilterableIter<Num> for I where I: Iterator<Item = Num> + Sized, Num: Filterable<Num> {
    fn dsp_filter<F>(self, f: F) -> FilterIter<Self, F, Num> where F: Filter<Num> {
        FilterIter::new(self, f)
    }
}
