use crate::filterable::Filter;

use num_complex::Complex32;

pub struct FmDemodulator<F: Filter<f32>> {
    loop_filter: F,
    last_phase: f32,
}

impl<F: Filter<f32>> FmDemodulator<F> {
    pub fn new(loop_filter: F) -> Self {
        FmDemodulator {
            loop_filter,
            last_phase: 0.0,
        }
    }

    pub fn process(&mut self, sample: Complex32) -> f32 {
        let phase = sample.arg();
        let mut delta = phase - self.last_phase;
        if delta > std::f32::consts::PI {
            delta -= 2.0 * std::f32::consts::PI;
        } else if delta < -std::f32::consts::PI {
            delta += 2.0 * std::f32::consts::PI;
        }
        self.last_phase = phase;
        self.loop_filter.process(delta)
    }
}

pub struct FmDemodIter<I, F: Filter<f32>> {
    iter: I,
    demodulator: FmDemodulator<F>,
}

impl<I, F: Filter<f32>> FmDemodIter<I, F> where I: Iterator<Item = Complex32> {
    pub fn new(iter: I, loop_filter: F) -> Self {
        FmDemodIter {
            iter,
            demodulator: FmDemodulator::new(loop_filter),
        }
    }
}

impl<I, F: Filter<f32>> Iterator for FmDemodIter<I, F> where I: Iterator<Item = Complex32> {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        let sample = self.iter.next()?;
        Some(self.demodulator.process(sample))
    }
}

pub trait FmDemodulatable {
    fn fm_demodulate<F: Filter<f32>>(self, loop_filter: F) -> FmDemodIter<Self, F> where Self: Sized;
}

impl<I> FmDemodulatable for I where I: Iterator<Item = Complex32> {
    fn fm_demodulate<F: Filter<f32>>(self, loop_filter: F) -> FmDemodIter<Self, F> {
        FmDemodIter::new(self, loop_filter)
    }
}
