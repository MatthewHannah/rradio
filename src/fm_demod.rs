use crate::biquad;

use num_complex::Complex32;

pub struct FmDemodulator {
    loop_filter: biquad::Biquad<f32>,
    last_phase: f32,
}

impl FmDemodulator {
    pub fn new(loop_filter: biquad::Biquad<f32>) -> Self {
        FmDemodulator {
            loop_filter: loop_filter,
            last_phase: 0.0,
        }
    }

    pub fn process(&mut self, sample: num_complex::Complex32) -> f32 {
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

pub struct FmDemodIter<I> {
    iter: I,
    demodulator: FmDemodulator,
}

impl<I> FmDemodIter<I> where I: Iterator<Item = Complex32> {
    pub fn new(iter: I, loop_filter: biquad::Biquad<f32>) -> Self {
        FmDemodIter {
            iter: iter,
            demodulator: FmDemodulator::new(loop_filter),
        }
    }
}

impl<I> Iterator for FmDemodIter<I> where I: Iterator<Item = Complex32> {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        let sample = self.iter.next()?;
        Some(self.demodulator.process(sample))
    }
}

pub trait FmDemodulatable {
    fn fm_demodulate(self, loop_filter: biquad::Biquad<f32>) -> FmDemodIter<Self> where Self: Sized;
}

impl<I> FmDemodulatable for I where I: Iterator<Item = Complex32> {
    fn fm_demodulate(self, loop_filter: biquad::Biquad<f32>) -> FmDemodIter<Self> {
        FmDemodIter::new(self, loop_filter)
    }
}