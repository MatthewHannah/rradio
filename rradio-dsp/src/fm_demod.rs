use num_complex::Complex32;

pub struct FmDemodulator {
    last_phase: f32,
}

impl FmDemodulator {
    pub fn new() -> Self {
        FmDemodulator {
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
        delta
    }
}

pub struct FmDemodIter<I> {
    iter: I,
    demodulator: FmDemodulator,
}

impl<I> FmDemodIter<I> where I: Iterator<Item = Complex32> {
    pub fn new(iter: I) -> Self {
        FmDemodIter {
            iter,
            demodulator: FmDemodulator::new(),
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
    fn fm_demodulate(self) -> FmDemodIter<Self> where Self: Sized;
}

impl<I> FmDemodulatable for I where I: Iterator<Item = Complex32> {
    fn fm_demodulate(self) -> FmDemodIter<Self> {
        FmDemodIter::new(self)
    }
}
