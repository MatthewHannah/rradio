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

    pub fn iter_mut(&mut self) -> &mut Self {
        self
    }
}

impl Iterator for Osc {
    type Item = Complex32;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.next())
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