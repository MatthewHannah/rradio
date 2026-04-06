use crate::filterable::{Filterable, Filter};

// Referencing biquad design & notes from https://arachnoid.com/BiQuadDesigner/index.html

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

    pub fn lowpass(fs: f32, cutoff: f32, q: f32) -> Self {
        let w0 = 2.0 * std::f32::consts::PI * cutoff / fs;
        let alpha = w0.sin() / (2.0 * q);
        let b0 = (1.0 - w0.cos()) / 2.0;
        let b1 = 1.0 - w0.cos();
        let b2 = (1.0 - w0.cos()) / 2.0;
        let a0 = 1.0 + alpha;
        let a1 = -2.0 * w0.cos();
        let a2 = 1.0 - alpha;

        Biquad::new(b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0)
    }

    pub fn highpass(fs: f32, cutoff: f32, q: f32) -> Self {
        let w0 = 2.0 * std::f32::consts::PI * cutoff / fs;
        let alpha = w0.sin() / (2.0 * q);
        let b0 = (1.0 + w0.cos()) / 2.0;
        let b1 = -(1.0 + w0.cos());
        let b2 = (1.0 + w0.cos()) / 2.0;
        let a0 = 1.0 + alpha;
        let a1 = -2.0 * w0.cos();
        let a2 = 1.0 - alpha;

        Biquad::new(b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0)
    }
}

impl<Num> Filter<Num> for Biquad<Num> where Num: Filterable<Num> {
    fn process(&mut self, x: Num) -> Num {
        self.process(x)
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

    #[test]
    fn test_lowpass_coefficients() {
        // 300kHz cutoff, 6MHz sample rate, Butterworth Q
        let bq: Biquad<f32> = Biquad::lowpass(6e6, 300000.0, 0.707);
        assert!((bq.b0 - 0.02008282).abs() < 1e-4, "b0 mismatch: {}", bq.b0);
        assert!((bq.b1 - 0.04016564).abs() < 1e-4, "b1 mismatch: {}", bq.b1);
        assert!((bq.b2 - 0.02008282).abs() < 1e-4, "b2 mismatch: {}", bq.b2);
        assert!((bq.a1 - -1.56097580).abs() < 1e-4, "a1 mismatch: {}", bq.a1);
        assert!((bq.a2 - 0.64130708).abs() < 1e-4, "a2 mismatch: {}", bq.a2);
    }

    #[test]
    fn test_highpass_coefficients() {
        // 23kHz cutoff, 240kHz sample rate, Butterworth Q
        let bq: Biquad<f32> = Biquad::highpass(240000.0, 23000.0, 0.707);
        assert!((bq.b0 - 0.65120842).abs() < 1e-4, "b0 mismatch: {}", bq.b0);
        assert!((bq.b1 - -1.30241684).abs() < 1e-4, "b1 mismatch: {}", bq.b1);
        assert!((bq.b2 - 0.65120842).abs() < 1e-4, "b2 mismatch: {}", bq.b2);
        assert!((bq.a1 - -1.17684383).abs() < 1e-4, "a1 mismatch: {}", bq.a1);
        assert!((bq.a2 - 0.42798985).abs() < 1e-4, "a2 mismatch: {}", bq.a2);
    }
}