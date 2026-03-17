/// Developed using Claude Opus 4.6

use crate::filterable::{Filterable, Filter};

#[derive(Debug, Clone)]
pub struct Fir<Num> {
    coeffs: Vec<f32>,
    delay_line: Vec<Num>,
    head: usize,
}

impl<Num: Filterable<Num>> Fir<Num> {
    pub fn new(coeffs: Vec<f32>) -> Self {
        let len = coeffs.len();
        Fir {
            coeffs,
            delay_line: vec![Num::zero(); len],
            head: 0,
        }
    }
}

impl<Num> Filter<Num> for Fir<Num> where Num: Filterable<Num> {
    fn process(&mut self, x: Num) -> Num {
        let len = self.coeffs.len();
        self.delay_line[self.head] = x;
        self.head = (self.head + 1) % len;

        let mut y = Num::zero();
        for i in 0..len {
            // Read backwards from head: most recent sample pairs with coeffs[0]
            let idx = (self.head + len - 1 - i) % len;
            y = y + self.delay_line[idx] * self.coeffs[i];
        }
        y
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fir() {
        let coeffs = vec![0.25, 0.25, 0.25, 0.25];
        let mut fir = Fir::new(coeffs);

        let input = vec![1.0, 2.0, 3.0, 4.0];
        let expected_output = vec![0.25, 0.75, 1.5, 2.5];

        for (x, expected) in input.into_iter().zip(expected_output.into_iter()) {
            let y = fir.process(x);
            assert!((y - expected).abs() < 1e-6);
        }
    }

    #[test]
    fn test_fir_impulse_response() {
        let coeffs = vec![0.25, 0.25, 0.25, 0.25];
        let mut fir = Fir::new(coeffs);

        let input = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let expected_output = vec![0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.0];

        for (x, expected) in input.into_iter().zip(expected_output.into_iter()) {
            let y = fir.process(x);
            assert!((y - expected).abs() < 1e-6);
        }
    }

    #[test]
    fn test_fir_no_op() {
        let coeffs = vec![1.0];
        let mut fir = Fir::new(coeffs);

        let input = vec![1.0, 2.0, 3.0, 4.0];
        let expected_output = vec![1.0, 2.0, 3.0, 4.0];

        for (x, expected) in input.into_iter().zip(expected_output.into_iter()) {
            let y = fir.process(x);
            assert!((y - expected).abs() < 1e-6);
        }
    }
}
