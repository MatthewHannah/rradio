/// Developed using Claude Opus 4.6

use crate::filterable::{Filterable, Filter};

#[derive(Debug, Clone)]
pub struct Fir<Num> {
    coeffs: Vec<f32>,
    delay_line: Vec<Num>,
    head: usize,
    len: usize,
}

impl<Num: Filterable<Num>> Fir<Num> {
    pub fn new(coeffs: Vec<f32>) -> Self {
        let len = coeffs.len();
        Fir {
            coeffs,
            delay_line: vec![Num::zero(); len * 2], // doubled for contiguous access
            head: 0,
            len,
        }
    }

    /// Push a sample into the delay line without computing output.
    #[inline]
    pub fn push(&mut self, x: Num) {
        self.delay_line[self.head] = x;
        self.delay_line[self.head + self.len] = x; // duplicate for wraparound
        self.head = (self.head + 1) % self.len;
    }

    /// Compute the filter output from the current delay line state.
    /// head points to the oldest sample; head+len-1 is the newest.
    #[inline]
    pub fn execute(&self) -> Num {
        self.delay_line[self.head..self.head + self.len]
            .iter()
            .zip(self.coeffs.iter().rev())
            .fold(Num::zero(), |acc, (&s, &c)| acc + s * c)
    }

    pub fn process(&mut self, x: Num) -> Num {
        self.push(x);
        self.execute()
    }
}

impl<Num> Filter<Num> for Fir<Num> where Num: Filterable<Num> {
    fn process(&mut self, x: Num) -> Num {
        self.process(x)
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
