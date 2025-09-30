use num_complex::Complex32;

pub struct Biquad {
    b0: f32,
    b1: f32,
    b2: f32,
    a1: f32,
    a2: f32,
    x1: Complex32,
    x2: Complex32,
    y1: Complex32,
    y2: Complex32,
}

impl Biquad {
    pub fn new(b0: f32, b1: f32, b2: f32, a1: f32, a2: f32) -> Self {
        Biquad {
            b0,
            b1,
            b2,
            a1,
            a2,
            x1: Complex32::new(0.0, 0.0),
            x2: Complex32::new(0.0, 0.0),
            y1: Complex32::new(0.0, 0.0),
            y2: Complex32::new(0.0, 0.0),
        }
    }

    pub fn process(&mut self, x: Complex32) -> Complex32 {
        let y = self.b0 * x + self.b1 * self.x1 + self.b2 * self.x2 - self.a1 * self.y1 - self.a2 * self.y2;
        self.x2 = self.x1;
        self.x1 = x;
        self.y2 = self.y1;
        self.y1 = y;
        y
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_biquad() {
        let mut bq = Biquad::new(0.02008282, 0.04016564, 0.02008282, -1.56097580, 0.64130708);
        let input = vec![1.0, 0.0, 0.0, 0.0, 0.0].iter().map(|&x| Complex32::new(x, 0.0)).collect::<Vec<Complex32>>();
        let output: Vec<Complex32> = input.iter().map(|&x| bq.process(x)).collect();
        let expected =  vec![0.02008282, 0.07151443602, 0.1188358693, 0.139637202, 0.1417600088].iter().map(|&x| Complex32::new(x, 0.0)).collect::<Vec<Complex32>>();
        for (o, e) in output.iter().zip(expected.iter()) {
            assert!((o - e).norm() < 1e-4);
        }
    }
}