use crate::filterable::{Filterable,Filter};

#[derive(Debug, Clone)]
pub struct Deemphasis<Num> {
    y1: Num,
    x1: Num,
    b0: f32,
    b1: f32,
    a1: f32
}

impl<Num> Deemphasis<Num> where Num: Filterable<Num> {
    pub fn new(fs: f32, tau: f32) -> Self {
        let a = ((tau * fs).recip() * 0.5).tan().recip();
        Deemphasis {
            y1: Num::zero(),
            x1: Num::zero(),
            b0: (1.0 + a).recip(),
            b1: (1.0 + a).recip(),
            a1: (1.0 - a) * (1.0 + a).recip()
        }
    }
}

impl<Num> Filter<Num> for Deemphasis<Num> where Num: Filterable<Num> {
    fn process(&mut self, x: Num) -> Num {
        let y = (x * self.b0 + self.x1 * self.b1) - self.y1 * self.a1;
        self.y1 = y;
        self.x1 = x;
        y
    }
}
