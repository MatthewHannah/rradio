use crate::filterable::{Filterable,Filter};

#[derive(Debug, Clone)]
pub struct Deemphasis<Num> {
    y1: Num,
    k: f32
}

impl<Num> Deemphasis<Num> where Num: Filterable<Num> {
    pub fn new(fs: f32, tau: f32) -> Self {
        Deemphasis {
            y1: Num::zero(),
            k: fs * tau,
        }
    }
}

impl<Num> Filter<Num> for Deemphasis<Num> where Num: Filterable<Num> {
    fn process(&mut self, x: Num) -> Num {
        let y = (x + self.y1 * self.k) * (1.0 + self.k).recip();
        self.y1 = y;
        y
    }
}
