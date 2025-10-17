use crate::filterable::Filter;

pub struct RealPll<F> where F: Filter<f32> {
    phase: f32,
    increment: f32,
    loop_gain: f32,
    loop_filter: F,
}

impl <F> RealPll<F> where F: Filter<f32> {
    pub fn new(freq: f32, sample_rate: f32, loop_gain: f32, loop_filter: F) -> Self {
        RealPll {
            phase: 0.0,
            increment: freq / sample_rate,
            loop_gain,
            loop_filter,
        }
    }

    pub fn process(&mut self, input: f32) -> f32 {
        let out = (2.0 * std::f32::consts::PI * self.phase).cos();
        let error = self.loop_filter.process(input * out) * self.loop_gain;
        self.phase = (self.phase + self.increment + error) % 1.0;
        out
    }
}