use crate::filterable::Filter;
use num_complex::Complex32;

pub struct PllOutput {
    pub out: f32,
    pub lock: f32,
}

pub struct ComplexPllOutput {
    pub out: Complex32,
    pub lock: f32,
}

pub struct RealPll<F> where F: Filter<f32> {
    phase: f32,
    increment: f32,
    loop_gain: f32,
    loop_filter: F,
    divisor: f32,
    lock_alpha: f32,
    lock_level: f32,
}

impl <F> RealPll<F> where F: Filter<f32> {
    pub fn new(freq: f32, sample_rate: f32, loop_gain: f32, loop_filter: F, divisor: f32) -> Self {
        // Single-pole IIR coefficient for ~10 Hz lock detector bandwidth
        let lock_alpha = 1.0 - (-2.0 * std::f32::consts::PI * 10.0 / sample_rate).exp();
        RealPll {
            phase: 0.0,
            increment: freq / sample_rate,
            loop_gain,
            loop_filter,
            divisor,
            lock_alpha,
            lock_level: 0.0,
        }
    }

    pub fn process(&mut self, input: f32) -> PllOutput {
        let divided_phase = 2.0 * std::f32::consts::PI * self.phase / self.divisor;
        let feedback = divided_phase.cos();
        let error = self.loop_filter.process(input * feedback) * self.loop_gain;

        // Quadrature lock detector: when locked, input correlates with sin(phase)
        let lock_inst = input * divided_phase.sin();
        self.lock_level += self.lock_alpha * (lock_inst - self.lock_level);

        // Full-rate output at divisor× frequency
        let out = (2.0 * std::f32::consts::PI * self.phase).cos();

        self.phase = (self.phase + self.increment + error) % self.divisor;
        PllOutput { out, lock: self.lock_level }
    }

    /// Same as `process` but outputs a complex sinusoid (cos + j*sin)
    /// at the full-rate frequency.
    pub fn process_complex(&mut self, input: f32) -> ComplexPllOutput {
        let divided_phase = 2.0 * std::f32::consts::PI * self.phase / self.divisor;
        let feedback = divided_phase.cos();
        let error = self.loop_filter.process(input * feedback) * self.loop_gain;

        let lock_inst = input * divided_phase.sin();
        self.lock_level += self.lock_alpha * (lock_inst - self.lock_level);

        let full_phase = 2.0 * std::f32::consts::PI * self.phase;
        let out = Complex32::new(full_phase.cos(), -full_phase.sin());

        self.phase = (self.phase + self.increment + error) % self.divisor;
        ComplexPllOutput { out, lock: self.lock_level }
    }
}
