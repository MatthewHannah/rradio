pub struct PiLoopFilter {
    pub integrator: f32,
    alpha: f32,
    beta: f32,

    min_integrator: f32,
    max_integrator: f32,
}

fn calculate_gains(loop_bw: f32, damping: f32, k: f32) -> (f32, f32) {
    let loop_bw_times_damping = loop_bw * damping;

    let k0 = 2.0 / k;
    let k1 = (-loop_bw_times_damping).exp();
    let sinh_loopbw_damping = loop_bw_times_damping.sinh();

    let cos_adjustment = if damping > 1.0 {
        (loop_bw * (damping * damping - 1.0).sqrt()).cosh()
    } else if (damping - 1.0).abs() < 1e-6 {
        1.0
    } else {
        (loop_bw * (1.0 - damping * damping).sqrt()).cos()
    };

    let alpha = k0 * k1 * sinh_loopbw_damping;
    let beta = k0 * (1.0 - k1 * (sinh_loopbw_damping + cos_adjustment));

    (alpha, beta)
}

impl PiLoopFilter {
    pub fn new(loop_bw: f32, damping: f32, k: f32, min_integrator: f32, max_integrator: f32) -> Self {
        let (alpha, beta) = calculate_gains(loop_bw, damping, k);
        PiLoopFilter {
            integrator: (min_integrator + max_integrator) / 2.0, // start in middle of allowed range
            alpha,
            beta,
            min_integrator,
            max_integrator,
        }
    }

    pub fn advance(&mut self, error: f32) -> f32 {
        self.integrator += self.beta * error;
        self.integrator = self.integrator.clamp(self.min_integrator, self.max_integrator);

        let output = self.integrator + self.alpha * error;
        output
    }
}