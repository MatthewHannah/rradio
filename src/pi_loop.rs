pub struct PiLoopFilter {
    pub integrator: f32,
    alpha: f32,
    beta: f32,

    min_integrator: f32,
    max_integrator: f32,

    // Optional extra pole (Type II Order 3): 1st-order IIR on error input.
    // Provides additional high-frequency noise rejection and loop delay
    // compensation. Disabled when pole_a1 = 0, pole_b0 = 1 (passthrough).
    pole_b0: f32,
    pole_a1: f32,
    pole_state: f32,
}

pub fn calculate_gains(loop_bw: f32, damping: f32, k: f32) -> (f32, f32) {
    // Compute in f64 to avoid catastrophic cancellation when ζ×ωn_norm is small.
    // The subtraction 1 - k1×(sinh+cos) loses all significance in f32 for
    // ωn_norm < ~0.001 (typical for timing loops with Bn→ωn conversion).
    let loop_bw = loop_bw as f64;
    let damping = damping as f64;
    let k = k as f64;

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

    (alpha as f32, beta as f32)
}

/// Compute IIR coefficients for an extra filtering pole.
///
/// Given ωp (rad/s) and the loop update rate (Hz), returns (b0, a1) for:
///     filtered[n] = b0 × input[n] + a1 × filtered[n-1]
///
/// DC gain = 1 (steady-state error passes through unchanged).
/// Set omega_p = 0 to disable (returns passthrough: b0=1, a1=0).
pub fn calculate_extra_pole(omega_p: f32, update_rate: f32) -> (f32, f32) {
    if omega_p <= 0.0 {
        return (1.0, 0.0); // passthrough — no filtering
    }
    let a1 = (-omega_p / update_rate).exp();
    let b0 = 1.0 - a1;
    (b0, a1)
}

impl PiLoopFilter {
    /// Create a PI loop filter (Type II, Order 2 — no extra pole).
    pub fn new(loop_bw: f32, damping: f32, k: f32, min_integrator: f32, max_integrator: f32) -> Self {
        Self::new_with_pole(loop_bw, damping, k, min_integrator, max_integrator, 1.0, 0.0)
    }

    /// Create a PI loop filter with an optional extra filtering pole
    /// (Type II, Order 3 when pole is active).
    ///
    /// * `pole_b0`, `pole_a1` — from `calculate_extra_pole()`. Use (1.0, 0.0) to disable.
    pub fn new_with_pole(loop_bw: f32, damping: f32, k: f32,
                         min_integrator: f32, max_integrator: f32,
                         pole_b0: f32, pole_a1: f32) -> Self {
        let (alpha, beta) = calculate_gains(loop_bw, damping, k);
        PiLoopFilter {
            integrator: (min_integrator + max_integrator) / 2.0,
            alpha,
            beta,
            min_integrator,
            max_integrator,
            pole_b0,
            pole_a1,
            pole_state: 0.0,
        }
    }

    /// Create a loop filter with raw coefficients (no Bn/ζ/K derivation).
    /// Use for matching external designs (e.g., liquid-dsp symsync).
    pub fn new_raw(alpha: f32, beta: f32,
                   min_integrator: f32, max_integrator: f32,
                   pole_b0: f32, pole_a1: f32) -> Self {
        PiLoopFilter {
            integrator: (min_integrator + max_integrator) / 2.0,
            alpha,
            beta,
            min_integrator,
            max_integrator,
            pole_b0,
            pole_a1,
            pole_state: 0.0,
        }
    }

    pub fn advance(&mut self, error: f32) -> f32 {
        // Extra pole: 1st-order IIR lowpass on the error signal
        let filtered = self.pole_b0 * error + self.pole_a1 * self.pole_state;
        self.pole_state = filtered;

        // Standard PI on the filtered error
        self.integrator += self.beta * filtered;
        self.integrator = self.integrator.clamp(self.min_integrator, self.max_integrator);

        self.integrator + self.alpha * filtered
    }
}