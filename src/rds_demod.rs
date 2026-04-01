use std::cell::RefCell;
use num_complex::Complex32;

use crate::{fir, resample::RationalResampler};

// ── Constants ──
const R_CHIP: f32 = 2375.0;
const F_BASE: f32 = 14250.0;         // Decimated rate (6 SPS)
const SPS: usize = 6;                // Samples per chip at F_BASE
const PRE_DECIMATE: usize = 12;      // 171000 / 14250 = 12
const RRC_ALPHA: f64 = 0.8;
const RRC_SPAN: usize = 3;

// Costas loop
const COSTAS_BN_HZ: f32 = 30.0;
const COSTAS_DAMPING: f32 = 0.707;
const COSTAS_K_DET: f32 = 0.761594;  // tanh(1) for tanh(I)×Q detector
const COSTAS_MAX_FREQ_HZ: f32 = 100.0;

// Timing loop
const TIMING_BN_HZ: f32 = 25.0;
const TIMING_DAMPING: f32 = 1.0;
const TIMING_K_TED: f32 = 0.6152;     // Gardner S-curve slope: 3.69 / SPS
const TIMING_NFILTERS: usize = 16;
const TIMING_TAPS_PER_ARM: usize = 8;

// ── Gain computation ──
fn compute_pi_gains(loop_bw: f32, damping: f32, k_det: f32, update_rate: f32) -> (f32, f32) {
    let loop_bw = loop_bw as f64;
    let damping = damping as f64;
    let k = k_det as f64;
    let update_rate = update_rate as f64;

    let omega_n = 2.0 * loop_bw / (damping + 1.0 / (4.0 * damping));
    let omega_n_norm = omega_n / update_rate;
    let zeta_wn = omega_n_norm * damping;

    let k0 = 2.0 / k;
    let k1 = (-zeta_wn).exp();
    let sinh_val = zeta_wn.sinh();

    let cos_val = if (damping - 1.0).abs() < 1e-6 {
        1.0
    } else if damping < 1.0 {
        (omega_n_norm * (1.0 - damping * damping).sqrt()).cos()
    } else {
        (omega_n_norm * (damping * damping - 1.0).sqrt()).cosh()
    };

    let alpha = k0 * k1 * sinh_val;
    let beta = k0 * (1.0 - k1 * (sinh_val + cos_val));
    (alpha as f32, beta as f32)
}

// ── Coarse NCO ──
struct CoarseNco {
    phase: f32,
    freq_incr: f32,
}

impl CoarseNco {
    fn new(freq: f32, sample_rate: f32) -> Self {
        CoarseNco {
            phase: 0.0,
            freq_incr: 2.0 * std::f32::consts::PI * freq / sample_rate,
        }
    }

    fn mix(&mut self, sample: f32) -> Complex32 {
        self.phase += self.freq_incr;
        let two_pi = 2.0 * std::f32::consts::PI;
        while self.phase >= two_pi { self.phase -= two_pi; }
        while self.phase < 0.0 { self.phase += two_pi; }
        sample * Complex32::new(self.phase.cos(), -self.phase.sin())
    }
}

// ── Fine Costas Loop (runs at 14250 Hz) ──
struct FineCostas {
    // RRC matched filter
    h_rrc: Vec<f32>,
    rrc_buf: Vec<Complex32>,
    rrc_idx: usize,

    // Fine NCO
    nco_phase: f32,
    nco_freq: f32,

    // PI loop
    alpha: f32,
    beta: f32,
    freq_integrator: f32,
    max_freq: f32,
}

impl FineCostas {
    fn new() -> Self {
        let h_rrc = fir::generate_rrc_taps(
            F_BASE as f64 * 1.0, // fs = F_BASE (SPS is already set by the rate)
            R_CHIP as f64,
            RRC_ALPHA,
            RRC_SPAN,
        );
        // Actually: generate_rrc_taps(fs, symbol_rate, beta, num_symbols)
        // fs = 14250, symbol_rate = 2375, beta = 0.8, num_symbols = 3
        // This gives taps at SPS = 14250/2375 = 6, span = ±3 chips
        // len = 2*3*6+1 = 37 taps
        let h_rrc = fir::generate_rrc_taps(
            F_BASE as f64, R_CHIP as f64, RRC_ALPHA, RRC_SPAN,
        );
        let rrc_len = h_rrc.len();

        let (alpha, beta) = compute_pi_gains(
            COSTAS_BN_HZ, COSTAS_DAMPING, COSTAS_K_DET, F_BASE,
        );
        let max_freq = 2.0 * std::f32::consts::PI * COSTAS_MAX_FREQ_HZ / F_BASE;

        FineCostas {
            rrc_buf: vec![Complex32::new(0.0, 0.0); rrc_len],
            rrc_idx: 0,
            h_rrc,
            nco_phase: 0.0,
            nco_freq: 0.0,
            alpha,
            beta,
            freq_integrator: 0.0,
            max_freq,
        }
    }

    /// Process one complex sample at 14250 Hz.
    /// Returns the carrier-corrected, MF-shaped sample.
    fn process(&mut self, sample: Complex32) -> Complex32 {
        // Fine NCO correction
        let corrected = sample * Complex32::new(self.nco_phase.cos(), -self.nco_phase.sin());
        self.nco_phase += self.nco_freq;
        let two_pi = 2.0 * std::f32::consts::PI;
        while self.nco_phase > std::f32::consts::PI { self.nco_phase -= two_pi; }
        while self.nco_phase < -std::f32::consts::PI { self.nco_phase += two_pi; }

        // Push into RRC filter
        self.rrc_buf[self.rrc_idx] = corrected;
        self.rrc_idx = (self.rrc_idx + 1) % self.rrc_buf.len();

        // Evaluate RRC (circular convolution)
        let rrc_len = self.h_rrc.len();
        let mut mf_out = Complex32::new(0.0, 0.0);
        for j in 0..rrc_len {
            let buf_idx = (self.rrc_idx + j) % rrc_len;
            mf_out += self.rrc_buf[buf_idx] * self.h_rrc[j];
        }

        // Phase error: tanh(I) × Q
        let phase_error = mf_out.re.tanh() * mf_out.im;

        // PI loop
        self.freq_integrator += self.beta * phase_error;
        self.freq_integrator = self.freq_integrator.clamp(-self.max_freq, self.max_freq);
        self.nco_freq = self.freq_integrator;
        self.nco_phase += self.alpha * phase_error;

        mf_out
    }
}

// ── AGC ──
struct Agc {
    rate: f32,
    target_power: f32,
    avg_power: f32,
    pub gain: f32,
    max_gain: f32,
}

impl Agc {
    fn new(rate: f32, target_rms: f32, max_gain: f32) -> Self {
        Agc {
            rate,
            target_power: target_rms * target_rms,
            avg_power: 1.0,
            gain: 1.0,
            max_gain,
        }
    }

    fn process(&mut self, sample: Complex32) -> Complex32 {
        let power = sample.re * sample.re + sample.im * sample.im;
        self.avg_power = (1.0 - self.rate) * self.avg_power + self.rate * power;
        if self.avg_power > 1e-20 {
            self.gain = (self.target_power / self.avg_power).sqrt().min(self.max_gain);
        }
        sample * self.gain
    }
}

// ── Polyphase Gardner Timing Recovery ──
struct PolyphaseGardner {
    // Interpolating filter (lowpass, not RRC)
    interp_arms: Vec<Vec<f32>>,  // nfilters × taps_per_arm
    nfilters: usize,
    taps_per_arm: usize,

    // Circular delay line
    buf: Vec<Complex32>,
    head: usize,
    buf_len: usize,

    // Phase accumulator
    phase: f32,
    fire_threshold: f32,   // nfilters × SPS

    // Period tracking
    avg_period: f32,
    inst_period: f32,
    inst_rate: f32,
    min_period: f32,
    max_period: f32,

    // PI loop
    alpha: f32,
    beta: f32,

    // AGC (post-MF, pre-TED)
    agc: Agc,
}

impl PolyphaseGardner {
    fn new() -> Self {
        let nfilters = TIMING_NFILTERS;
        let taps_per_arm = TIMING_TAPS_PER_ARM;
        let sps = SPS;

        // Design lowpass interpolating filter (sinc + Blackman)
        let proto_len = nfilters * taps_per_arm;
        let center = (proto_len as f32 - 1.0) / 2.0;
        let mut proto = vec![0.0_f32; proto_len];
        for i in 0..proto_len {
            let n = i as f32 - center;
            // Sinc
            let sinc_val = if n.abs() < 1e-6 { 1.0 } else {
                let x = std::f32::consts::PI * n / nfilters as f32;
                x.sin() / x
            };
            // Blackman window
            let w = 0.42 - 0.5 * (2.0 * std::f32::consts::PI * i as f32 / (proto_len - 1) as f32).cos()
                + 0.08 * (4.0 * std::f32::consts::PI * i as f32 / (proto_len - 1) as f32).cos();
            proto[i] = sinc_val * w;
        }

        // Decompose into polyphase arms
        let mut arms = vec![vec![0.0_f32; taps_per_arm]; nfilters];
        for k in 0..nfilters {
            for j in 0..taps_per_arm {
                let idx = k + j * nfilters;
                if idx < proto_len {
                    arms[k][j] = proto[idx];
                }
            }
            // Normalize arm to unity DC gain
            let sum: f32 = arms[k].iter().sum();
            if sum.abs() > 1e-10 {
                for c in arms[k].iter_mut() { *c /= sum; }
            }
        }

        let fire_threshold = (nfilters * sps) as f32;
        let avg_period = sps as f32;
        let inst_rate = fire_threshold / avg_period;

        let (alpha, beta) = compute_pi_gains(
            TIMING_BN_HZ, TIMING_DAMPING, TIMING_K_TED, R_CHIP,
        );

        let buf_len = taps_per_arm + sps + 2;

        PolyphaseGardner {
            interp_arms: arms,
            nfilters,
            taps_per_arm,
            buf: vec![Complex32::new(0.0, 0.0); buf_len * 2],
            head: 0,
            buf_len,
            phase: 0.0,
            fire_threshold,
            avg_period,
            inst_period: avg_period,
            inst_rate,
            min_period: sps as f32 - 1.0,
            max_period: sps as f32 + 1.0,
            alpha,
            beta,
            agc: Agc::new(0.01, 1.0, 1e5),
        }
    }

    /// Interpolate at a fractional delay from the newest sample.
    fn interp(&self, frac_delay: f32) -> Complex32 {
        let int_delay = frac_delay.floor() as usize;
        let frac = frac_delay - int_delay as f32;
        let arm = ((frac * self.nfilters as f32).round() as usize) % self.nfilters;

        let start = (self.head as isize - 1 - int_delay as isize
                     - self.taps_per_arm as isize + 1)
            .rem_euclid(self.buf_len as isize) as usize;

        let mut result = Complex32::new(0.0, 0.0);
        for j in 0..self.taps_per_arm {
            let idx = (start + j) % self.buf_len;
            result += self.buf[idx] * self.interp_arms[arm][j];
        }
        result
    }

    /// Push one complex sample at 14250 Hz. Returns a chip when fire occurs.
    fn push_sample(&mut self, sample: Complex32) -> Option<Complex32> {
        // AGC on input
        let sample = self.agc.process(sample);

        // Push into delay line
        self.buf[self.head] = sample;
        if self.head + self.buf_len < self.buf.len() {
            self.buf[self.head + self.buf_len] = sample;
        }
        self.head = (self.head + 1) % self.buf_len;

        // Advance phase
        self.phase += self.inst_rate;

        if self.phase >= self.fire_threshold {
            self.phase -= self.fire_threshold;

            // Compute fractional delay for the three Gardner samples
            let overshoot = self.phase / self.nfilters as f32;

            let y_now = self.interp(overshoot);
            let y_mid = self.interp(overshoot + self.inst_period / 2.0);
            let y_prev = self.interp(overshoot + self.inst_period);

            // Gardner error (negated: positive → decrease period)
            let error = -(y_mid.conj() * (y_now - y_prev)).re;
            let error = error.clamp(-4.0, 4.0);

            // PI loop
            self.avg_period += self.beta * error;
            self.avg_period = self.avg_period.clamp(self.min_period, self.max_period);
            self.inst_period = self.avg_period + self.alpha * error;
            if self.inst_period <= 0.0 {
                self.inst_period = self.avg_period;
            }
            self.inst_rate = self.fire_threshold / self.inst_period;

            return Some(y_now);
        }

        None
    }
}

// ── Top-level demodulator ──
pub struct RdsDemod {
    coarse_nco: CoarseNco,
    downsampler: RationalResampler<Complex32>,
    costas: FineCostas,
    agc_pre: Agc,    // pre-MF AGC for Costas
    gardner: PolyphaseGardner,
}

impl RdsDemod {
    pub fn new() -> Self {
        let downsample_filter = fir::generate_lowpass_taps(
            171e3, 2500.0, 1001, &fir::WindowType::Blackman,
        );

        RdsDemod {
            coarse_nco: CoarseNco::new(57e3, 171e3),
            downsampler: RationalResampler::new(downsample_filter, 1, PRE_DECIMATE),
            costas: FineCostas::new(),
            agc_pre: Agc::new(0.001, 1.0, 1e5),
            gardner: PolyphaseGardner::new(),
        }
    }

    /// Process input samples at 171 kHz. Returns one chip when available.
    pub fn next(&mut self, iter: &mut impl Iterator<Item = f32>) -> Option<Complex32> {
        let coarse_nco = &mut self.coarse_nco;
        let coarse_nco = RefCell::new(coarse_nco);
        let downsampler = &mut self.downsampler;
        let costas = &mut self.costas;
        let agc_pre = &mut self.agc_pre;
        let gardner = &mut self.gardner;

        loop {
            // Pull one decimated sample: NCO mix → LPF → decimate 12×
            let dec_sample = {
                let pipeline = iter.map(|sample| coarse_nco.borrow_mut().mix(sample));
                match downsampler.iter(pipeline).next() {
                    Some(s) => s,
                    None => return None,
                }
            };

            // AGC → Costas (carrier recovery + MF) at 14250 Hz
            let agc_sample = agc_pre.process(dec_sample);
            let mf_sample = costas.process(agc_sample);

            // Gardner timing recovery → fires at chip rate
            if let Some(chip) = gardner.push_sample(mf_sample) {
                return Some(chip);
            }
        }
    }
}
