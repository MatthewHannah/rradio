use std::cell::RefCell;

use num_complex::Complex32;

use crate::{pi_loop::{self, PiLoopFilter}, rds_taps, resample::RationalResampler, symbolsync::SymbolSync};

struct NcoState {
    phase: f32,
    base_freq_incr: f32,
    freq_adj: f32,     // ongoing frequency correction (rad/sample at F_S)
    phase_adj: f32,    // one-shot phase correction (radians, applied then cleared)
}

impl NcoState {
    fn new(freq: f32, sample_rate: f32) -> Self {
        NcoState {
            phase: 0.0,
            base_freq_incr: freq / sample_rate * 2.0 * std::f32::consts::PI,
            freq_adj: 0.0,
            phase_adj: 0.0,
        }
    }

    fn mix(&mut self, sample: f32) -> Complex32 {
        let step = self.base_freq_incr + self.freq_adj;
        self.phase += step + self.phase_adj;
        self.phase_adj = 0.0; // one-shot: applied on first sample, then zero
        // Wrap to [0, 2π) — Rust's % can produce negative remainders
        let two_pi = 2.0 * std::f32::consts::PI;
        while self.phase >= two_pi { self.phase -= two_pi; }
        while self.phase < 0.0 { self.phase += two_pi; }
        let osc_sample = Complex32::new(self.phase.cos(), -self.phase.sin());
        sample * osc_sample
    }
}

#[derive(Default)]
pub struct RdsDemodDebug {
    pub filtered_samples_hist: Vec<Complex32>,
    pub agc_hist: Vec<f32>,
    pub phase_adj_hist: Vec<f32>,
    pub freq_avg_hist: Vec<f32>,
}

/// Per-chip diagnostic record for dropout analysis
#[derive(Clone, Copy)]
pub struct ChipDiag {
    pub mf_re: f32,
    pub mf_im: f32,
    pub costas_phase_err: f32,
    pub costas_freq: f32,
    pub timing_period: f32,
    pub timing_avg_period: f32,
    pub agc_gain: f32,
    pub input_power: f32,   // 171 kHz input power (before NCO)
    pub ds_power: f32,      // 7125 Hz power (after NCO + downsample, before MF)
}

/// Per-chip biphase diagnostic (captured in pipeline after demod)
pub struct BiphaseDiag {
    pub biphase_energy: f32,
    pub even_sum: f32,
    pub odd_sum: f32,
    pub polarity: u8,
}

pub struct RdsDemod {
    nco: NcoState,
    symbol_sync: SymbolSync,
    costas_loop_filter: PiLoopFilter,
    downsampler: RationalResampler<Complex32>,

    // Power tracking (exponential averages, updated at sample rate)
    pub input_power_avg: f32,   // at 171 kHz, before NCO
    pub ds_power_avg: f32,      // at 7125 Hz, after NCO+downsample, before MF

    pub debug: RdsDemodDebug,
}

// Shared constants (must match Python bpsk_pfb_receiver.py)
const R_CHIP: f32 = 2375.0;
const TOTAL_DECIMATE: f32 = 72.0; // PRE_DECIMATE(24) × SPS(3)

// Costas carrier recovery (post-MF, chip rate)
const COSTAS_BN_HZ: f32 = 30.0;
const COSTAS_DAMPING: f32 = 0.707;
const COSTAS_K_DET: f32 = 0.761594; // tanh(1): soft-decision tanh(I)·Q detector at unit power

// Symbol timing recovery (PFB + ML TED)
const TIMING_BN_HZ: f32 = 25.0;
const TIMING_DAMPING: f32 = 1.0;    // critically damped
const TIMING_K_TED: f32 = 0.040125; // ML TED S-curve slope per acc unit (post-AGC, from filter)
const TIMING_NFILTERS: usize = 32;

const COSTAS_MAX_FREQ_HZ: f32 = 10.0;

// Extra filtering pole multiplier (× ωn). Set to 0.0 to disable.
const COSTAS_POLE_MULT: f32 = 0.0;
const TIMING_POLE_MULT: f32 = 0.0;

// Timing loop mode selection
// 0 = analytical PI (clock_tracking_loop gains, optional extra pole)
// 1 = liquid-dsp empirical coefficients
// 2 = pole-placed IIR+PI (principled, accounts for K_TED)
const TIMING_MODE: u8 = 0;

// liquid-dsp symsync coefficients (from bt = 2200/171000):
//   IIR prefilter: b0 = 0.005589, a1 = 0.964850
//   PI: alpha = 1.0 (full filtered error), beta = 0.006433 (rate_adj)
const LIQUID_BT: f32 = 2200.0 / 171000.0;

// Pole-placed IIR+PI coefficients (from liquid_loop_derivation.py):
// Design at K=1 (clipped TED regime, governs acquisition)
// Specs: T_acq=100 chips (42ms), T_freq=1500 chips (632ms)
// Poles at K=1: z_fast=0.970446, z_slow=0.999334
// IIR pole: z_fast + 0.003 = 0.973446
// Noise BW: 0.65 Hz (K=1), ~6 Hz (K=5.14 locked)
const PP_ALPHA: f32 = 1.0;
const PP_BETA: f32 = 0.005401;
const PP_POLE_B0: f32 = 0.003647;
const PP_POLE_A1: f32 = 0.973446;

/// Bn → ωn: noise bandwidth to natural frequency (rad/s)
fn bn_to_omega_n(bn_hz: f32, zeta: f32) -> f32 {
    2.0 * bn_hz / (zeta + 1.0 / (4.0 * zeta))
}

/// Compute liquid-dsp symsync loop coefficients from normalized bandwidth bt.
fn liquid_symsync_coeffs(bt: f32) -> (f32, f32, f32, f32) {
    let alpha_iir = 1.0 - bt;
    let beta_c = 0.22 * bt;
    let a = 0.5_f32;
    let b = 0.495_f32;
    let a0 = 1.0 - a * alpha_iir;

    let pole_b0 = beta_c / a0;
    let pole_a1 = (b * alpha_iir) / a0;  // positive feedback coeff
    let pi_alpha = 1.0;                   // full filtered error → output
    let pi_beta = 0.5 * bt;              // rate_adjustment

    (pi_alpha, pi_beta, pole_b0, pole_a1)
}

impl RdsDemod {

    pub fn new() -> Self {
        let mf = rds_taps::generate_rrc_taps(2375.0*3.0*32.0, 2375.0, 0.8, 3);

        let costas_omega_n = bn_to_omega_n(COSTAS_BN_HZ, COSTAS_DAMPING);
        let costas_omega_n_norm = costas_omega_n / R_CHIP;
        let costas_max = 2.0 * std::f32::consts::PI * COSTAS_MAX_FREQ_HZ / R_CHIP;
        let costas_omega_p = COSTAS_POLE_MULT * costas_omega_n;
        let (costas_pole_b0, costas_pole_a1) = pi_loop::calculate_extra_pole(costas_omega_p, R_CHIP);

        let samples_per_symbol = 3;
        let symbol_max_period_deviation = 1.0;

        let downsample_filter = rds_taps::generate_lowpass_taps(171e3, 2500.0, 1001, &rds_taps::WindowType::Blackman);

        // Pre-NCO bandpass filter: reject audio/stereo/pilot before 57 kHz mixing.
        // Design: lowpass at 3 kHz, modulated to 57 kHz center.
        // 501 taps at 171 kHz gives ~60 dB rejection of audio transients.

        let symbol_sync = match TIMING_MODE {
            1 => {
                // Liquid-dsp empirical coefficients
                let (alpha, beta, pole_b0, pole_a1) = liquid_symsync_coeffs(LIQUID_BT);
                SymbolSync::new_raw(mf, samples_per_symbol, symbol_max_period_deviation,
                                    TIMING_NFILTERS, alpha, beta, pole_b0, pole_a1)
            }
            2 => {
                // Pole-placed IIR+PI (principled, K_TED-aware)
                SymbolSync::new_raw(mf, samples_per_symbol, symbol_max_period_deviation,
                                    TIMING_NFILTERS, PP_ALPHA, PP_BETA, PP_POLE_B0, PP_POLE_A1)
            }
            _ => {
                // Standard PI from clock_tracking_loop
                let timing_omega_n = bn_to_omega_n(TIMING_BN_HZ, TIMING_DAMPING);
                let timing_omega_n_norm = timing_omega_n / R_CHIP;
                let timing_omega_p = TIMING_POLE_MULT * timing_omega_n;
                let (timing_pole_b0, timing_pole_a1) = pi_loop::calculate_extra_pole(timing_omega_p, R_CHIP);
                let k_ted_per_sample = TIMING_K_TED * TIMING_NFILTERS as f32;
                SymbolSync::new(mf, samples_per_symbol, symbol_max_period_deviation,
                                TIMING_NFILTERS, timing_omega_n_norm, TIMING_DAMPING,
                                k_ted_per_sample, timing_pole_b0, timing_pole_a1)
            }
        };

        RdsDemod {
            nco: NcoState::new(57e3, 171e3),
            symbol_sync,
            costas_loop_filter: PiLoopFilter::new_with_pole(costas_omega_n_norm, COSTAS_DAMPING, COSTAS_K_DET, -costas_max, costas_max, costas_pole_b0, costas_pole_a1),
            downsampler: RationalResampler::new(downsample_filter, 1, 24),
            input_power_avg: 0.0,
            ds_power_avg: 0.0,
            debug: RdsDemodDebug::default(),
        }
    }

    // expects to receive 171 kHz sample rate
    pub fn next(&mut self, iter: &mut impl Iterator<Item = f32>) -> Option<Complex32> {
        // Destructure self for split borrows (Rust 2021 closure capture rules)
        let nco = &mut self.nco;
        let nco = RefCell::new(nco);
        let downsampler = &mut self.downsampler;
        let debug = &mut self.debug;
        let symbol_sync = &mut self.symbol_sync;
        let costas = &mut self.costas_loop_filter;

        // Pipeline: BPF 57kHz → NCO mix (171kHz) → downsample 24× (7125Hz) → SymbolSync
        let mf_out = {
            let pipeline = iter.map(|sample| {
                
                
                
                nco.borrow_mut().mix(sample)
            });
            let mut symsync_input = downsampler.iter(pipeline);
            symbol_sync.next(&mut symsync_input)?
        }; // iterator dropped

        // debug.filtered_samples_hist.push(mf_out);
        // debug.agc_hist.push(symbol_sync.agc.gain);

        // Costas carrier recovery (once per chip, post-MF)
        // Soft-decision detector: tanh(I)·Q, K_det = tanh(1) ≈ 0.762
        let phase_error = mf_out.re.tanh() * mf_out.im;
        let pi_output = costas.advance(phase_error);

        // Multi-rate NCO update:
        //   freq (ongoing): integrator in rad/chip → ÷72 → rad/F_S_sample
        //   phase (one-shot): proportional term (α×error) in radians
        {
            let mut nco_ref = nco.borrow_mut();
            nco_ref.freq_adj = costas.integrator / TOTAL_DECIMATE;
            nco_ref.phase_adj = pi_output - costas.integrator; // = α × error
        }

        // debug.phase_adj_hist.push(pi_output);
        // debug.freq_avg_hist.push(costas.integrator);

        Some(mf_out)
    }

    /// Same as next() but also returns per-chip diagnostic data.
    pub fn next_diag(&mut self, iter: &mut impl Iterator<Item = f32>) -> Option<(Complex32, ChipDiag)> {
        let nco = &mut self.nco;
        let nco = RefCell::new(nco);
        let downsampler = &mut self.downsampler;
        let symbol_sync = &mut self.symbol_sync;
        let costas = &mut self.costas_loop_filter;
        let input_pwr = RefCell::new(&mut self.input_power_avg);
        let ds_pwr = RefCell::new(&mut self.ds_power_avg);

        let mf_out = {
            let pipeline = iter.map(|sample| {
                let p = sample * sample;
                let avg = **input_pwr.borrow();
                **input_pwr.borrow_mut() = 0.9999 * avg + 0.0001 * p;
                
                
                
                nco.borrow_mut().mix(sample)
            });
            let mut symsync_input = downsampler.iter(pipeline).map(|s| {
                let p = s.re * s.re + s.im * s.im;
                let avg = **ds_pwr.borrow();
                **ds_pwr.borrow_mut() = 0.999 * avg + 0.001 * p;
                s
            });
            symbol_sync.next(&mut symsync_input)?
        };

        let phase_error = mf_out.re.tanh() * mf_out.im;
        let pi_output = costas.advance(phase_error);

        let diag = ChipDiag {
            mf_re: mf_out.re,
            mf_im: mf_out.im,
            costas_phase_err: phase_error,
            costas_freq: costas.integrator,
            timing_period: symbol_sync.inst_period,
            timing_avg_period: symbol_sync.avg_period,
            agc_gain: symbol_sync.agc.gain,
            input_power: **input_pwr.borrow(),
            ds_power: **ds_pwr.borrow(),
        };

        {
            let mut nco_ref = nco.borrow_mut();
            nco_ref.freq_adj = costas.integrator / TOTAL_DECIMATE;
            nco_ref.phase_adj = pi_output - costas.integrator;
        }

        Some((mf_out, diag))
    }
}
