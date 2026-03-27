use std::cell::RefCell;

use num_complex::Complex32;

use crate::{pi_loop::PiLoopFilter, rds_taps, resample::RationalResampler, symbolsync::SymbolSync};

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

pub struct RdsDemod {
    nco: NcoState,
    symbol_sync: SymbolSync,
    costas_loop_filter: PiLoopFilter,
    downsampler: RationalResampler<Complex32>,

    pub debug: RdsDemodDebug,
}

// Shared constants (must match Python bpsk_pfb_receiver.py)
const R_CHIP: f32 = 2375.0;
const TOTAL_DECIMATE: f32 = 72.0; // PRE_DECIMATE(24) × SPS(3)

// Costas carrier recovery (post-MF, chip rate)
const COSTAS_BN_HZ: f32 = 5.0;
const COSTAS_DAMPING: f32 = 0.707;
const COSTAS_K_DET: f32 = 0.761594; // tanh(1): soft-decision tanh(I)·Q detector at unit power

// Symbol timing recovery (PFB + ML TED)
const TIMING_BN_HZ: f32 = 0.5;
const TIMING_DAMPING: f32 = 1.0;    // critically damped
const TIMING_K_TED: f32 = 0.160671; // ML TED S-curve slope per acc unit (post-AGC, from filter)

const COSTAS_MAX_FREQ_HZ: f32 = 10.0;

/// Bn → ωn → ωn_norm: noise bandwidth to normalized natural frequency
fn bn_to_omega_n_norm(bn_hz: f32, zeta: f32, update_rate: f32) -> f32 {
    let omega_n = 2.0 * bn_hz / (zeta + 1.0 / (4.0 * zeta));
    omega_n / update_rate
}

impl RdsDemod {

    pub fn new() -> Self {
        let mf = rds_taps::generate_rrc_taps(2375.0*3.0*12.0, 2375.0, 0.8, 8);

        let costas_omega_n_norm = bn_to_omega_n_norm(COSTAS_BN_HZ, COSTAS_DAMPING, R_CHIP);
        let costas_max = 2.0 * std::f32::consts::PI * COSTAS_MAX_FREQ_HZ / R_CHIP;

        let symbol_omega_n_norm = bn_to_omega_n_norm(TIMING_BN_HZ, TIMING_DAMPING, R_CHIP);

        let samples_per_symbol = 3;
        let symbol_nfilt = 12;
        let symbol_max_period_deviation = 1.0;

        let downsample_filter = rds_taps::generate_lowpass_taps(171e3, 2500.0, 1001, &rds_taps::WindowType::Blackman);

        RdsDemod {
            nco: NcoState::new(57e3, 171e3),
            symbol_sync: SymbolSync::new(mf, samples_per_symbol, symbol_max_period_deviation, symbol_nfilt, symbol_omega_n_norm, TIMING_DAMPING, TIMING_K_TED),
            costas_loop_filter: PiLoopFilter::new(costas_omega_n_norm, COSTAS_DAMPING, COSTAS_K_DET, -costas_max, costas_max),
            downsampler: RationalResampler::new(downsample_filter, 1, 24),
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

        // Pipeline: NCO mix (171kHz) → downsample 24× (7125Hz) → SymbolSync
        // SymbolSync includes internal post-MF AGC (applied to both MF and dMF
        // before TED, matching Python placement). Output is unit-power.
        let mf_out = {
            let pipeline = iter.map(|sample| nco.borrow_mut().mix(sample));
            let mut symsync_input = downsampler.iter(pipeline);
            symbol_sync.next(&mut symsync_input)?
        }; // iterator dropped

        debug.filtered_samples_hist.push(mf_out);
        debug.agc_hist.push(symbol_sync.agc.gain);

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

        debug.phase_adj_hist.push(pi_output);
        debug.freq_avg_hist.push(costas.integrator);

        Some(mf_out)
    }
}
