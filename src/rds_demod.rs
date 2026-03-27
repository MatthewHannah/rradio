use std::cell::RefCell;

use num_complex::Complex32;

use crate::{fir::Fir, pi_loop::PiLoopFilter, rds_taps, resample::{Downsampleable, RationalResampler}, symbolsync::SymbolSync};

struct NcoAdjustment {
    freq_adj: f32,
    phase_adj: f32,
}

struct AgcState {
    rate: f32,
    reference: f32,
    gain: f32,
    max_gain: f32,
    initial_gain_detected: bool,
    initial_gain_num_accum: usize,
    initial_gain_accum: f32,
}

impl AgcState {
    fn new(rate: f32, reference: f32, max_gain: f32) -> Self {
        AgcState {
            rate,
            reference,
            gain: 1.0,
            max_gain,
            initial_gain_detected: false,
            initial_gain_num_accum: 0,
            initial_gain_accum: 0.0,
        }
    }

    fn process(&mut self, input: Complex32) -> Complex32 {
        if !self.initial_gain_detected {
            self.initial_gain_accum += input.norm();
            self.initial_gain_num_accum += 1;
            if self.initial_gain_num_accum >= 100 {
                let initial_gain = (self.reference) / (self.initial_gain_accum / self.initial_gain_num_accum as f32);
                self.gain = initial_gain.min(self.max_gain);
                self.initial_gain_detected = true;
                println!("AGC initial gain detected: {}, Accum: {}", self.gain, self.initial_gain_accum);
            }
            input
        } else {
            let out = input * self.gain;
            self.gain += self.rate * (self.reference - out.norm());
            self.gain = self.gain.clamp(1e-6, self.max_gain); // prevent negative gain
            out
        }
    }
}

struct NcoState {
    phase: f32,
    base_freq_incr: f32,
    adj: f32,
}

impl NcoState {
    fn new(freq: f32, sample_rate: f32) -> Self {
        NcoState {
            phase: 0.0,
            base_freq_incr: freq / sample_rate * 2.0 * std::f32::consts::PI,
            adj: 0.0,
        }
    }

    fn mix(&mut self, sample: f32) -> Complex32 {
        let step = self.base_freq_incr + self.adj;
        self.phase = (self.phase + step) % (2.0 * std::f32::consts::PI);
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
    agc: AgcState,
    symbol_sync: SymbolSync,
    costas_loop_filter: PiLoopFilter,
    downsampler: RationalResampler<Complex32>,

    pub debug: RdsDemodDebug,
}

impl RdsDemod {

    pub fn new() -> Self {
        let mf = rds_taps::generate_rrc_taps(2375.0*3.0*12.0, 2375.0, 0.8, 8);

        let costas_operating_freq = 7125.0;
        let costas_bw = 100.0; // Hz
        let costas_bw = 2.0 * std::f32::consts::PI * costas_bw / costas_operating_freq; // convert to normalized rad/sample at chip rate
        let costas_min_integrator = -100.0; // Hz
        let costas_min_integrator = 2.0 * std::f32::consts::PI * costas_min_integrator / costas_operating_freq; // convert to normalized rad/sample at full
        let costas_max_integrator = 100.0; // Hz
        let costas_max_integrator = 2.0 * std::f32::consts::PI * costas_max_integrator / costas_operating_freq; // convert to normalized rad/sample at full
        let costas_damping = 0.707;
        let costas_gain = 0.5;

        let symbol_bw = 10.0; // Hz
        let symbol_bw = 2.0 * std::f32::consts::PI * symbol_bw / 2375.0; // convert to normalized rad/sample at chip rate

        let samples_per_symbol = 3;
        let symbol_nfilt = 12;
        let symbol_damping = 0.707;
        let symbol_max_period_deviation = 1.0; // units of expected symbol period
        let k_ted = 1.48; // empirically derived TED gain per chip for BPSK with this matched filter

        let downsample_filter = rds_taps::generate_lowpass_taps(171e3, 2500.0, 1001, &rds_taps::WindowType::Blackman);

        RdsDemod {
            nco: NcoState::new(57e3, 171e3),
            agc: AgcState::new(1e-2, 1.0, 1e5),
            symbol_sync: SymbolSync::new(mf, samples_per_symbol, symbol_max_period_deviation, symbol_nfilt, symbol_bw, symbol_damping, k_ted),
            costas_loop_filter: PiLoopFilter::new(costas_bw, costas_damping, costas_gain, costas_min_integrator, costas_max_integrator),
            downsampler: RationalResampler::new(downsample_filter, 1, 24),
            debug: RdsDemodDebug::default(),
        }
    }

    // expects to receive 171 kHz sample rate
    pub fn next(&mut self, iter: &mut impl Iterator<Item = f32>) -> Option<Complex32> {
        // Destructure self to allow split borrows across the closure and symbol_sync
        let nco = &mut self.nco;
        let nco = RefCell::new(nco);
        let downsampler = &mut self.downsampler;
        let agc = &mut self.agc;
        let debug = &mut self.debug;
        let symbol_sync = &mut self.symbol_sync;

        let pipeline_into_downsampler = iter.map(|sample| nco.borrow_mut().mix(sample));

        let mut pipeline_into_symsync = downsampler
            .iter(pipeline_into_downsampler)
            .map(|sample| {
                debug.agc_hist.push(agc.gain);
                let processed = agc.process(sample);
                debug.filtered_samples_hist.push(processed);

                let phase_error = processed.re * processed.im;
                let _out = self.costas_loop_filter.advance(phase_error) / 24.0;
                let freq_avg = self.costas_loop_filter.integrator;

                nco.borrow_mut().adj = self.costas_loop_filter.integrator / 24.0;
                debug.phase_adj_hist.push(_out);
                debug.freq_avg_hist.push(freq_avg);

                processed
            })
            ;

        // estimate of symbol from timing loop
        //let sample = symbol_sync.next(&mut pipeline_into_symsync)?;
        let sample = pipeline_into_symsync.next()?;



        Some(sample)
    }
}
