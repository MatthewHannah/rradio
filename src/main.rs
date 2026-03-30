mod sigmf;
mod osc;
mod biquad;
mod resample;
mod fm_demod;
mod playback;
mod deemphasis;
mod filterable;
mod spy;
mod pll;
mod interleaver;
mod wideband_fm_audio;
mod pluto;
mod buffer;
mod soapy;
mod fir;
mod gardner_clock_recovery;
mod manchester;
mod rds_block_sync;
mod rds_config;
mod rds_taps;
mod agc;
mod costas;
mod nco;
mod biphase;
mod psk_modem;
mod symsync;
mod symbolsync;
mod rds_demod;
mod pi_loop;

use std::sync::atomic;
use std::sync::Arc;

use num_complex::Complex;
use num_integer::Integer;
use plotly::Histogram;
use plotly::ImageFormat;
use plotly::{HeatMap, Plot, Scatter};
use num_complex::Complex32;
use rustfft::{num_traits::Zero, FftPlanner};

use crate::costas::CostasDemodulable;
use crate::filterable::Filter;
use crate::filterable::FilterableIter;
use crate::fm_demod::FmDemodulatable;
use crate::gardner_clock_recovery::GardnerClockRecoverable;
use crate::interleaver::InterleaveableIter;
use crate::manchester::ManchesterDecodable;
use crate::rds_block_sync::RdsBlockSyncable;
use crate::rds_config::CostasConfig;
use crate::resample::{Downsampleable, RationalResampleable};
use crate::spy::SpyableIter;
use crate::wideband_fm_audio::WidebandFmAudioIterable;
use crate::osc::Mixable;

#[allow(dead_code)]
fn spectrogram<T>(window_size: usize, overlap: usize, fs: f32, samples: &[T]) where Complex32: From<T>, T: Copy {
    let mut ffts: Vec<Vec<f32>> = vec![];
    let mut start: usize = 0;
    let mut stop: usize = start + window_size;
    let incr = window_size - overlap;

    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(window_size);
    let mut holder = vec![Complex32::zero(); window_size];

    while stop < samples.len() {
        let mut working: Vec<Complex32> = samples[start..stop].iter().cloned().map(|x| Complex32::from(x)).collect();
        fft.process(&mut working);

        holder[window_size/2..].copy_from_slice(&working[..window_size/2]);
        holder[..window_size/2].copy_from_slice(&working[window_size/2..]);
        ffts.push(holder.iter().map(|s| s.norm().log10() * 20.0).collect());

        start = start + incr;
        stop = start + window_size;
    }

    if stop != samples.len() {
        let mut remainder: Vec<Complex32> = samples[start..].iter().cloned().map(|x| Complex32::from(x)).collect();
        remainder.resize(window_size, Complex32::zero());
        fft.process(&mut remainder);

        holder[window_size/2..].copy_from_slice(&remainder[..window_size/2]);
        holder[..window_size/2].copy_from_slice(&remainder[window_size/2..]);
        ffts.push(holder.iter().map(|s| s.norm().log10() * 20.0).collect());
    }
    println!("FFT done!");

    let x: Vec<f32> = (0..window_size).map(|i| i as f32 * fs / window_size as f32 - fs / 2.0).collect();
    let y: Vec<f32> = (0..ffts.len()).map(|i| i as f32).collect();

    let mut plot = Plot::new();
    let trace = HeatMap::new(x.clone(), y.clone(), ffts.clone());
    plot.add_trace(trace);
    plot.show();
}

/// Save a PSD (power spectral density) plot to a PNG file.
#[allow(dead_code)]
fn save_psd_png(samples: &[Complex32], fs: f32, title: &str, path: &str) {
    let mut planner = FftPlanner::<f32>::new();
    let n = samples.len().min(8192);
    let fft = planner.plan_fft_forward(n);

    let mut working: Vec<Complex32> = samples[..n].to_vec();
    fft.process(&mut working);

    let mut holder = vec![Complex32::zero(); n];
    holder[n/2..].copy_from_slice(&working[..n/2]);
    holder[..n/2].copy_from_slice(&working[n/2..]);
    let magnitudes: Vec<f32> = holder.iter().map(|s| s.norm().log10() * 20.0).collect();
    let freqs: Vec<f32> = (0..n).map(|i| i as f32 * fs / n as f32 - fs / 2.0).collect();

    let mut plot = Plot::new();
    let trace = Scatter::new(freqs, magnitudes);
    plot.add_trace(trace);
    plot.set_layout(plotly::Layout::new().title(title));
    plot.write_image(path, plotly::ImageFormat::PNG, 1200, 600, 1.0);
    eprintln!("Saved PSD: {}", path);
}

/// Save a PSD of real-valued samples to a PNG file.
#[allow(dead_code)]
fn save_real_psd_png(samples: &[f32], fs: f32, title: &str, path: &str) {
    let complex: Vec<Complex32> = samples.iter().map(|&s| Complex32::new(s, 0.0)).collect();
    save_psd_png(&complex, fs, title, path);
}

#[allow(dead_code)]
fn plot_freq_resp(samples: &[Complex32], fs: f32) {
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(samples.len());
    let mut holder = vec![Complex32::zero(); samples.len()];

    let mut working: Vec<Complex32> = samples.iter().cloned().collect();
    fft.process(&mut working);
    let holder_len = holder.len();
    holder[holder_len/2..].copy_from_slice(&working[..holder_len/2]);
    holder[..holder_len/2].copy_from_slice(&working[holder_len/2..]);
    let magnitudes = holder.iter().map(|s| s.norm().log10() * 20.0).collect::<Vec<_>>();

    let freqs: Vec<f32> = (0..samples.len()).map(|i| i as f32 * fs / samples.len() as f32 - fs / 2.0).collect();

    let mut plot = Plot::new();
    let trace = Scatter::new(freqs, magnitudes);
    plot.add_trace(trace);
    plot.show();
}

#[allow(dead_code)]
fn plot_re_im(samples: &[Complex32]) {
    let sample_count: Vec<usize> = (0..samples.len()).collect();
    let r: Vec<f32> = samples.iter().map(|s| s.re).collect();
    let i: Vec<f32> = samples.iter().map(|s| s.im).collect();
    let mut plot = Plot::new();
    let rtrace: Box<Scatter<usize, f32>> = Scatter::new(sample_count.clone(), r).into();
    let itrace: Box<Scatter<usize, f32>> = Scatter::new(sample_count.clone(), i).into();
    plot.add_trace(rtrace);
    plot.add_trace(itrace);
    plot.show();
}

#[allow(dead_code)]
fn plot_re(samples: &[f32], title: &str) {
    let sample_count: Vec<usize> = (0..samples.len()).collect();
    let mut plot = Plot::new();
    let rtrace: Box<Scatter<usize, f32>> = Scatter::new(sample_count.clone(), samples.to_vec()).into();
    plot.add_trace(rtrace);
    plot.set_layout(plotly::Layout::new().title(title));
    plot.show();
}

fn buffer_pluto(done: &atomic::AtomicBool, config: &pluto::SdrConfig, mut out: buffer::SendBuf<Vec<Complex32>>) {
    const MAX_RETRIES: u32 = 3;

    while !done.load(atomic::Ordering::SeqCst) {
        // (Re)connect and configure SDR
        let mut sdr = match pluto::PlutoSdr::connect(config) {
            Ok(sdr) => {
                eprintln!("SDR connected");
                sdr
            }
            Err(e) => {
                eprintln!("SDR connect failed: {:?}, retrying in 1s...", e);
                std::thread::sleep(std::time::Duration::from_secs(1));
                continue;
            }
        };

        let mut streamer = match sdr.start_iq() {
            Ok(s) => s,
            Err(e) => {
                eprintln!("SDR start_iq failed: {:?}, retrying in 1s...", e);
                std::thread::sleep(std::time::Duration::from_secs(1));
                continue;
            }
        };

        // Stream until error or shutdown
        let mut consecutive_errors: u32 = 0;
        while !done.load(atomic::Ordering::SeqCst) {
            if let Some(mut buf) = out.get() {
                match streamer.collect_iq(&mut buf) {
                    Ok(()) => {
                        consecutive_errors = 0;
                        out.commit(buf);
                    }
                    Err(e) => {
                        consecutive_errors += 1;
                        if consecutive_errors >= MAX_RETRIES {
                            eprintln!("SDR error: {:?} ({} consecutive failures), reconnecting...", e, consecutive_errors);
                            break;
                        }
                        eprintln!("SDR error: {:?} (attempt {}/{})", e, consecutive_errors, MAX_RETRIES);
                        std::thread::sleep(std::time::Duration::from_millis(100));
                    }
                }
            } else {
                // SendBuf returned None — receiver is gone
                return;
            }
        }

        // Tear down before reconnecting
        sdr.stop_iq(streamer);

        if !done.load(atomic::Ordering::SeqCst) {
            eprintln!("SDR connection lost, reconnecting in 1s...");
            std::thread::sleep(std::time::Duration::from_secs(1));
        }
    }
}

fn buffer_soapy(done: &atomic::AtomicBool, config: &soapy::SoapyConfig, mut out: buffer::SendBuf<Vec<Complex32>>) {
    const MAX_RETRIES: u32 = 3;

    while !done.load(atomic::Ordering::SeqCst) {
        let sdr = match soapy::SoapySdr::connect(config) {
            Ok(sdr) => {
                eprintln!("SoapySDR connected");
                sdr
            }
            Err(e) => {
                eprintln!("SoapySDR connect failed: {}, retrying in 1s...", e);
                std::thread::sleep(std::time::Duration::from_secs(1));
                continue;
            }
        };

        let mut streamer = match sdr.start_iq() {
            Ok(s) => s,
            Err(e) => {
                eprintln!("SoapySDR start_iq failed: {}, retrying in 1s...", e);
                std::thread::sleep(std::time::Duration::from_secs(1));
                continue;
            }
        };

        let mut consecutive_errors: u32 = 0;
        while !done.load(atomic::Ordering::SeqCst) {
            if let Some(mut buf) = out.get() {
                match streamer.collect_iq(&mut buf) {
                    Ok(()) => {
                        consecutive_errors = 0;
                        out.commit(buf);
                    }
                    Err(e) => {
                        consecutive_errors += 1;
                        if consecutive_errors >= MAX_RETRIES {
                            eprintln!("SoapySDR error: {} ({} consecutive failures), reconnecting...", e, consecutive_errors);
                            break;
                        }
                        eprintln!("SoapySDR error: {} (attempt {}/{})", e, consecutive_errors, MAX_RETRIES);
                        std::thread::sleep(std::time::Duration::from_millis(100));
                    }
                }
            } else {
                return;
            }
        }

        if let Err(e) = streamer.deactivate() {
            eprintln!("SoapySDR deactivate error: {}", e);
        }

        if !done.load(atomic::Ordering::SeqCst) {
            eprintln!("SoapySDR connection lost, reconnecting in 1s...");
            std::thread::sleep(std::time::Duration::from_secs(1));
        }
    }
}

fn buffer_sigmf(done: &atomic::AtomicBool, mut streamer: sigmf::SigmfStreamer, mut out: buffer::SendBuf<Vec<Complex32>>, tune_offset: f32, fs: f32) {
    const CHUNK: usize = 4*1024*1024;
    while !done.load(atomic::Ordering::SeqCst) {
        if let Some(mut buf) = out.get() {
            buf.clear();
            buf.extend((&mut streamer).take(CHUNK).mix(tune_offset, fs));
            if buf.is_empty() {
                break;
            }
            out.commit(buf);
        }
    }
}

struct AudioPipelineObservationSettings {
    spy_iq: bool,
    spy_demoded: bool,
    spy_audio: bool,
}

struct SignalPipelineSettings {
    iq_downsample: usize,
    fm_demod_downsample: usize,
}

fn signal_pipeline(
    done: &atomic::AtomicBool,
    fs: f32,
    inbuf: buffer::RecvBuf<Vec<Complex32>>,
    mut audio_out: buffer::SendBuf<Vec<(f32, f32)>>,
    mut rds_out: buffer::SendBuf<Vec<f32>>,
    settings: SignalPipelineSettings,
    obs_settings: AudioPipelineObservationSettings,
    mpx_path: Option<String>,
) {
    let samples = buffer::RecvBufIter::new(inbuf);

    let fs_spy = fs;
    let samples = samples.maybe_spy(6000000, move |iq_samples| {
        spectrogram(8192, 512, fs_spy, &iq_samples);
    }, obs_settings.spy_iq);

    // IQ filtering + downsample
    let iq_filter_stages = 2;
    let iq_filters: Vec<biquad::Biquad<Complex32>> = (0..iq_filter_stages)
        .map(|_| biquad::Biquad::lowpass(fs, 300000.0, 0.707))
        .collect();
    let mut filtered: Box<dyn Iterator<Item = Complex32>> = Box::new(samples);
    for filt in iq_filters {
        filtered = Box::new(filtered.dsp_filter(filt));
    }
    let fs: f32 = fs / (settings.iq_downsample as f32);
    let resampled = filtered.downsample(settings.iq_downsample);

    let fs_spy = fs;
    let resampled = resampled.maybe_spy(6000000, move |audio_samples| {
        spectrogram(8192, 512, fs_spy, &audio_samples);
    }, obs_settings.spy_iq);

    // FM demodulation + downsample
    let fm_filt = biquad::Biquad::lowpass(fs, 80000.0, 0.707);
    let fm_loop_filt = biquad::Biquad::lowpass(fs, 80000.0, 0.707);
    let demoded = resampled.dsp_filter(fm_filt).fm_demodulate(fm_loop_filt).downsample(settings.fm_demod_downsample);
    let _fs = fs / (settings.fm_demod_downsample as f32);

    let fs_spy = _fs;
    let demoded = demoded.maybe_spy(48000, move |audio_samples| {
        spectrogram(1024, 512, fs_spy, &audio_samples);
    }, obs_settings.spy_demoded);

    // Optional MPX output: write FM-demodulated baseband to WAV
    let mut mpx_writer = mpx_path.map(|path| {
        let spec = hound::WavSpec {
            channels: 1,
            sample_rate: _fs as u32,
            bits_per_sample: 16,
            sample_format: hound::SampleFormat::Int,
        };
        hound::WavWriter::create(&path, spec)
            .unwrap_or_else(|e| panic!("Failed to create MPX file {}: {}", path, e))
    });

    // Tap demoded signal for MPX writing, then feed to wfm_audio
    let demoded = demoded.inspect(move |&s| {
        if let Some(ref mut writer) = mpx_writer {
            // Scale to i16 range — FM demod output is roughly ±1
            let sample = (s * 16000.0).clamp(-32767.0, 32767.0) as i16;
            let _ = writer.write_sample(sample);
        }
    });

    // Wideband FM audio (stereo + RDS extraction) — tee to two consumers
    let mut wfm = demoded.wfm_audio(_fs);

    let mut audio_batch: Option<buffer::BufToken<Vec<(f32, f32)>>> = None;
    let mut rds_batch: Option<buffer::BufToken<Vec<f32>>> = None;

    while !done.load(atomic::Ordering::SeqCst) {
        // Ensure we have output buffers
        if audio_batch.is_none() {
            audio_batch = match audio_out.get() {
                Some(mut tok) => { tok.clear(); Some(tok) }
                None => break,
            };
        }
        if rds_batch.is_none() {
            rds_batch = rds_out.get().map(|mut tok| { tok.clear(); tok });
        }

        let sample = match wfm.next() {
            Some(s) => s,
            None => break,
        };

        // Tee: audio gets (left, right), RDS gets raw MPX
        if let Some(ref mut buf) = audio_batch {
            buf.push((sample.left, sample.right));
            if buf.len() >= 4096 {
                audio_out.commit(audio_batch.take().unwrap());
            }
        }
        if let Some(ref mut buf) = rds_batch {
            buf.push(sample.mpx);
            if buf.len() >= 4096 {
                rds_out.commit(rds_batch.take().unwrap());
            }
        }
    }

    // Flush remaining
    if let Some(buf) = audio_batch {
        if !buf.is_empty() { audio_out.commit(buf); }
    }
    if let Some(buf) = rds_batch {
        if !buf.is_empty() { rds_out.commit(buf); }
    }
}

fn compute_pipeline_settings(fs: f32) -> SignalPipelineSettings {
    // Total downsample from IQ to wfm_audio stage: must reach 240 kHz
    // wfm_audio runs at fs / iq_downsample / fm_demod_downsample
    // Audio consumer does ÷5 → 48 kHz, so wfm rate must be 240 kHz
    let total_downsample = fs / 48000.0;
    let total_downsample_int = total_downsample.round() as usize;
    if (total_downsample - total_downsample_int as f32).abs() > 0.5 {
        panic!("Sample rate {} does not divide evenly to 48 kHz (ratio {})", fs, total_downsample);
    }
    // Total = iq_downsample × fm_demod_downsample × audio_downsample(5)
    // We need total / 5 for the signal pipeline portion (iq × fm_demod)
    if total_downsample_int % 5 != 0 {
        panic!("Sample rate {} requires total downsample of {} which is not divisible by 5", fs, total_downsample_int);
    }
    let signal_downsample = total_downsample_int / 5; // iq × fm_demod
    if signal_downsample % 5 != 0 {
        panic!("Signal downsample {} is not divisible by 5 for fm_demod stage", signal_downsample);
    }
    let iq_downsample = signal_downsample / 5;
    if iq_downsample == 0 {
        panic!("Sample rate {} is too low for this pipeline (need at least 1.2 MHz)", fs);
    }
    SignalPipelineSettings {
        iq_downsample,
        fm_demod_downsample: 5,
    }
}

enum IqSource {
    Pluto { config: pluto::SdrConfig },
    Soapy { config: soapy::SoapyConfig },
    Sigmf { streamer: sigmf::SigmfStreamer, tune_offset: f32 },
}

enum AudioOutput {
    Playback,
    Wav(String),
}

fn write_wav<I>(samples: I, path: &str, fs: u32) where I: Iterator<Item = f32> {
    let spec = hound::WavSpec {
        channels: 2,
        sample_rate: fs,
        bits_per_sample: 32,
        sample_format: hound::SampleFormat::Float,
    };
    let mut writer = hound::WavWriter::create(path, spec)
        .unwrap_or_else(|e| panic!("Failed to create WAV file {}: {}", path, e));

    let mut count: u64 = 0;
    for sample in samples {
        writer.write_sample(sample).expect("Failed to write WAV sample");
        count += 1;
    }
    writer.finalize().expect("Failed to finalize WAV file");
    eprintln!("Wrote {} samples ({:.1}s) to {}", count, count as f64 / (fs as f64 * 2.0), path);
}

fn rds_pipeline_v4(done: &atomic::AtomicBool, rds_rx: buffer::RecvBuf<Vec<f32>>, wfm_fs: f32, config: &rds_config::RdsConfig, debug: bool, metrics: bool, diag_path: Option<&str>) {
    let iterable = buffer::RecvBufIter::new(rds_rx);

    let stage1_target_fs: f32 = 57e3 * 3.0;
    let wfm_fs_u64 = wfm_fs as u64;
    let lcm = (stage1_target_fs as u64).lcm(&wfm_fs_u64) as f64;
    let stage1_up = (lcm / (wfm_fs as f64)) as usize;
    let stage1_down = (lcm / (stage1_target_fs as f64)) as usize;
    let up_fs = wfm_fs * (stage1_up as f32);

    println!("Stage 1 resample: {} → {} (up {} down {})", wfm_fs, stage1_target_fs, stage1_up, stage1_down);

    let mut base = iterable
        .resample(rds_taps::generate_lowpass_taps(up_fs as f64, 80e3, 255, &rds_taps::WindowType::Blackman), stage1_up, stage1_down);

    let mut rds = rds_demod::RdsDemod::new();

    let mut biphase_decoder = biphase::BiphaseDecoder::new();
    let mut delta_decoder = biphase::DeltaDecoder::new();

    let mut block_sync = rds_block_sync::RdsBlockSync::new(
        std::iter::empty::<u8>(), &config.sync, debug);
    let mut decoder = rds_block_sync::RdsDecoder::new();
    let mut display = rds_block_sync::RdsDisplay::new();
    let start_time = std::time::Instant::now();

    let mut constellation: Vec<Complex32> = vec![];
    let mut constellation_num = 0;

    // Diagnostic output: per-chip CSV of loop internals
    let mut diag_writer = diag_path.map(|path| {
        let mut w = std::io::BufWriter::new(std::fs::File::create(path).expect("Failed to create diag file"));
        use std::io::Write;
        writeln!(w, "chip,mf_re,mf_im,costas_phase_err,costas_freq,timing_period,timing_avg_period,agc_gain,input_power,ds_power,bit,data_bit").unwrap();
        w
    });
    let mut chip_idx: u64 = 0;

    loop {
        if done.load(atomic::Ordering::SeqCst) {
            break;
        }

        let (sym, diag) = if diag_writer.is_some() {
            match rds.next_diag(&mut base) {
                Some((s, d)) => (s, Some(d)),
                None => break,
            }
        } else {
            match rds.next(&mut base) {
                Some(s) => (s, None),
                None => break,
            }
        };

        // Feedforward timing slip info to biphase decoder
        biphase_decoder.notify_timing_adjust(
            rds.symbol_sync.last_samples_consumed,
            rds.symbol_sync.nominal_sps,
        );

        // Biphase + delta decode
        let result = biphase_decoder.push(sym);
        let data_bit = if result.has_value {
            Some(delta_decoder.decode(result.bit))
        } else {
            None
        };

        if let (Some(w), Some(d)) = (&mut diag_writer, diag) {
            use std::io::Write;
            let biphase_bit: i8 = if result.has_value { if result.bit { 1 } else { 0 } } else { -1 };
            let data_bit_val: i8 = data_bit.map(|b| b as i8).unwrap_or(-1);
            writeln!(w, "{},{:.6},{:.6},{:.6},{:.8},{:.6},{:.6},{:.4},{:.8},{:.8},{},{}",
                chip_idx, d.mf_re, d.mf_im, d.costas_phase_err, d.costas_freq,
                d.timing_period, d.timing_avg_period, d.agc_gain, d.input_power, d.ds_power,
                biphase_bit, data_bit_val).unwrap();
        }
        chip_idx += 1;

        // Block sync
        if let Some(data_bit) = data_bit {
            if let Some(event) = block_sync.push_bit(data_bit as u8) {
                match event {
                    rds_block_sync::SyncEvent::Group(group) => {
                        let state = decoder.process(&group);
                        let elapsed = start_time.elapsed().as_secs_f64();
                        if metrics {
                            let pi_str = if group.pi_code != 0 { format!("{:04X}", group.pi_code) } else { "0000".to_string() };
                            eprintln!("RDSMETRIC {{\"t\":{:.3},\"bler\":{:.4},\"groups\":{},\"pi\":\"{}\"}}", elapsed, group.rolling_bler, state.groups_decoded, pi_str);
                        }
                        if debug {
                            let v = if group.version { "B" } else { "A" };
                            eprintln!("RDS Group {}{}: PI=0x{:04X}  BLER={:.1}%", group.group_type, v, group.pi_code, group.rolling_bler * 100.0);
                        } else if !metrics {
                            display.set_synced(true);
                            display.render(&state);
                        }
                    }
                    rds_block_sync::SyncEvent::Locked => {
                        if !metrics && !debug { display.set_synced(true); display.render(&decoder.display_state()); }
                    }
                    rds_block_sync::SyncEvent::LostSync | rds_block_sync::SyncEvent::Searching => {
                        if !metrics && !debug { display.set_synced(false); display.render(&decoder.display_state()); }
                    }
                }
            }
        }

    }
    // plot_re(&rds.debug.agc_hist[101..], "AGC history");
    // plot_re(&rds.debug.freq_avg_hist, "Costas frequency adjustment history");
    // plot_re(&rds.debug.phase_adj_hist, "Costas phase adjustment history");
    // save_psd_png(&rds.debug.filtered_samples_hist, 171e3/24.0, "Filtered Samples", "filt.png");

    if metrics {
        let elapsed = start_time.elapsed().as_secs_f64();
        let state = decoder.display_state();
        eprintln!("RDSSUMMARY {{\"groups\":{},\"duration\":{:.1},\"final_bler\":{:.4}}}",
            state.groups_decoded, elapsed, 0.0);
    }

}

fn rds_pipeline_v3(done: &atomic::AtomicBool, rds_rx: buffer::RecvBuf<Vec<f32>>, wfm_fs: f32, config: &rds_config::RdsConfig, debug: bool, metrics: bool) {
    let iterable = buffer::RecvBufIter::new(rds_rx);

    use plotly::{Plot, Scatter, common::Mode};

    // need to mix 57kHz to baseband
    let stage1_target_fs: f32 = 57e3 * 2.0;
    let wfm_fs_u64 = wfm_fs as u64;
    let lcm = (stage1_target_fs as u64).lcm(&wfm_fs_u64) as f64;
    let stage1_up = (lcm / (wfm_fs as f64)) as usize;
    let stage1_down = (lcm / (stage1_target_fs as f64)) as usize;
    let proto_cutoff = 100e3;
    let stage1_taps = rds_taps::generate_lowpass_taps(lcm, proto_cutoff, 570, &rds_taps::WindowType::Blackman);
    println!("Stage 1 resample: {} → {} (up {} down {})", wfm_fs, stage1_target_fs, stage1_up, stage1_down);

    // then get down to ~3 samples per chip
    let stage2_target_fs = 2375.0 * 3.0;
    let stage2_up = 1usize;
    let stage2_down = (stage1_target_fs / (stage2_target_fs)) as usize;
    let stage2_taps = rds_taps::generate_lowpass_taps(stage1_target_fs as f64, 2400.0, 255, &rds_taps::WindowType::Blackman);
    let stage2_target_fs = stage2_target_fs as f32;
    println!("Stage 2 resample: {} → {} (up {} down {})", stage1_target_fs, stage2_target_fs, stage2_up, stage2_down);

    let agc_bw = 500.0 / 171e3;
    let agc_initial_gain = 0.08;

    let symsync_bw: f32 = 0.013;
    const SYMSYNC_DELAY: usize = 8;
    const SYMSYNC_BETA: f32 = 0.8;
    const SYMSYNC_NPFB: usize = 32;

    let start_time = std::time::Instant::now();

    let mut symsync = symsync::SymSync::new(3, SYMSYNC_DELAY, SYMSYNC_BETA, SYMSYNC_NPFB);
    symsync.set_bandwidth(symsync_bw);
    symsync.set_output_rate(1);

    let _stage1_fs = stage1_target_fs as f32;
    let _stage2_fs = stage2_target_fs as f32;

    let mut base = iterable
        .resample(stage1_taps, stage1_up, stage1_down)
        .mix(-57e3, stage1_target_fs)
        .resample(stage2_taps, stage2_up, stage2_down)
        .dsp_filter(agc::Agc::new_liquid(agc_bw, agc_initial_gain));

    let costas_config = rds_config::CostasConfig { loop_bw: 0.08 };
    let mut costas = costas::CostasLoop::new(&costas_config);

    let mut biphase_decoder = biphase::BiphaseDecoder::new();
    let mut delta_decoder = biphase::DeltaDecoder::new();

    let mut block_sync = rds_block_sync::RdsBlockSync::new(
        std::iter::empty::<u8>(), &config.sync, debug);
    let mut decoder = rds_block_sync::RdsDecoder::new();
    let mut display = rds_block_sync::RdsDisplay::new();
    let start_time = std::time::Instant::now();

    let mut total_bits: usize = 0;

    // let mut constellation_pre_pre_vec = vec![];
    // let mut constellation_pre_vec = vec![];
    // let mut constellation_post_vec = vec![];

    loop {
        if done.load(atomic::Ordering::Relaxed) { break; }

        let sample = match base.next() {
            Some(s) => s,
            None => break,
        };

        // if constellation_pre_pre_vec.len() < 5000 {
        //     constellation_pre_pre_vec.push(sample);
        //     if constellation_pre_pre_vec.len() == 5000 {
        //         let mut plot = Plot::new();
        //         let trace = Scatter::new(constellation_pre_pre_vec.iter().map(|s| s.re).collect(), constellation_pre_pre_vec.iter().map(|s| s.im).collect()).mode(Mode::Markers);
        //         plot.add_trace(trace);
        //         plot.set_layout(plotly::Layout::new().title("Constellation Pre Sym Sync"));
        //         plot.write_image("constellation_pre_sync.png", plotly::ImageFormat::PNG, 600, 600, 1.0);
        //     }
        // }


        let syms = symsync.execute(sample);
        if syms.is_empty() { continue; }

        // if constellation_pre_vec.len() < 5000 {
        //     constellation_pre_vec.push(sample);
        //     if constellation_pre_vec.len() == 5000 {
        //         let mut plot = Plot::new();
        //         let trace = Scatter::new(constellation_pre_vec.iter().map(|s| s.re).collect(), constellation_pre_vec.iter().map(|s| s.im).collect()).mode(Mode::Markers);
        //         plot.add_trace(trace);
        //         plot.set_layout(plotly::Layout::new().title("Constellation Pre Costas"));
        //         plot.write_image("constellation_pre.png", plotly::ImageFormat::PNG, 600, 600, 1.0);
        //     }
        // }

        for &sym in syms {
            let sym = costas.process(sym);

            // if constellation_post_vec.len() < 5000 {
            //     constellation_post_vec.push(sym);
            //     if constellation_post_vec.len() == 5000 {
            //         let mut plot = Plot::new();
            //         let trace = Scatter::new(constellation_post_vec.iter().map(|s| s.re).collect(), constellation_post_vec.iter().map(|s| s.im).collect()).mode(Mode::Markers);
            //         plot.add_trace(trace);
            //         plot.set_layout(plotly::Layout::new().title("Constellation Post Costas"));
            //         plot.write_image("constellation_post.png", plotly::ImageFormat::PNG, 600, 600, 1.0);
            //     }
            // }


            let result = biphase_decoder.push(sym);
            if result.has_value {
                let data_bit = delta_decoder.decode(result.bit);

                total_bits += 1;

                if let Some(event) = block_sync.push_bit(data_bit as u8) {
                    match event {
                        rds_block_sync::SyncEvent::Group(group) => {
                            let state = decoder.process(&group);
                            let elapsed = start_time.elapsed().as_secs_f64();
                            if metrics {
                                let pi_str = if group.pi_code != 0 { format!("{:04X}", group.pi_code) } else { "0000".to_string() };
                                eprintln!("RDSMETRIC {{\"t\":{:.3},\"bler\":{:.4},\"groups\":{},\"pi\":\"{}\"}}", elapsed, group.rolling_bler, state.groups_decoded, pi_str);
                            }
                            if debug {
                                let v = if group.version { "B" } else { "A" };
                                eprintln!("RDS Group {}{}: PI=0x{:04X}  BLER={:.1}%", group.group_type, v, group.pi_code, group.rolling_bler * 100.0);
                            } else if !metrics {
                                display.set_synced(true);
                                display.render(&state);
                            }
                        }
                        rds_block_sync::SyncEvent::Locked => {
                            if !metrics && !debug { display.set_synced(true); display.render(&decoder.display_state()); }
                        }
                        rds_block_sync::SyncEvent::LostSync | rds_block_sync::SyncEvent::Searching => {
                            if !metrics && !debug { display.set_synced(false); display.render(&decoder.display_state()); }
                        }
                    }
                }

            }
        }
    }

    if metrics {
        let elapsed = start_time.elapsed().as_secs_f64();
        let state = decoder.display_state();
        eprintln!("RDSSUMMARY {{\"groups\":{},\"duration\":{:.1},\"final_bler\":{:.4}}}",
            state.groups_decoded, elapsed, 0.0);
    }

}

const AUDIO_DOWNSAMPLE: usize = 5;
const RDS_FIRST_DECIMATE: usize = 12; // 240 kHz → 20 kHz
const RDS_UPSAMPLE: usize = 19;       // 20 kHz → 380 kHz (virtual)
const RDS_SECOND_DECIMATE: usize = 20; // 380 kHz → 19 kHz
const RDS_FS: f32 = 19000.0;          // final RDS sample rate


/// RDS decoding pipeline.
///
/// The approach here — using a Manchester-shaped RRC matched filter with
/// symbol sync at bit rate (1187.5 Hz) and a Costas loop for BPSK phase
/// recovery — is informed by Bastian Bloessl's gr-rds flowgraph for
/// GNU Radio, which was shown to be the best-performing RDS decoder across
/// 38 FM stations in a comparative analysis by site2241.net (January 2024).
///
/// References:
///   - Bloessl's gr-rds: https://github.com/bastibl/gr-rds
///   - Comparative analysis: https://www.site2241.net/january2024.htm
///   - Andy Walls, "Symbol Clock Recovery", GRCon 2017 (AGC requirements for TEDs)
///   - RDS standard: IEC 62106
/// Redsea-style RDS pipeline.
///
/// Input: Complex32 baseband at wfm_fs (240 kHz) from PLL-driven 57 kHz mixer.
///
/// Architecture (mimicking windytan/redsea):
///   FIR LPF 2400 Hz → ÷ to ~3 samp/chip → AGC → SymSync (complex clock recovery)
///   → PSK2 modem → phase error → NCO PLL feedback → biphase decode → delta decode
///
/// Key differences from our old pipeline:
///   - Tight 2400 Hz LPF (was 4000 Hz)
///   - Chip-rate clock recovery at 2375 Hz (was bit-rate at 1187.5 Hz)
///   - Decision-directed carrier tracking via modem phase error (was Costas loop)
///   - Explicit biphase Manchester decode (was Manchester-RRC matched filter)
///
/// TODO: For true Redsea replication, the pipeline should receive real MPX (f32)
/// and do its own 57 kHz NCO mixing with PLL feedback. Currently we receive
/// post-PLL complex baseband, so the phase error feedback path is not yet closed.
/// Redsea-compatible RDS pipeline.
///
/// Exact replication of windytan/redsea's DSP chain using liquid-dsp algorithms:
///   Resample -> NCO 57kHz (PLL) -> FIR LPF 2400Hz -> div24 -> AGC -> SymSync
///   -> Modem phase error -> PLL feedback -> Biphase -> Delta -> Block sync
fn rds_pipeline(done: &atomic::AtomicBool, rds_rx: buffer::RecvBuf<Vec<f32>>, wfm_fs: f32, config: &rds_config::RdsConfig, debug: bool, metrics: bool, diag_path: Option<&str>) {
    // === REDSEA CONSTANTS ===
    const TARGET_FS: f32 = 171_000.0;
    const BITS_PER_SECOND: f32 = 1187.5;
    const CHIP_RATE: f32 = 2375.0;
    const SAMPLES_PER_SYMBOL: usize = 3;
    const DECIMATE_RATIO: usize = 24; // 171000 / 2375 / 3 = 24
    const AGC_BW: f32 = 500.0 / TARGET_FS;
    const AGC_INITIAL_GAIN: f32 = 0.08;
    const LPF_CUTOFF: f32 = 2400.0 / TARGET_FS;
    const LPF_LEN: usize = 255;
    const SYMSYNC_BW: f32 = 2200.0 / TARGET_FS;
    const SYMSYNC_DELAY: usize = 3;
    const SYMSYNC_BETA: f32 = 0.8;
    const SYMSYNC_NPFB: usize = 32;
    const PLL_BW_HZ: f32 = 0.03;
    const PLL_MULTIPLIER: f32 = 12.0;

    // === RESAMPLE 240 kHz -> 171 kHz (L=57, M=80) ===
    let resamp_l = 57_usize;
    let resamp_m = 80_usize;
    let proto_fs = resamp_l as f64 * wfm_fs as f64;
    let proto_cutoff = TARGET_FS as f64 / 2.0;
    let mut resamp_taps = rds_taps::generate_lowpass_taps(proto_fs, proto_cutoff, 570, &rds_taps::WindowType::Blackman);
    for t in resamp_taps.iter_mut() { *t *= resamp_l as f32; }
    let mut resampler = resample::RationalResampler::<f32>::new(resamp_taps, resamp_l, resamp_m);

    // === NCO at 57 kHz with PLL (liquid-dsp: alpha=bw, beta=sqrt(bw)) ===
    let mut nco = nco::Nco::new(57000.0, TARGET_FS, PLL_BW_HZ, CHIP_RATE);

    // === Polyphase decimator: LPF + ÷24 in one step (L=1, M=24) ===
    let mut lpf_taps = rds_taps::generate_lowpass_taps(TARGET_FS as f64, 2400.0, LPF_LEN, &rds_taps::WindowType::Blackman);
    let fir_scale = 2.0 * LPF_CUTOFF;
    for t in lpf_taps.iter_mut() { *t *= fir_scale; }  // bake scale into taps
    let mut decimator = resample::RationalResampler::<Complex32>::new(lpf_taps, 1, DECIMATE_RATIO);

    // === AGC (liquid-dsp compatible) ===
    let mut agc = agc::Agc::new_liquid(AGC_BW, AGC_INITIAL_GAIN);

    // === SymSync: polyphase filterbank clock recovery ===
    let mut symsync = symsync::SymSync::new(SAMPLES_PER_SYMBOL, SYMSYNC_DELAY, SYMSYNC_BETA, SYMSYNC_NPFB);
    symsync.set_bandwidth(SYMSYNC_BW);
    symsync.set_output_rate(1);

    // === Modem + Biphase + Delta ===
    let modem = psk_modem::BpskModem::new();
    let mut biphase_decoder = biphase::BiphaseDecoder::new();
    let mut delta_decoder = biphase::DeltaDecoder::new();

    // === Block Sync ===
    let mut block_sync = rds_block_sync::RdsBlockSync::new(
        std::iter::empty::<u8>(), &config.sync, debug);
    let mut decoder = rds_block_sync::RdsDecoder::new();
    let mut display = rds_block_sync::RdsDisplay::new();
    let start_time = std::time::Instant::now();

    let mut total_chips: u64 = 0;
    let mut total_bits: u64 = 0;

    // Diagnostic collectors (only when running from recording, not live)
    let collect_diag = !metrics; // Skip diagnostics in metrics mode for cleaner benchmarks
    let mut diag_phase_err: Vec<f32> = Vec::new();
    let mut diag_strobe_re: Vec<f32> = Vec::new();
    let mut diag_strobe_im: Vec<f32> = Vec::new();

    // Chip-level CSV diagnostic output
    let mut chip_writer = diag_path.map(|path| {
        let mut w = std::io::BufWriter::new(std::fs::File::create(path).expect("Failed to create diag file"));
        use std::io::Write;
        writeln!(w, "chip,re,im").unwrap();
        w
    });

    let decimated_fs = TARGET_FS / DECIMATE_RATIO as f32;
    eprintln!("Redsea pipeline: resample {:.0}->{:.0} -> NCO 57kHz -> polyphase LPF+div{} ({}t) -> {:.0} Hz -> AGC -> SymSync(k={},b={},npfb={}) -> BPSK -> PLL(bw={:.2e},x{})",
        wfm_fs, TARGET_FS, DECIMATE_RATIO, LPF_LEN, decimated_fs,
        SAMPLES_PER_SYMBOL, SYMSYNC_BETA, SYMSYNC_NPFB, PLL_BW_HZ, PLL_MULTIPLIER);

    use crate::filterable::Filter;
    let mut rds_iter = buffer::RecvBufIter::new(rds_rx);

    let mut constellation_pre_symsync_vec = vec![];
    let mut constellation_post_symsync_vec = vec![];

    // === MAIN PROCESSING LOOP ===
    // Chain: rds_iter → resample 240k→171k → NCO mix → polyphase LPF+÷24 → 7125 Hz
    loop {
        if done.load(atomic::Ordering::Relaxed) { break; }

        // Polyphase decimator pulls from NCO-mixed stream, which pulls from resampler
        let filtered = match decimator.process(&mut std::iter::from_fn(|| {
            let mpx_sample = resampler.process(&mut rds_iter)?;
            let baseband = nco.mix_down(mpx_sample);
            nco.step();
            Some(baseband)
        })) {
            Some(s) => s,
            None => break,
        };

        // --- Running at 7125 Hz (3 samp/chip) ---

        // 4. AGC
        let agc_out = agc.process(filtered);

        use plotly::{Plot, Scatter, Layout, common::Mode};


        if constellation_pre_symsync_vec.len() < 5000 {
            constellation_pre_symsync_vec.push(agc_out);
            if constellation_pre_symsync_vec.len() == 500 {
                let mut plot = Plot::new();
                let trace = Scatter::new(constellation_pre_symsync_vec.iter().map(|s| s.re).collect(), constellation_pre_symsync_vec.iter().map(|s| s.im).collect()).mode(Mode::Markers);
                plot.add_trace(trace);
                plot.set_layout(plotly::Layout::new().title("Constellation Pre Symsync"));
                plot.write_image("constellation_pre.png", plotly::ImageFormat::PNG, 600, 600, 1.0);
            }
        }

        // 5. SymSync: polyphase filterbank clock recovery
        let symbols = symsync.execute(agc_out);

        for &symbol in symbols.iter() {
            if constellation_post_symsync_vec.len() < 5000 {
                constellation_post_symsync_vec.push(symbol);
                if constellation_post_symsync_vec.len() == 5000 {
                    let mut plot = Plot::new();
                    let trace = Scatter::new(constellation_post_symsync_vec.iter().map(|s| s.re).collect(), constellation_post_symsync_vec.iter().map(|s| s.im).collect()).mode(Mode::Markers);
                    plot.add_trace(trace);
                    plot.set_layout(plotly::Layout::new().title("Constellation Post Symsync"));
                    plot.write_image("constellation_post.png", plotly::ImageFormat::PNG, 600, 600, 1.0);
                }
            }
            // --- Running at ~2375 Hz (chip rate) ---
            total_chips += 1;

            // Dump chip to CSV if diagnostic enabled
            if let Some(ref mut w) = chip_writer {
                use std::io::Write;
                writeln!(w, "{},{:.6},{:.6}", total_chips, symbol.re, symbol.im).unwrap();
            }

            // 6. BPSK modem: phase error (decision ignored, only error used)
            let demod = modem.demodulate(symbol);

            // Collect diagnostics (skip in metrics mode)
            if collect_diag {
                //diag_phase_err.push(demod.phase_error);
                //diag_strobe_re.push(symbol.re);
                //diag_strobe_im.push(symbol.im);
            }

            // 7. PLL feedback: phase error * multiplier -> NCO
            let clamped = demod.phase_error.clamp(-std::f32::consts::PI, std::f32::consts::PI);
            nco.step_pll(clamped * PLL_MULTIPLIER);

            // 8. Biphase (Manchester) decode: 2375 -> 1187.5 Hz
            let biphase_result = biphase_decoder.push(symbol);
            if biphase_result.has_value {
                // 9. Delta (differential) decode
                let data_bit = delta_decoder.decode(biphase_result.bit);
                total_bits += 1;

                // Dump bit to chip CSV (reuse same file)
                if let Some(ref mut w) = chip_writer {
                    use std::io::Write;
                    writeln!(w, "BIT,{},{}", total_bits, data_bit as u8).unwrap();
                }

                // 10. Block sync + group decode
                if let Some(event) = block_sync.push_bit(data_bit as u8) {
                    match event {
                        rds_block_sync::SyncEvent::Group(group) => {
                            let state = decoder.process(&group);
                            let elapsed = start_time.elapsed().as_secs_f64();
                            if metrics {
                                let pi_str = if group.pi_code != 0 { format!("{:04X}", group.pi_code) } else { "0000".to_string() };
                                eprintln!("RDSMETRIC {{\"t\":{:.3},\"bler\":{:.4},\"groups\":{},\"pi\":\"{}\"}}", elapsed, group.rolling_bler, state.groups_decoded, pi_str);
                            }
                            if debug {
                                let v = if group.version { "B" } else { "A" };
                                eprintln!("RDS Group {}{}: PI=0x{:04X}  BLER={:.1}%", group.group_type, v, group.pi_code, group.rolling_bler * 100.0);
                            } else if !metrics {
                                display.set_synced(true);
                                display.render(&state);
                            }
                        }
                        rds_block_sync::SyncEvent::Locked => {
                            if !metrics && !debug { display.set_synced(true); display.render(&decoder.display_state()); }
                        }
                        rds_block_sync::SyncEvent::LostSync | rds_block_sync::SyncEvent::Searching => {
                            if !metrics && !debug { display.set_synced(false); display.render(&decoder.display_state()); }
                        }
                    }
                }
            }
        }
    }

    eprintln!("Pipeline stats: {} chips, {} bits", total_chips, total_bits);

    // === Diagnostic plots ===
    if !diag_phase_err.is_empty() {
        use plotly::{Plot, Scatter, Layout};
        use plotly::layout::Axis;
        let time: Vec<f32> = (0..diag_phase_err.len()).map(|i| i as f32 / CHIP_RATE).collect();

        let mut plot = Plot::new();
        plot.add_trace(Scatter::new(time, diag_phase_err.clone()).name("Phase error (rad)"));
        plot.set_layout(Layout::new()
            .title("07 PLL Phase Error vs Time")
            .x_axis(Axis::new().title("Time (s)"))
            .y_axis(Axis::new().title("rad").range(vec![-3.2, 3.2])));
        plot.write_image("07_pll_phase_error.png", plotly::ImageFormat::PNG, 1400, 400, 1.0);
        eprintln!("Saved: 07_pll_phase_error.png");

        let mut plot = Plot::new();
        use plotly::common::Mode;
        plot.add_trace(Scatter::new(diag_strobe_re.clone(), diag_strobe_im.clone())
            .mode(Mode::Markers).name("Symbol IQ"));
        plot.set_layout(Layout::new().title("09 Constellation (chip strobes)"));
        plot.write_image("09_constellation.png", plotly::ImageFormat::PNG, 600, 600, 1.0);
        eprintln!("Saved: 09_constellation.png");
    }

    if metrics {
        let elapsed = start_time.elapsed().as_secs_f64();
        let state = decoder.display_state();
        eprintln!("RDSSUMMARY {{\"groups\":{},\"duration\":{:.1},\"final_bler\":{:.4}}}",
            state.groups_decoded, elapsed, 0.0);
    }
}

fn run(iq_source: IqSource, audio_output: AudioOutput, done_sig: Arc<atomic::AtomicBool>, obs_settings: AudioPipelineObservationSettings, rds_debug: bool, rds_metrics: bool, record_path: Option<String>, rds_config: rds_config::RdsConfig, mpx_path: Option<String>, diag_path: Option<String>) {
    let fs = match &iq_source {
        IqSource::Pluto { config } => config.fs,
        IqSource::Soapy { config } => config.fs,
        IqSource::Sigmf { streamer, .. } => streamer.sample_rate(),
    };
    let station_freq = match &iq_source {
        IqSource::Pluto { config } => config.station as f64,
        IqSource::Soapy { config } => config.station as f64,
        IqSource::Sigmf { .. } => 0.0,
    };

    let settings = compute_pipeline_settings(fs);
    let wfm_fs = fs / (settings.iq_downsample as f32) / (settings.fm_demod_downsample as f32);
    eprintln!("Pipeline: fs={} Hz, iq_ds={}, fm_ds={}, wfm_fs={} Hz",
        fs, settings.iq_downsample, settings.fm_demod_downsample, wfm_fs);
    eprintln!("  Audio: wfm @ {} Hz → ÷{} → {} Hz stereo",
        wfm_fs, AUDIO_DOWNSAMPLE, wfm_fs / AUDIO_DOWNSAMPLE as f32);
    eprintln!("  RDS:   wfm @ {} Hz → ÷{} → ↑{} ↓{} → {} Hz baseband",
        wfm_fs, RDS_FIRST_DECIMATE, RDS_UPSAMPLE, RDS_SECOND_DECIMATE, RDS_FS);

    // Buffer pairs
    let (iq_tx, iq_rx) = buffer::buf_pair::<Vec<Complex32>>(8);
    let (audio_tx, audio_rx) = buffer::buf_pair::<Vec<(f32, f32)>>(8);
    let (rds_tx, rds_rx) = buffer::buf_pair::<Vec<f32>>(8);

    // Optional IQ recording: splitter tees raw IQ to both signal pipeline and recorder
    let (pipeline_rx, record_thread) = if let Some(ref path) = record_path {
        let (pipeline_tx, pipeline_rx) = buffer::buf_pair::<Vec<Complex32>>(8);
        let (record_tx, record_rx) = buffer::buf_pair::<Vec<Complex32>>(4);

        let hw = match &iq_source {
            IqSource::Pluto { .. } => "PlutoSDR",
            IqSource::Soapy { .. } => "RTL-SDR via SoapySDR",
            IqSource::Sigmf { .. } => "SigMF playback",
        };
        let mut writer = sigmf::SigmfWriter::new(path, fs as f64, station_freq, hw)
            .expect("Failed to create SigMF recording");

        // Splitter thread: reads IQ, writes to pipeline + recorder
        let done_ref = done_sig.clone();
        let mut pipeline_out = pipeline_tx;
        let mut record_out = record_tx;
        let splitter = std::thread::spawn(move || {
            let mut rx = iq_rx;
            while !done_ref.load(atomic::Ordering::SeqCst) {
                let token = match rx.get() {
                    Some(t) => t,
                    None => break,
                };

                // Write to recorder (blocking — never drop)
                if let Some(mut rec_tok) = record_out.get() {
                    rec_tok.clear();
                    rec_tok.extend_from_slice(&token);
                    record_out.commit(rec_tok);
                }

                // Write to pipeline (blocking)
                if let Some(mut pipe_tok) = pipeline_out.get() {
                    pipe_tok.clear();
                    pipe_tok.extend_from_slice(&token);
                    pipeline_out.commit(pipe_tok);
                } else {
                    break;
                }

                rx.release(token);
            }
        });

        // Recorder thread: writes IQ to disk
        let done_ref = done_sig.clone();
        let rec_path = path.clone();
        let recorder = std::thread::spawn(move || {
            let iter = buffer::RecvBufIter::new(record_rx);
            let mut batch = Vec::with_capacity(8192);
            for sample in iter {
                if done_ref.load(atomic::Ordering::SeqCst) { break; }
                batch.push(sample);
                if batch.len() >= 8192 {
                    if let Err(e) = writer.write_samples(&batch) {
                        eprintln!("SigMF write error: {}", e);
                        break;
                    }
                    batch.clear();
                }
            }
            if !batch.is_empty() {
                let _ = writer.write_samples(&batch);
            }
            if let Err(e) = writer.finalize() {
                eprintln!("SigMF finalize error: {}", e);
            }
        });

        // Splitter runs in background, we need to join both later
        // Store splitter handle alongside recorder
        let combined = std::thread::spawn(move || {
            splitter.join().unwrap();
            recorder.join().unwrap();
        });

        (pipeline_rx, Some(combined))
    } else {
        (iq_rx, None)
    };

    // Thread 1: IQ source
    let done_ref = done_sig.clone();
    let iq_thread = std::thread::spawn(move || {
        match iq_source {
            IqSource::Pluto { config } => buffer_pluto(&done_ref, &config, iq_tx),
            IqSource::Soapy { config } => buffer_soapy(&done_ref, &config, iq_tx),
            IqSource::Sigmf { streamer, tune_offset } => buffer_sigmf(&done_ref, streamer, iq_tx, tune_offset, fs),
        }
    });

    // Thread 2: Signal pipeline (FM demod + stereo/RDS extraction → tee)
    let done_ref = done_sig.clone();
    let signal_thread = std::thread::spawn(move || {
        signal_pipeline(&done_ref, fs, pipeline_rx, audio_tx, rds_tx, settings, obs_settings, mpx_path);
    });

    // Thread 3: RDS consumer
    let done_ref = done_sig.clone();
    let rds_config2 = rds_config;
    let rds_thread = std::thread::spawn(move || {
        // Toggle: use rds_pipeline (redsea) or rds_pipeline_v4 (ours)
        let use_redsea = std::env::var("USE_REDSEA").is_ok();
        if use_redsea {
            rds_pipeline(&done_ref, rds_rx, wfm_fs, &rds_config2, rds_debug, rds_metrics, diag_path.as_deref());
        } else {
            rds_pipeline_v4(&done_ref, rds_rx, wfm_fs, &rds_config2, rds_debug, rds_metrics, diag_path.as_deref());
        }
    });

    // Main thread: Audio consumer (downsample + interleave + output)
    let audio_iter = buffer::RecvBufIter::new(audio_rx)
        .downsample(AUDIO_DOWNSAMPLE)
        .interleave();

    match audio_output {
        AudioOutput::Playback => {
            playback::playback_iter(audio_iter, 48000);
        }
        AudioOutput::Wav(path) => {
            write_wav(audio_iter, &path, 48000);
        }
    }
    done_sig.store(true, atomic::Ordering::SeqCst);

    iq_thread.join().unwrap();
    signal_thread.join().unwrap();
    rds_thread.join().unwrap();
    if let Some(rt) = record_thread {
        rt.join().unwrap();
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();

    let done_sig = Arc::new(atomic::AtomicBool::new(false));

    let done_ref = done_sig.clone();
    ctrlc::set_handler(move || {
        done_ref.store(true, atomic::Ordering::SeqCst);
    }).expect("Error setting CTRL-C handler");

    let obs_settings = AudioPipelineObservationSettings {
        spy_iq: false,
        spy_demoded: false,
        spy_audio: false,
    };

    // Extract -- args
    let mut positional = Vec::new();
    let mut wav_path: Option<String> = None;
    let mut rds_debug = false;
    let mut rds_metrics = false;
    let mut rds_config_path: Option<String> = None;
    let mut record_path: Option<String> = None;
    let mut mpx_path: Option<String> = None;
    let mut diag_path: Option<String> = None;
    let mut duration_secs: Option<f64> = None;
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--wav" {
            wav_path = Some(args.get(i + 1).expect("Usage: --wav <path>").clone());
            i += 2;
        } else if args[i] == "--rds-debug" {
            rds_debug = true;
            i += 1;
        } else if args[i] == "--rds-metrics" {
            rds_metrics = true;
            i += 1;
        } else if args[i] == "--rds-config" {
            rds_config_path = Some(args.get(i + 1).expect("Usage: --rds-config <path>").clone());
            i += 2;
        } else if args[i] == "--record" {
            record_path = Some(args.get(i + 1).expect("Usage: --record <path>").clone());
            i += 2;
        } else if args[i] == "--mpx" {
            mpx_path = Some(args.get(i + 1).expect("Usage: --mpx <path.wav>").clone());
            i += 2;
        } else if args[i] == "--diag" {
            diag_path = Some(args.get(i + 1).expect("Usage: --diag <path.csv>").clone());
            i += 2;
        } else if args[i] == "--duration" {
            duration_secs = Some(args.get(i + 1).expect("Usage: --duration <seconds>")
                .parse().expect("--duration must be a number"));
            i += 2;
        } else {
            positional.push(args[i].clone());
            i += 1;
        }
    }

    // Load RDS config
    let rds_config = match rds_config_path {
        Some(path) => rds_config::RdsConfig::from_file(&path).expect("Failed to load RDS config"),
        None => rds_config::RdsConfig::default(),
    };

    // Duration timer: spawn a thread that sets done after the specified time
    if let Some(secs) = duration_secs {
        let done_ref = done_sig.clone();
        std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_secs_f64(secs));
            eprintln!("Duration {:.1}s reached, shutting down...", secs);
            done_ref.store(true, atomic::Ordering::SeqCst);
        });
    }

    let audio_output = match wav_path {
        Some(path) => AudioOutput::Wav(path),
        None => AudioOutput::Playback,
    };

    let mut pos = positional.iter().map(|s| s.as_str());
    match pos.next() {
        Some("sigmf") => {
            let path = pos.next().expect("Usage: rradio sigmf <path> [tune_offset_khz]");
            let tune_offset: f32 = pos.next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0)
                * 1e3;
            let streamer = sigmf::SigmfStreamer::new(path).expect("Failed to open SigMF file");
            let source = IqSource::Sigmf { streamer, tune_offset };
            run(source, audio_output, done_sig, obs_settings, rds_debug, rds_metrics, record_path, rds_config.clone(), mpx_path.clone(), diag_path.clone());
        }
        Some("soapy") => {
            let filter = pos.next().expect("Usage: rradio soapy <filter> [station_mhz]");
            let station: f32 = pos.next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(96.1)
                * 1e6;
            let config = soapy::SoapyConfig {
                filter: filter.to_string(),
                station,
                bw: 600e6,
                fs: 2.4e6,
            };
            run(IqSource::Soapy { config }, audio_output, done_sig, obs_settings, rds_debug, rds_metrics, record_path, rds_config.clone(), mpx_path.clone(), diag_path.clone());
        }
        Some("pluto") => {
            let station: f32 = pos.next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(96.1)
                * 1e6;
            let config = pluto::SdrConfig {
                uri: "ip:pluto.local".to_string(),
                station,
                bw: 600e6,
                fs: 2.4e6,
            };
            run(IqSource::Pluto { config }, audio_output, done_sig, obs_settings, rds_debug, rds_metrics, record_path, rds_config.clone(), mpx_path.clone(), diag_path.clone());
        }
        _ => {
            eprintln!("Usage: rradio <source> [options] [--wav <output.wav>]");
            eprintln!("  rradio pluto [station_mhz]");
            eprintln!("  rradio soapy <filter> [station_mhz]");
            eprintln!("  rradio sigmf <path.sigmf-meta> [tune_offset_khz]");
            std::process::exit(1);
        }
    }
}