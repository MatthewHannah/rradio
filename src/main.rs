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
mod agc;
mod costas;

use std::sync::atomic;
use std::sync::Arc;

use plotly::{HeatMap, Plot, Scatter};
use num_complex::Complex32;
use rustfft::{num_traits::Zero, FftPlanner};

use crate::costas::CostasDemodulable;
use crate::filterable::FilterableIter;
use crate::fm_demod::FmDemodulatable;
use crate::gardner_clock_recovery::GardnerClockRecoverable;
use crate::interleaver::InterleaveableIter;
use crate::manchester::ManchesterDecodable;
use crate::rds_block_sync::RdsBlockSyncable;
use crate::resample::{Downsampleable, Upsampleable};
use crate::spy::SpyableIter;
use crate::wideband_fm_audio::WidebandFmAudioIterable;

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
    let mut osc = osc::Osc::new(tune_offset, fs);
    while !done.load(atomic::Ordering::SeqCst) {
        if let Some(mut buf) = out.get() {
            buf.clear();
            buf.extend((&mut streamer).take(CHUNK).map(|s| s * osc.next()));
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
    mut rds_out: buffer::SendBuf<Vec<Complex32>>,
    settings: SignalPipelineSettings,
    obs_settings: AudioPipelineObservationSettings,
) {
    let samples = buffer::RecvBufIter::new(inbuf);

    let fs_spy = fs;
    let samples = samples.maybe_spy(6000000, move |iq_samples| {
        spectrogram(8192, 512, fs_spy, &iq_samples);
    }, obs_settings.spy_iq);

    // IQ filtering + downsample
    let filt1 = biquad::Biquad::lowpass(fs, 300000.0, 0.707);
    let filt2 = filt1.clone();
    let filtered = samples.dsp_filter(filt1).dsp_filter(filt2);
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

    // Wideband FM audio (stereo + RDS extraction) — tee to two consumers
    let mut wfm = demoded.wfm_audio(_fs);

    let mut audio_batch: Option<buffer::BufToken<Vec<(f32, f32)>>> = None;
    let mut rds_batch: Option<buffer::BufToken<Vec<Complex32>>> = None;

    while !done.load(atomic::Ordering::SeqCst) {
        // Ensure we have output buffers
        if audio_batch.is_none() {
            audio_batch = match audio_out.get() {
                Some(mut tok) => { tok.clear(); Some(tok) }
                None => break,
            };
        }
        if rds_batch.is_none() {
            rds_batch = rds_out.try_get().map(|mut tok| { tok.clear(); tok });
        }

        let sample = match wfm.next() {
            Some(s) => s,
            None => break,
        };

        // Tee: audio gets (left, right), RDS gets baseband
        if let Some(ref mut buf) = audio_batch {
            buf.push((sample.left, sample.right));
            if buf.len() >= 4096 {
                audio_out.commit(audio_batch.take().unwrap());
            }
        }
        if let Some(ref mut buf) = rds_batch {
            buf.push(sample.rds);
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

const AUDIO_DOWNSAMPLE: usize = 5;
const RDS_FIRST_DECIMATE: usize = 12; // 240 kHz → 20 kHz
const RDS_UPSAMPLE: usize = 19;       // 20 kHz → 380 kHz (virtual)
const RDS_SECOND_DECIMATE: usize = 20; // 380 kHz → 19 kHz
const RDS_FS: f32 = 19000.0;          // final RDS sample rate

/// RDS matched filter taps: RRC with α=1 at 2375 Hz chip rate.
/// Generated by py/rds_matched_filter.py --fs 19000 --spans 8 --window blackman
const RDS_MATCHED_FILTER_TAPS: [f32; 129] = [
    6.1299193e-21, -6.9963545e-08, 2.3289892e-20, 6.7633324e-07,
    1.7681212e-06, 2.0356788e-06, 2.1531523e-19, -4.3595446e-06,
    -8.4422146e-06, -7.9353035e-06, -7.2503226e-19, 1.3144058e-05,
    2.3349695e-05, 2.0482121e-05, -1.5158848e-18, -3.0587236e-05,
    -5.2217386e-05, -4.4271209e-05, 3.1813753e-18, 6.2561361e-05,
    1.0440662e-04, 8.6754307e-05, 5.1026858e-18, -1.1848720e-04,
    -1.9489120e-04, -1.5983412e-04, -4.5568355e-34, 2.1343907e-04,
    3.4770565e-04, 2.8270348e-04, 0.0000000e+00, -3.7205590e-04,
    -6.0248867e-04, -4.8734863e-04, 0.0000000e+00, 6.3645518e-04,
    1.0280190e-03, 8.3019361e-04, 0.0000000e+00, -1.0837966e-03,
    -1.7530309e-03, -1.4193314e-03, 1.3369382e-33, 1.8698438e-03,
    3.0451715e-03, 2.4867866e-03, 1.5782676e-33, -3.3540995e-03,
    -5.5482995e-03, -4.6167697e-03, 0.0000000e+00, 6.5439316e-03,
    1.1184980e-02, 9.6846800e-03, 0.0000000e+00, -1.5316045e-02,
    -2.8272031e-02, -2.7055719e-02, 0.0000000e+00, 5.9373683e-02,
    1.4825785e-01, 2.5335044e-01, 3.5349513e-01, 4.2560313e-01,
    4.5186651e-01, 4.2560313e-01, 3.5349513e-01, 2.5335044e-01,
    1.4825785e-01, 5.9373683e-02, 0.0000000e+00, -2.7055719e-02,
    -2.8272031e-02, -1.5316045e-02, 0.0000000e+00, 9.6846800e-03,
    1.1184980e-02, 6.5439316e-03, 0.0000000e+00, -4.6167697e-03,
    -5.5482995e-03, -3.3540995e-03, 1.5782676e-33, 2.4867866e-03,
    3.0451715e-03, 1.8698438e-03, 1.3369382e-33, -1.4193314e-03,
    -1.7530309e-03, -1.0837966e-03, 0.0000000e+00, 8.3019361e-04,
    1.0280190e-03, 6.3645518e-04, 0.0000000e+00, -4.8734863e-04,
    -6.0248867e-04, -3.7205590e-04, 0.0000000e+00, 2.8270348e-04,
    3.4770565e-04, 2.1343907e-04, -4.5568355e-34, -1.5983412e-04,
    -1.9489120e-04, -1.1848720e-04, 5.1026858e-18, 8.6754307e-05,
    1.0440662e-04, 6.2561361e-05, 3.1813753e-18, -4.4271209e-05,
    -5.2217386e-05, -3.0587236e-05, -1.5158848e-18, 2.0482121e-05,
    2.3349695e-05, 1.3144058e-05, -7.2503226e-19, -7.9353035e-06,
    -8.4422146e-06, -4.3595446e-06, 2.1531523e-19, 2.0356788e-06,
    1.7681212e-06, 6.7633324e-07, 2.3289892e-20, -6.9963545e-08,
    6.1299193e-21,
];

fn rds_pipeline(done: &atomic::AtomicBool, rds_rx: buffer::RecvBuf<Vec<Complex32>>, wfm_fs: f32) {
    let chip_rate = 2375.0_f32;
    let rds_iter = buffer::RecvBufIter::new(rds_rx);

    // Anti-alias filter before ÷12 decimation (complex), cutoff at 9.5 kHz
    // (Nyquist of the 19 kHz target rate)
    let aa1: biquad::Biquad<Complex32> = biquad::Biquad::lowpass(wfm_fs, 9000.0, 0.707);
    let aa2 = aa1.clone();
    let aa3 = aa1.clone();

    // Rational resample: ÷12 → 20 kHz → ↑19 ↓20 → 19 kHz
    // The anti-alias filter for the ↑19↓20 stage runs at the virtual 380 kHz rate.
    // We use a biquad at 9 kHz (below 9.5 kHz Nyquist) at the 20 kHz intermediate rate.
    let intermediate_fs = wfm_fs / RDS_FIRST_DECIMATE as f32; // 20 kHz
    let resamp_aa1: biquad::Biquad<Complex32> = biquad::Biquad::lowpass(intermediate_fs * RDS_UPSAMPLE as f32, 9000.0, 0.707);
    let resamp_aa2 = resamp_aa1.clone();

    let matched_filter: fir::Fir<Complex32> = fir::Fir::new(RDS_MATCHED_FILTER_TAPS.to_vec());
    let mut groups = rds_iter
        .dsp_filter(aa1).dsp_filter(aa2).dsp_filter(aa3)
        .downsample(RDS_FIRST_DECIMATE)
        .upsample(RDS_UPSAMPLE)
        .dsp_filter(resamp_aa1).dsp_filter(resamp_aa2)
        .downsample(RDS_SECOND_DECIMATE)
        .dsp_filter(matched_filter)
        .costas_demod(0.05)
        .dsp_filter(agc::Agc::new(RDS_FS, 10.0, 0.001))
        .clock_recover(chip_rate, RDS_FS)
        .manchester_decode()
        .rds_block_sync();

    let mut decoder = rds_block_sync::RdsDecoder::new();

    while !done.load(atomic::Ordering::SeqCst) {
        match groups.next() {
            Some(group) => {
                decoder.process(&group);
            }
            None => break,
        }
    }
}

fn run(iq_source: IqSource, audio_output: AudioOutput, done_sig: Arc<atomic::AtomicBool>, obs_settings: AudioPipelineObservationSettings) {
    let fs = match &iq_source {
        IqSource::Pluto { config } => config.fs,
        IqSource::Soapy { config } => config.fs,
        IqSource::Sigmf { streamer, .. } => streamer.sample_rate(),
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
    let (iq_tx, iq_rx) = buffer::buf_pair(8);
    let (audio_tx, audio_rx) = buffer::buf_pair::<Vec<(f32, f32)>>(8);
    let (rds_tx, rds_rx) = buffer::buf_pair::<Vec<Complex32>>(8);

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
        signal_pipeline(&done_ref, fs, iq_rx, audio_tx, rds_tx, settings, obs_settings);
    });

    // Thread 3: RDS consumer
    let done_ref = done_sig.clone();
    let rds_thread = std::thread::spawn(move || {
        rds_pipeline(&done_ref, rds_rx, wfm_fs);
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
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--wav" {
            wav_path = Some(args.get(i + 1).expect("Usage: --wav <path>").clone());
            i += 2;
        } else {
            positional.push(args[i].clone());
            i += 1;
        }
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
            run(source, audio_output, done_sig, obs_settings);
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
            run(IqSource::Soapy { config }, audio_output, done_sig, obs_settings);
        }
        Some("pluto") => {
            let station: f32 = pos.next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(96.1)
                * 1e6;
            let config = pluto::SdrConfig {
                uri: "ip:192.168.2.1".to_string(),
                station,
                bw: 600e6,
                fs: 2.4e6,
            };
            run(IqSource::Pluto { config }, audio_output, done_sig, obs_settings);
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