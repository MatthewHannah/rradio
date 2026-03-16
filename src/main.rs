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

use std::sync::atomic;
use std::sync::Arc;

use plotly::{HeatMap, Plot, Scatter};
use num_complex::Complex32;
use rustfft::{num_traits::Zero, FftPlanner};

use crate::filterable::FilterableIter;
use crate::fm_demod::FmDemodulatable;
use crate::interleaver::InterleaveableIter;
use crate::pluto::SdrConfig;
use crate::resample::Downsampleable;
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

fn buffer_iq(done: &atomic::AtomicBool, config: &SdrConfig, mut out: buffer::SendBuf<Vec<Complex32>>) {
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

fn buffer_sigmf(done: &atomic::AtomicBool, mut streamer: sigmf::SigmfStreamer, mut out: buffer::SendBuf<Vec<Complex32>>, tune_offset: f32) {
    const CHUNK: usize = 4*1024*1024;
    let mut osc = osc::Osc::new(tune_offset, 6e6);
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

struct AudioPipelineSettings {
    iq_downsample: usize,
    fm_demod_downsample: usize,
    audio_downsample: usize,
}

fn audio_pipeline(done: &atomic::AtomicBool, fs: f32, inbuf: buffer::RecvBuf<Vec<Complex32>>, mut outbuf: buffer::SendBuf<Vec<f32>>, settings: AudioPipelineSettings, obs_settings: AudioPipelineObservationSettings) {
    let samples = buffer::RecvBufIter::new(inbuf);

    let fs_spy = fs;
    let samples = samples.maybe_spy(6000000, move |iq_samples| {
        println!("Got {} IQ samples for spectrogram", iq_samples.len());
        spectrogram(8192, 512, fs_spy, &iq_samples);
        plot_re_im(&iq_samples);
    }, obs_settings.spy_iq);

    // audio pipeline
    let filt1 = biquad::Biquad::lowpass(fs, 300000.0, 0.707);
    let filt2 = filt1.clone();
    let filtered = samples.dsp_filter(filt1).dsp_filter(filt2);
    let fs: f32 = fs / (settings.iq_downsample as f32);
    let resampled = filtered.downsample(settings.iq_downsample);

    let fs_spy = fs;
    let resampled = resampled.maybe_spy(6000000, move |audio_samples| {
        println!("Got {} audio samples for spectrogram", audio_samples.len());
        spectrogram(8192, 512, fs_spy, &audio_samples);
        plot_re_im(&audio_samples);
    }, obs_settings.spy_audio);

    let fm_filt = biquad::Biquad::lowpass(fs, 80000.0, 0.707);
    let fm_loop_filt = biquad::Biquad::lowpass(fs, 80000.0, 0.707);
    let demoded = resampled.dsp_filter(fm_filt).fm_demodulate(fm_loop_filt).downsample(settings.fm_demod_downsample);
    let fs = fs / (settings.fm_demod_downsample as f32);

    let fs_spy = fs;
    let demoded = demoded.maybe_spy(48000, move |audio_samples| {
        println!("Got {} audio samples for spectrogram", audio_samples.len());
        spectrogram(1024, 512, fs_spy, &audio_samples);
    }, obs_settings.spy_demoded);

    let mut lr = demoded.wfm_audio(fs).downsample(settings.audio_downsample).interleave();

    while !done.load(atomic::Ordering::SeqCst) {
        // buffering into outbuf
        let out = outbuf.get();
        if out.is_none() {
            break;
        }
        let mut out = out.unwrap();
        out.clear();
        out.extend((&mut lr).take(8192));

        if out.is_empty() {
            break;
        }
        outbuf.commit(out);
    }
}

fn run_from_sdr(station: f32, done_sig: Arc<atomic::AtomicBool>, obs_settings: AudioPipelineObservationSettings) {
    let config = SdrConfig {
        uri: "ip:192.168.2.1".to_string(),
        station,
        bw: 600e6,
        fs: 2.4e6,
    };
    let fs = config.fs;

    let settings = AudioPipelineSettings {
        iq_downsample: 2,
        fm_demod_downsample: 5,
        audio_downsample: 5,
    };
    if (fs / settings.iq_downsample as f32) / (settings.fm_demod_downsample as f32) / (settings.audio_downsample as f32) - 48000.0 > 100.0 {
        println!("Warning: output sample rate is {}, expected 48000", (fs / settings.iq_downsample as f32) / settings.fm_demod_downsample as f32 / settings.audio_downsample as f32);
    }

    let (iq_tx, iq_rx) = buffer::buf_pair(8);
    let (audio_tx, audio_rx) = buffer::buf_pair(8);

    let done_ref = done_sig.clone();
    let iq_config = config.clone();
    let iq_thread = std::thread::spawn(move || {
        buffer_iq(&done_ref, &iq_config, iq_tx);
    });

    let done_ref = done_sig.clone();
    let audio_thread = std::thread::spawn(move || {
        audio_pipeline(&done_ref, fs, iq_rx, audio_tx, settings, obs_settings);
    });

    let audio_rx_iter = buffer::RecvBufIter::new(audio_rx);
    playback::playback_iter(audio_rx_iter, 48000);

    iq_thread.join().unwrap();
    audio_thread.join().unwrap();
}

fn run_from_soapy(filter: &str, station: f32, done_sig: Arc<atomic::AtomicBool>, obs_settings: AudioPipelineObservationSettings) {
    let config = soapy::SoapyConfig {
        filter: filter.to_string(),
        station,
        bw: 600e6,
        fs: 1.2e6,
    };
    let fs = config.fs;

    let settings = AudioPipelineSettings {
        iq_downsample: 1,
        fm_demod_downsample: 5,
        audio_downsample: 5,
    };
    if (fs / settings.iq_downsample as f32) / (settings.fm_demod_downsample as f32) / (settings.audio_downsample as f32) - 48000.0 > 100.0 {
        println!("Warning: output sample rate is {}, expected 48000", (fs / settings.iq_downsample as f32) / settings.fm_demod_downsample as f32 / settings.audio_downsample as f32);
    }

    let (iq_tx, iq_rx) = buffer::buf_pair(8);
    let (audio_tx, audio_rx) = buffer::buf_pair(8);

    let done_ref = done_sig.clone();
    let iq_config = config.clone();
    let iq_thread = std::thread::spawn(move || {
        buffer_soapy(&done_ref, &iq_config, iq_tx);
    });

    let done_ref = done_sig.clone();
    let audio_thread = std::thread::spawn(move || {
        audio_pipeline(&done_ref, fs, iq_rx, audio_tx, settings, obs_settings);
    });

    let audio_rx_iter = buffer::RecvBufIter::new(audio_rx);
    playback::playback_iter(audio_rx_iter, 48000);

    iq_thread.join().unwrap();
    audio_thread.join().unwrap();
}

fn run_from_sigmf(path: &str, done_sig: Arc<atomic::AtomicBool>, obs_settings: AudioPipelineObservationSettings) {
    let streamer = sigmf::SigmfStreamer::new(path).expect("Failed to open SigMF file");

    let fs = 6e6;
    let settings = AudioPipelineSettings {
        iq_downsample: 2,
        fm_demod_downsample: 5,
        audio_downsample: 5,
    };
    if (fs / settings.iq_downsample as f32) / (settings.fm_demod_downsample as f32) / (settings.audio_downsample as f32) - 48000.0 > 100.0 {
        println!("Warning: output sample rate is {}, expected 48000", (fs / settings.iq_downsample as f32) / settings.fm_demod_downsample as f32 / settings.audio_downsample as f32);
    }

    let (iq_tx, iq_rx) = buffer::buf_pair(8);
    let (audio_tx, audio_rx) = buffer::buf_pair(8);

    let done_ref = done_sig.clone();
    let iq_thread = std::thread::spawn(move || {
        buffer_sigmf(&done_ref, streamer, iq_tx, 400e3);
    });

    let done_ref = done_sig.clone();
    let audio_thread = std::thread::spawn(move || {
        audio_pipeline(&done_ref, fs, iq_rx, audio_tx, settings, obs_settings);
    });

    let audio_rx_iter = buffer::RecvBufIter::new(audio_rx);
    // let audio_samples: Vec<f32> = audio_rx_iter.take(48000 * 10).collect();
    // playback::playback_buffer(audio_samples, 48000);
    playback::playback_iter(audio_rx_iter, 48000);
    done_sig.store(true, atomic::Ordering::SeqCst);

    iq_thread.join().unwrap();
    audio_thread.join().unwrap();
}

fn main() {
    let mut args = std::env::args().skip(1);

    let done_sig = Arc::new(atomic::AtomicBool::new(false));

    let done_ref = done_sig.clone();
    ctrlc::set_handler(move || {
        done_ref.store(true, atomic::Ordering::SeqCst);
    }).expect("Error setting CTRL-C handler");

    let settings = AudioPipelineObservationSettings {
        spy_iq: false,
        spy_demoded: false,
        spy_audio: false,
    };

    match args.next().as_deref() {
        Some("sigmf") => {
            let path = args.next().expect("Usage: rradio sigmf <path>");
            run_from_sigmf(&path, done_sig, settings);
        }
        Some("soapy") => {
            let filter = args.next().expect("Usage: rradio soapy <filter> [station_mhz]");
            let default_station = 96.1;
            let station: f32 = args.next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(default_station)
                * 1e6;
            run_from_soapy(&filter, station, done_sig, settings);
        }
        station_arg => {
            let default_station = 96.1;
            let station: f32 = station_arg
                .and_then(|s| s.parse().ok())
                .unwrap_or(default_station)
                * 1e6;
            run_from_sdr(station, done_sig, settings);
        }
    }
}