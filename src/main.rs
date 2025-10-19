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

use std::sync::atomic;
use std::sync::Arc;

use plotly::{HeatMap, Plot, Scatter};
use plotly::common::Mode;
use num_complex::Complex32;
use rustfft::{num_traits::Zero, FftPlanner};

use crate::filterable::{Filter, FilterableIter};
use crate::fm_demod::FmDemodulatable;
use crate::interleaver::InterleaveableIter;
use crate::pluto::PlutoSdrIqStreamer;
use crate::resample::Upsampleable;
use crate::resample::Downsampleable;
use crate::spy::SpyableIter;
use crate::wideband_fm_audio::WidebandFmAudioIterable;

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

fn play_mono<I>(demoded: I, fs: f32) where I: Iterator<Item = f32> + Send + 'static {
    // grab mono and downsample to 48 KHz
    let mono_filt = biquad::Biquad::new(0.04125202, 0.08250404, 0.04125202, -1.34891824, 0.51392633);
    let mono_filt_stage2 = mono_filt.clone();
    let pilot_notch = biquad::Biquad::new(0.69904438, 1.10917838, 0.69904438, 1.10917838, 0.39808875);
    let mono = demoded.dsp_filter(mono_filt).dsp_filter(mono_filt_stage2).downsample(5).dsp_filter(pilot_notch);

    let fs = fs / 5.0;

    let deemph = deemphasis::Deemphasis::new(fs, 75e-6);
    let mono = mono.dsp_filter(deemph);

    let stereo_mono = mono.upsample(2);

    let now = std::time::Instant::now();
    let seconds_to_play = 5;
    let mono_audio = stereo_mono.take((fs as usize) * seconds_to_play * 2).collect::<Vec<f32>>();
    println!("Took {:?} to do {} of audio", now.elapsed(), seconds_to_play);

    playback::playback_buffer(mono_audio, fs as u32);
    //playback::playback_iter(stereo_mono, fs as u32);
}

fn buffer_iq(done: &atomic::AtomicBool, streamer: &mut PlutoSdrIqStreamer, mut out: buffer::SendBuf<Vec<Complex32>>) {
    while !done.load(atomic::Ordering::SeqCst) {
        if let Some(mut buf) = out.get() {
            if let Ok(()) = streamer.collect_iq(&mut buf) {
                out.commit(buf);
            }
        }
    }
}

fn audio_pipeline(done: &atomic::AtomicBool, fs: f32, mut inbuf: buffer::RecvBuf<Vec<Complex32>>, mut outbuf: buffer::SendBuf<Vec<f32>>) {
    let down = 5;
    let samples = buffer::RecvBufIter::new(inbuf);

    // audio pipeline
    let filt1 = biquad::Biquad::new(0.02008282, 0.04016564, 0.02008282, -1.56097580, 0.64130708);
    let filt2 = filt1.clone();
    let filtered = samples.dsp_filter(filt1).dsp_filter(filt2);
    let fs: f32 = fs / (down as f32);
    let resampled = filtered.downsample(down);

    let fm_filt = biquad::Biquad::new(0.03357068, 0.06714135, 0.03357068, -1.41893478, 0.55321749);
    let fm_loop_filt = biquad::Biquad::new(0.03357068, 0.06714135, 0.03357068, -1.41893478, 0.55321749);
    let demoded = resampled.dsp_filter(fm_filt).fm_demodulate(fm_loop_filt).downsample(down);
    let fs = fs / (down as f32);

    let mut lr = demoded.wfm_audio(fs).interleave();
    let fs = fs / 5.0;

    while !done.load(atomic::Ordering::SeqCst) {
        // buffering into outbuf
        let mut out = outbuf.get().unwrap();
        out.resize(8192, 0.0);
        (&mut lr).take(8192).zip(out.iter_mut()).for_each(|(s, o)| *o = s);
        outbuf.commit(out);
    }
}

fn main() {

    let default_station = 96.1;
    let station: f32 = std::env::args().skip(1).next().map(|x| x.parse()).unwrap_or(Ok(default_station)).unwrap_or(default_station) * 1e6;    

    let fs = 6e6f32;

    let ctx = industrial_io::Context::from_uri("ip:192.168.2.1").unwrap();
    let mut sdr = pluto::PlutoSdr::new(ctx).unwrap();
    let bw: f32 = 6e6;

    sdr.set_center(station).unwrap();
    sdr.set_rf_bandwidth(bw).unwrap();
    sdr.set_sampling_freq(fs).unwrap();
       
    let done_sig = Arc::new(atomic::AtomicBool::new(false));

    let (iq_tx, iq_rx) = buffer::buf_pair(8);
    let (audio_tx, audio_rx) = buffer::buf_pair(8);

    let done_ref = done_sig.clone();
    ctrlc::set_handler(move || {
        done_ref.store(true, atomic::Ordering::SeqCst);
    }).expect("Error setting CTRL-C handler");

    let done_ref = done_sig.clone();
    let iq_thread = std::thread::spawn(move || {
        let mut iq_streamer= sdr.start_iq().unwrap();
        buffer_iq(&done_ref, &mut iq_streamer, iq_tx);
        sdr.stop_iq(iq_streamer);
    });

    let done_ref = done_sig.clone();
    let audio_thread = std::thread::spawn(move || {
        audio_pipeline(&done_ref, fs, iq_rx, audio_tx);
    });

    let audio_rx_iter = buffer::RecvBufIter::new(audio_rx);

    let fs = 48000.0;
    // let now = std::time::Instant::now();
    // let seconds_to_play = 10;
    // let audio = audio_rx_iter.take((fs as usize) * seconds_to_play).collect::<Vec<f32>>();
    // println!("Audio len {}", audio.len());
    // println!("Took {:?} to do {} of audio", now.elapsed(), seconds_to_play);

    playback::playback_iter(audio_rx_iter, fs as u32);

    iq_thread.join().unwrap();
    audio_thread.join().unwrap();

    //playback_iter(lr, fs as u32);
}