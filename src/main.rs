mod sigmf;
mod osc;
mod biquad;
mod resample;
mod fm_demod;
mod playback;
mod deemphasis;
mod filterable;
mod spy;

use plotly::{HeatMap, Plot, Scatter};
use plotly::common::Mode;
use num_complex::Complex32;
use rustfft::{num_traits::Zero, FftPlanner};

use crate::filterable::FilterableIter;
use crate::fm_demod::FmDemodulatable;
use crate::resample::Upsampleable;
use crate::resample::Downsampleable;
use crate::spy::SpyableIter;

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

fn main() {
    let fs = 6e6f32;

    let up = 1;
    let down = 5;

    let file = sigmf::SigmfStreamer::new("./res/fm_radio_20250920_6msps.sigmf-data").unwrap();

    let bw = 6e6;

    // select channel
    let center_freq = 94.5e6;
    let station = 96.1e6;
    let freq_shift: f32 = center_freq - station;

    if freq_shift.abs() > ((bw / 2.0) - 200e3) {
        panic!("Station {} MHz is out of the capture bandwidth {} MHz centered at {} MHz", station / 1e6, bw / 1e6, center_freq / 1e6);
    }

    let o = osc::Osc::new(freq_shift, fs);
    let shifted = file.zip(o.into_iter()).map(|(s, w)| s * w);

    // filter and resample to 1.2 Msps
    let filt1 = biquad::Biquad::new(0.02008282, 0.04016564, 0.02008282, -1.56097580, 0.64130708);
    let filt2 = filt1.clone();
    let filtered = shifted.dsp_filter(filt1).dsp_filter(filt2);
    let fs = fs * (up as f32) / (down as f32);
    let resampled = filtered.upsample(up).downsample(down);

    // demodulate FM
    let fm_filt = biquad::Biquad::new(0.03357068, 0.06714135, 0.03357068, -1.41893478, 0.55321749);
    let fm_loop_filt = biquad::Biquad::new(0.03357068, 0.06714135, 0.03357068, -1.41893478, 0.55321749);
    let demoded = resampled.dsp_filter(fm_filt).fm_demodulate(fm_loop_filt).downsample(down);
    let fs = fs / (down as f32);
    //spectrogram(8192, 256, fs, demoded.take(1000000).map(|x| Complex32::new(x, 0.0)).collect::<Vec<Complex32>>().as_slice());

    // pilot recovery
    // filter for pilot (maybe unneeded)
    // 240kHz sample, 19kHz center, Q = 8
    let pilot_filt = biquad::Biquad::new(0.02895880, 0.0, -0.02895880, -1.70673525, 0.94208240);
    let pilot = demoded.dsp_filter(pilot_filt);

    spectrogram(8192, 256, fs, pilot.take(1000000).collect::<Vec<f32>>().as_slice());

    // spin up pll on pilot
    // multiple output by itself to get 38 kHz
    // lowpass to get 38 ref
    // multiply by 38 kHz ref to get stereo difference
    // lowpass to get stereo difference

    //play_mono(demoded, fs);
}