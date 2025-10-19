use std::time::Duration;
use std::sync::{Arc, Mutex, Condvar};

use cpal::{OutputCallbackInfo, SampleFormat, StreamError};
use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};


fn playback<D,E,S>(fs: u32, data_callback: D, error_callback: E, stop_callback: S) where
    D: FnMut(&mut [f32], &OutputCallbackInfo) + Send + 'static,
    E: FnMut(StreamError) + Send + 'static,
    S: FnOnce(),
{
    let host = cpal::default_host();
    let device = host.default_output_device().expect("no output device available");
    let mut supported_configs_range = device.supported_output_configs()
        .expect("error while querying configs");

    let supported_config = supported_configs_range.find(|config| {
        config.max_sample_rate().0 >= fs && config.min_sample_rate().0 <= fs && config.sample_format() == SampleFormat::F32 && config.channels() == 2
    }).expect("no supported config with required sample rate").with_sample_rate(cpal::SampleRate(fs));

    println!("Output config: {:?}", supported_config);

    let sample_format = supported_config.sample_format();
    let config = supported_config.into();

    let stream = match sample_format {
        SampleFormat::F32 => device.build_output_stream(&config, data_callback, error_callback, None),
        sample_format => panic!("Unsupported sample format '{sample_format}'")
    }.unwrap();

    stream.play().unwrap();
    stop_callback();
}

#[allow(dead_code)]
pub fn playback_sine(freq: f32, fs: u32, duration: Duration) {
    let mut playable = PlayableOsc::new(freq, fs as f32);
    let data_fn = move |data: &mut [f32], _: &cpal::OutputCallbackInfo| {
        playable.write_samples(data);
    };
    let err_fn = |err| eprintln!("an error occurred on the output audio stream: {}", err);
    let stop_fn = || std::thread::sleep(duration);
    playback(fs, data_fn, err_fn, stop_fn);
}

pub fn playback_iter<I>(mut samples: I, fs: u32) where
    I: Iterator<Item=f32> + Send + 'static
{
    let cv = Condvar::new();
    let done = Mutex::new(false);
    let data_done_signal =Arc::new((cv, done));
    let wait_done_signal = data_done_signal.clone();

    let data_fn = move |data: &mut [f32], _: &cpal::OutputCallbackInfo| {
        let mut done = false;
        for sample in data.iter_mut() {
            *sample = match samples.next() {
                Some(s) => s,
                None => {
                    done = true;
                    0.0
                }
            }
        }

        if done {
            let (cv, lock) = &*data_done_signal;
            let mut done = lock.lock().unwrap();
            *done = true;
            cv.notify_all();
        }
    };
    let err_fn = |err| eprintln!("an error occurred on the output audio stream: {}", err);
    let stop_fn = || {
        let (cv, lock) = &*wait_done_signal;
        let guard = lock.lock().unwrap();
        drop(cv.wait_while(guard, |d| !*d).unwrap());
    };
    playback(fs, data_fn, err_fn, stop_fn);
}

pub fn playback_buffer(samples: Vec<f32>, fs: u32) {
    playback_iter(samples.into_iter(), fs);
}

#[derive(Clone)]
struct PlayableOsc {
    osc: crate::osc::Osc,
}

impl PlayableOsc {
    fn new(freq: f32, sample_rate: f32) -> Self {
        PlayableOsc {
            osc: crate::osc::Osc::new(freq, sample_rate),
        }
    }

    fn next_sample(&mut self) -> f32 {
        let c = self.osc.next();
        c.re
    }

    fn write_samples(&mut self, data: &mut [f32]) {
        let mut l = true;
        let mut curr = 0.0;
        for sample in data.iter_mut() {
            if l {
                curr = self.next_sample();
            }
            *sample = curr;
            l = !l;
        }
    }
}
