use rradio_dsp::biquad::Biquad;
use rradio_dsp::deemphasis::Deemphasis;
use rradio_dsp::filterable::Filter;
use rradio_dsp::pll::RealPll;
use num_complex::Complex32;

pub struct WidebandFmAudio<I> {
    input: I, downmix: StereoDownmixer,
    l_deemph: Deemphasis<f32>, r_deemph: Deemphasis<f32>,
    l_audio_filt: Biquad<f32>, r_audio_filt: Biquad<f32>,
}

#[derive(Debug, Copy, Clone)]
pub struct WidebandFmAudioOutput {
    pub left: f32,
    pub right: f32,
    pub stereo_lock: f32,
    pub mpx: f32,
}

impl<I> WidebandFmAudio<I> where I: Iterator<Item = f32> {
    fn new(fs: f32, i: I) -> WidebandFmAudio<I> {
        WidebandFmAudio {
            input: i, downmix: StereoDownmixer::new(fs),
            l_deemph: Deemphasis::new(fs, 75e-6), r_deemph: Deemphasis::new(fs, 75e-6),
            l_audio_filt: Biquad::lowpass(fs, 17000.0, 0.707),
            r_audio_filt: Biquad::lowpass(fs, 17000.0, 0.707),
        }
    }
    fn process(&mut self, s: f32) -> WidebandFmAudioOutput {
        let mono = s;
        let StereoDownmixerOutput { stereo, stereo_lock } = self.downmix.process(s);
        let l = self.l_audio_filt.process(self.l_deemph.process(mono + stereo));
        let r = self.r_audio_filt.process(self.r_deemph.process(mono - stereo));
        WidebandFmAudioOutput { left: l, right: r, stereo_lock, mpx: s }
    }
}
impl<I> Iterator for WidebandFmAudio<I> where I: Iterator<Item = f32> {
    type Item = WidebandFmAudioOutput;
    fn next(&mut self) -> Option<Self::Item> { let s = self.input.next()?; Some(self.process(s)) }
}
pub trait WidebandFmAudioIterable {
    fn wfm_audio(self, fs: f32) -> WidebandFmAudio<Self> where Self: Sized;
}
impl <I> WidebandFmAudioIterable for I where I: Iterator<Item = f32> {
    fn wfm_audio(self, fs: f32) -> WidebandFmAudio<Self> where Self: Sized { WidebandFmAudio::new(fs, self) }
}

struct StereoDownmixer {
    pll: RealPll<Biquad<f32>>,
}
struct StereoDownmixerOutput { stereo: f32, stereo_lock: f32 }

impl StereoDownmixer {
    fn new(fs: f32) -> Self {
        let pll_loop_filt = Biquad::lowpass(fs, 1000.0, 0.707);
        let pll = RealPll::new(38e3, fs, 0.05, pll_loop_filt, 2.0);
        StereoDownmixer { pll }
    }
    fn process(&mut self, s: f32) -> StereoDownmixerOutput {
        let pll_out = self.pll.process(s);
        StereoDownmixerOutput { stereo: pll_out.out * s, stereo_lock: pll_out.lock }
    }
}
