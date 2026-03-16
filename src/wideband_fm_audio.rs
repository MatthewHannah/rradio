use crate::biquad::Biquad;
use crate::deemphasis::Deemphasis;
use crate::filterable::Filter;
use crate::pll::RealPll;

pub struct WidebandFmAudio<I> {
    input: I,
    downmix: StereoDownmixer,
    l_deemph: Deemphasis<f32>,
    r_deemph: Deemphasis<f32>,
    l_audio_filt: Biquad<f32>,
    r_audio_filt: Biquad<f32>,
}

impl<I> WidebandFmAudio<I> where I: Iterator<Item = f32> {
    fn new(fs: f32, i: I) -> WidebandFmAudio<I> {
        WidebandFmAudio {
            input: i,
            downmix: StereoDownmixer::new(fs),
            l_deemph: Deemphasis::new(fs, 75e-6),
            r_deemph: Deemphasis::new(fs, 75e-6),
            l_audio_filt: Biquad::lowpass(fs, 17000.0, 0.707),
            r_audio_filt: Biquad::lowpass(fs, 17000.0, 0.707),
        }
    }

    fn process(&mut self, s: f32) -> (f32, f32) {
        let mono = s;
        let stereo = self.downmix.process(s);

        let l = self
            .l_audio_filt
            .process(self.l_deemph.process(mono + stereo));
        let r = self
            .r_audio_filt
            .process(self.r_deemph.process(mono - stereo));

        (l, r)
    }
}

impl<I> Iterator for WidebandFmAudio<I> where I: Iterator<Item = f32> {
    type Item = (f32, f32);

    fn next(&mut self) -> Option<Self::Item> {
        let s = self.input.next()?;
        Some(self.process(s))
    }
}

pub trait WidebandFmAudioIterable {
    fn wfm_audio(self, fs: f32) -> WidebandFmAudio<Self> where Self: Sized;
}

impl <I> WidebandFmAudioIterable for I where I: Iterator<Item = f32> {
    fn wfm_audio(self, fs: f32) -> WidebandFmAudio<Self> where Self: Sized {
        WidebandFmAudio::new(fs, self)
    }
}

struct StereoDownmixer {
    pll: RealPll<Biquad<f32>>,
    hp_filt: Biquad<f32>,
}

impl StereoDownmixer {
    fn new(fs: f32) -> Self {
        let pll_loop_filt = Biquad::lowpass(fs, 1000.0, 0.707);
        let pll = RealPll::new(19e3, fs, 0.1, pll_loop_filt);

        let hp_filt = Biquad::highpass(fs, 23000.0, 0.707);

        StereoDownmixer {
            pll,
            hp_filt,
        }
    }

    fn process(&mut self, s: f32) -> f32 {
        // pilot extract
        let tone = self.pll.process(s);

        // tone = cos(w)
        // tone^2 = (1 + cos(2w)) / 2
        // real tone^2 has DC that we have to kill
        // otherwise we pass mono combined with stereo downmix
        // mix with incoming to get stereo diff
        // multiply by 2 to get back to unity after ^2
        self.hp_filt.process(tone * tone * 2.0) * s
    }
}
