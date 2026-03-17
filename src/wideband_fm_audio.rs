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

#[derive(Debug, Copy, Clone)]
pub struct WidebandFmAudioOutput {
    pub left: f32,
    pub right: f32,
    pub rds: f32,
    pub stereo_lock: f32,
    pub rds_lock: f32,
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

    fn process(&mut self, s: f32) -> WidebandFmAudioOutput {
        let mono = s;
        let StereoDownmixerOutput { stereo, rds, stereo_lock, rds_lock } = self.downmix.process(s);

        let l = self
            .l_audio_filt
            .process(self.l_deemph.process(mono + stereo));
        let r = self
            .r_audio_filt
            .process(self.r_deemph.process(mono - stereo));

        WidebandFmAudioOutput {
            left: l,
            right: r,
            rds,
            stereo_lock,
            rds_lock,
        }
    }
}

impl<I> Iterator for WidebandFmAudio<I> where I: Iterator<Item = f32> {
    type Item = WidebandFmAudioOutput;

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
    pll_rds: RealPll<Biquad<f32>>,
    rds_filt: Biquad<f32>,
    rds_filt2: Biquad<f32>,
}

struct StereoDownmixerOutput {
    stereo: f32,
    rds: f32,
    stereo_lock: f32,
    rds_lock: f32,
}

impl StereoDownmixer {
    fn new(fs: f32) -> Self {
        let pll_loop_filt = Biquad::lowpass(fs, 1000.0, 0.707);
        let pll = RealPll::new(38e3, fs, 0.05, pll_loop_filt, 2.0);

        let pll_rds_loop_filt = Biquad::lowpass(fs, 1000.0, 0.707);
        let pll_rds = RealPll::new(57e3, fs, 0.05, pll_rds_loop_filt, 3.0);

        let rds_filt = Biquad::lowpass(fs, 3e3, 0.707);
        let rds_filt2 = rds_filt.clone();
        StereoDownmixer {
            pll,
            pll_rds,
            rds_filt,
            rds_filt2,
        }
    }

    fn process(&mut self, s: f32) -> StereoDownmixerOutput {
        // pilot extract
        let pll_out = self.pll.process(s);
        let rds_out = self.pll_rds.process(s);

        let rds = self.rds_filt.process(rds_out.out * s);
        let rds = self.rds_filt2.process(rds);

        StereoDownmixerOutput {
            stereo: pll_out.out * s,
            rds: rds,
            stereo_lock: pll_out.lock,
            rds_lock: rds_out.lock,
        }
    }
}
