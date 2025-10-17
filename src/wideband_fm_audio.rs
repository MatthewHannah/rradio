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
        if (fs - 240000.0).abs() > 1.0 {
            panic!("Hardcoded filter parameters for 240kHz sample rate");
        }

        WidebandFmAudio {
            input: i,
            downmix: StereoDownmixer::new(fs),
            l_deemph: Deemphasis::new(fs, 75e-6),
            r_deemph: Deemphasis::new(fs, 75e-6),
            l_audio_filt: Biquad::new(0.04125202, 0.08250404, 0.04125202, -1.34891824, 0.51392633),
            r_audio_filt: Biquad::new(0.04125202, 0.08250404, 0.04125202, -1.34891824, 0.51392633),
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
        // in place downsample
        let init = self.input.next()?;
        let val = self.process(init);
        for _ in 0..4 {
            let next = self.input.next();
            if next.is_some() {
                self.process(next?);
            }
        }
        Some(val)
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
    lp_filt1: Biquad<f32>,
    lp_filt2: Biquad<f32>,
}

impl StereoDownmixer {
    fn new(fs: f32) -> Self {
        // 1kHz cutoff, 240kHz sample, butterworth
        let pll_loop_filt = Biquad::new(
            1.68223247e-4,
            3.36446494e-4,
            1.68223247e-4,
            -1.96297470,
            0.96364759,
        );
        let pll = RealPll::new(19e3, fs, 1.0, pll_loop_filt);

        // 23kHz cutoff, 240kHz sample, butterworth
        let lp_filt1 = Biquad::new(0.65120842, -1.30241684, 0.65120842, -1.17684383, 0.42798985);
        let lp_filt2 = lp_filt1.clone();

        StereoDownmixer {
            pll,
            lp_filt1,
            lp_filt2,
        }
    }

    fn process(&mut self, s: f32) -> f32 {
        let tone = self.pll.process(s);
        // square tone, then lowpass to get clean 38 kHz ref
        let mix = self.lp_filt1.process(tone * tone);

        // lp incoming to kill mono
        // mix with incoming to get stereo diff
        self.lp_filt2.process(s) * mix
    }
}

