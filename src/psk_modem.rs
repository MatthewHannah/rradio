/// BPSK modem: hard decision + phase error extraction.
///
/// Makes a hard binary decision on a complex PSK symbol and returns the
/// phase error between the received symbol and the nearest constellation
/// point. The phase error is used for decision-directed carrier tracking.

use num_complex::Complex32;

pub struct BpskModem;

pub struct DemodResult {
    pub decision: bool,
    pub phase_error: f32,
}

impl BpskModem {
    pub fn new() -> Self {
        BpskModem
    }

    /// Demodulate one complex symbol.
    ///
    /// Returns the hard decision (true = +1, false = -1) and the phase error
    /// in radians between the received symbol and the nearest BPSK constellation
    /// point on the real axis.
    #[inline]
    pub fn demodulate(&self, symbol: Complex32) -> DemodResult {
        let decision = symbol.re >= 0.0;

        // Phase error relative to nearest constellation point:
        // If decision is +1 (re >= 0): error = atan2(im, re)
        // If decision is -1 (re < 0):  error = atan2(-im, -re)  (rotate to +1 reference)
        let phase_error = if decision {
            symbol.im.atan2(symbol.re)
        } else {
            (-symbol.im).atan2(-symbol.re)
        };

        DemodResult { decision, phase_error }
    }
}
