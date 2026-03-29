/// Biphase (Manchester) decoder — port of Redsea's BiphaseDecoder.
///
/// Takes PSK symbols at chip rate (2375 Hz) and outputs data bits at
/// bit rate (1187.5 Hz). A biphase symbol consists of two PSK symbols
/// with opposite phase. The decoder tracks the even/odd alignment by
/// comparing energy sums over a 128-symbol window.

use num_complex::Complex32;

pub struct BiphaseDecoder {
    prev_psk_symbol: Complex32,
    clock_history: [f32; 128],
    clock: usize,
    clock_polarity: usize,

    // Diagnostics: last computed energy sums and current state
    pub last_even_sum: f32,
    pub last_odd_sum: f32,
    pub polarity: usize,
    pub biphase_energy: f32,  // current chip's |biphase.re|
}

/// Result from pushing a PSK symbol — may or may not produce a bit.
pub struct BiphaseResult {
    pub bit: bool,
    pub has_value: bool,
}

impl BiphaseDecoder {
    pub fn new() -> Self {
        BiphaseDecoder {
            prev_psk_symbol: Complex32::new(0.0, 0.0),
            clock_history: [0.0; 128],
            clock: 0,
            clock_polarity: 0,
            last_even_sum: 0.0,
            last_odd_sum: 0.0,
            polarity: 0,
            biphase_energy: 0.0,
        }
    }

    /// Push one PSK symbol (at 2375 Hz). Returns a bit when available (at 1187.5 Hz).
    pub fn push(&mut self, psk_symbol: Complex32) -> BiphaseResult {
        // A biphase symbol = difference of consecutive PSK symbols
        let biphase_symbol = (psk_symbol - self.prev_psk_symbol) * 0.5;
        let bit = biphase_symbol.re >= 0.0;
        let has_value = (self.clock % 2) == self.clock_polarity;
        self.prev_psk_symbol = psk_symbol;
        self.biphase_energy = biphase_symbol.re.abs();

        // Track which alignment (even/odd) has more energy
        self.clock_history[self.clock] = self.biphase_energy;
        self.clock += 1;

        // Every 128 PSK symbols, check whether to shift the decoding window
        if self.clock == 128 {
            let mut even_sum = 0.0_f32;
            let mut odd_sum = 0.0_f32;
            for i in (0..128).step_by(2) {
                even_sum += self.clock_history[i];
                odd_sum += self.clock_history[i + 1];
            }

            if even_sum > odd_sum {
                self.clock_polarity = 0;
            } else if odd_sum > even_sum {
                self.clock_polarity = 1;
            }
            self.last_even_sum = even_sum;
            self.last_odd_sum = odd_sum;
            self.polarity = self.clock_polarity;

            self.clock_history = [0.0; 128];
            self.clock = 0;
        }

        BiphaseResult { bit, has_value }
    }
}

/// Delta (differential) decoder: XOR of consecutive bits.
pub struct DeltaDecoder {
    prev: bool,
}

impl DeltaDecoder {
    pub fn new() -> Self {
        DeltaDecoder { prev: false }
    }

    pub fn decode(&mut self, input: bool) -> bool {
        let output = input != self.prev;
        self.prev = input;
        output
    }
}
