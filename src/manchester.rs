/// Developed using Claude Opus 4.6
/// Manchester (biphase) decoder with dual-alignment tracking
/// and differential decoding for RDS.
///
/// Runs both chip alignments (even/odd) in parallel, continuously
/// tracking the valid-pair ratio for each. Emits bits from whichever
/// alignment currently has the better score.

const SCORE_WINDOW: usize = 50; // pairs to track for scoring

/// One alignment's pairing state
struct AlignmentState {
    pending: Option<i8>,   // buffered first chip of current pair
    prev_diff_bit: u8,     // differential decode state
    // Sliding window score tracking
    history: Vec<bool>,    // ring buffer: true = valid pair
    history_idx: usize,
    valid_count: usize,
}

impl AlignmentState {
    fn new() -> Self {
        AlignmentState {
            pending: None,
            prev_diff_bit: 0,
            history: vec![false; SCORE_WINDOW],
            history_idx: 0,
            valid_count: 0,
        }
    }

    /// Feed a chip. Returns Some(data_bit) when a pair completes validly.
    fn push_chip(&mut self, chip: i8) -> Option<u8> {
        match self.pending.take() {
            None => {
                self.pending = Some(chip);
                None
            }
            Some(c0) => {
                let valid = (c0 == 1 && chip == -1) || (c0 == -1 && chip == 1);

                // Update sliding window score
                let old = self.history[self.history_idx];
                if old { self.valid_count -= 1; }
                self.history[self.history_idx] = valid;
                if valid { self.valid_count += 1; }
                self.history_idx = (self.history_idx + 1) % SCORE_WINDOW;

                if valid {
                    let diff_bit = if c0 == -1 && chip == 1 { 1 } else { 0 };
                    let data_bit = diff_bit ^ self.prev_diff_bit;
                    self.prev_diff_bit = diff_bit;
                    Some(data_bit)
                } else {
                    None
                }
            }
        }
    }

    fn score(&self) -> usize {
        self.valid_count
    }
}

pub struct ManchesterDecoder<I: Iterator<Item = f32>> {
    iter: I,
    align: [AlignmentState; 2],
    // Buffer: alignment 1 is one chip behind alignment 0, so when
    // alignment 1 produces a bit, alignment 0 may also have just
    // produced one. We buffer output from the active alignment.
    output_buf: Option<u8>,
    total_chips: usize,
}

fn threshold(sample: f32) -> i8 {
    if sample > 0.0 { 1 } else { -1 }
}

impl<I: Iterator<Item = f32>> ManchesterDecoder<I> {
    pub fn new(iter: I) -> Self {
        // Alignment 1 starts with a pending chip (offset by 1)
        let mut align1 = AlignmentState::new();
        align1.pending = None; // will get its first chip naturally

        ManchesterDecoder {
            iter,
            align: [AlignmentState::new(), AlignmentState::new()],
            output_buf: None,
            total_chips: 0,
        }
    }

    fn active_alignment(&self) -> usize {
        if self.align[1].score() > self.align[0].score() { 1 } else { 0 }
    }
}

impl<I: Iterator<Item = f32>> Iterator for ManchesterDecoder<I> {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        // Return buffered output first
        if let Some(bit) = self.output_buf.take() {
            return Some(bit);
        }

        loop {
            let sample = self.iter.next()?;
            let chip = threshold(sample);
            self.total_chips += 1;

            // Alignment 0: pairs chips [0,1], [2,3], [4,5], ...
            let result0 = self.align[0].push_chip(chip);

            // Alignment 1: pairs chips [1,2], [3,4], [5,6], ...
            // Skip the first chip so it's offset by one
            let result1 = if self.total_chips >= 2 {
                self.align[1].push_chip(chip)
            } else {
                None
            };

            // During warmup, default to alignment 0
            if self.total_chips < SCORE_WINDOW * 2 + 2 {
                if let Some(bit) = result0 {
                    return Some(bit);
                }
                continue;
            }

            // After warmup, emit from the alignment with better score
            let active = self.active_alignment();

            match active {
                0 => {
                    if let Some(bit) = result0 {
                        return Some(bit);
                    }
                }
                1 => {
                    if let Some(bit) = result1 {
                        return Some(bit);
                    }
                }
                _ => unreachable!(),
            }
        }
    }
}

pub trait ManchesterDecodable {
    fn manchester_decode(self) -> ManchesterDecoder<Self>
    where
        Self: Sized + Iterator<Item = f32>;
}

impl<I: Iterator<Item = f32>> ManchesterDecodable for I {
    fn manchester_decode(self) -> ManchesterDecoder<Self> {
        ManchesterDecoder::new(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Encode data bits → differential → Manchester → f32 chips
    fn encode_rds_bits(data_bits: &[u8]) -> Vec<f32> {
        let mut diff_bits = Vec::with_capacity(data_bits.len());
        let mut prev = 0u8;
        for &b in data_bits {
            let d = b ^ prev;
            diff_bits.push(d);
            prev = d;
        }
        let mut chips = Vec::with_capacity(diff_bits.len() * 2);
        for &d in &diff_bits {
            if d == 0 {
                chips.push(1.0_f32);
                chips.push(-1.0_f32);
            } else {
                chips.push(-1.0_f32);
                chips.push(1.0_f32);
            }
        }
        chips
    }

    #[test]
    fn test_roundtrip() {
        let data_bits: Vec<u8> = (0..300).map(|i| ((i * 7 + 3) % 11 > 5) as u8).collect();
        let chips = encode_rds_bits(&data_bits);

        let decoded: Vec<u8> = chips.into_iter().manchester_decode().collect();

        let mut best_matches = 0;
        for offset in 0..5 {
            let check_len = decoded.len().min(data_bits.len() - offset);
            let matches: usize = (1..check_len)
                .filter(|&i| decoded[i] == data_bits[i + offset])
                .count();
            if matches > best_matches { best_matches = matches; }
        }

        let check_len = decoded.len().min(data_bits.len()) - 1;
        let accuracy = best_matches as f32 / check_len as f32;
        assert!(
            accuracy > 0.95,
            "Expected >95% roundtrip accuracy, got {:.1}%",
            accuracy * 100.0
        );
    }

    #[test]
    fn test_inverted_carrier_phase() {
        let data_bits: Vec<u8> = (0..300).map(|i| ((i * 13 + 5) % 7 > 3) as u8).collect();
        let chips: Vec<f32> = encode_rds_bits(&data_bits)
            .into_iter()
            .map(|c| -c)
            .collect();

        let decoded: Vec<u8> = chips.into_iter().manchester_decode().collect();

        let mut best_matches = 0;
        for offset in 0..5 {
            let check_len = decoded.len().min(data_bits.len() - offset);
            let matches: usize = (1..check_len)
                .filter(|&i| decoded[i] == data_bits[i + offset])
                .count();
            if matches > best_matches { best_matches = matches; }
        }

        let check_len = decoded.len().min(data_bits.len()) - 1;
        let accuracy = best_matches as f32 / check_len as f32;
        assert!(
            accuracy > 0.95,
            "Inverted phase: expected >95% accuracy, got {:.1}%",
            accuracy * 100.0
        );
    }

    #[test]
    fn test_offset_alignment() {
        let data_bits: Vec<u8> = (0..300).map(|i| (i % 3 > 0) as u8).collect();
        let mut chips = encode_rds_bits(&data_bits);
        // Insert one extra chip to force odd alignment
        chips.insert(0, 0.5);

        let decoded: Vec<u8> = chips.into_iter().manchester_decode().collect();
        assert!(decoded.len() > 100, "Should decode despite offset start");
    }

    #[test]
    fn test_selects_better_alignment() {
        // With proper alignment, all pairs are valid.
        // The dual tracker should converge to the correct one.
        let data_bits: Vec<u8> = (0..500).map(|i| ((i * 3 + 1) % 5 > 2) as u8).collect();
        let chips = encode_rds_bits(&data_bits);

        let decoded: Vec<u8> = chips.into_iter().manchester_decode().collect();

        // Should decode most bits correctly
        let mut best_matches = 0;
        for offset in 0..5 {
            let check_len = decoded.len().min(data_bits.len() - offset);
            let matches: usize = (1..check_len)
                .filter(|&i| decoded[i] == data_bits[i + offset])
                .count();
            if matches > best_matches { best_matches = matches; }
        }

        let check_len = decoded.len().min(data_bits.len()) - 1;
        let accuracy = best_matches as f32 / check_len as f32;
        assert!(
            accuracy > 0.90,
            "Expected >90% accuracy with dual alignment, got {:.1}%",
            accuracy * 100.0
        );
    }
}
