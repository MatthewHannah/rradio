/// Combined biphase + block synchronization for RDS.
///
/// Takes raw chips (Complex32 at 2375 Hz) and produces SyncEvents.
/// During Searching, two parallel biphase decoders (phase 0 and phase 1)
/// race to find a CRC match. Once sync is acquired, the winning phase is
/// frozen — no energy-based polarity re-evaluation.
///
/// Also contains CRC primitives, offset words, and error correction for RDS blocks.

use num_complex::Complex32;

// ── CRC and block-level constants ──────────────────────────────────────────

/// CRC-10 generator polynomial: x^10 + x^8 + x^7 + x^5 + x^4 + x^3 + 1
/// Binary (degree 10 to 0): 10110111001 = 0x5B9
const CRC_POLY: u32 = 0x5B9; // 11-bit polynomial including leading x^10 term

/// Offset words for each block type (10 bits each)
const OFFSET_A:  u16 = 0x0FC;
const OFFSET_B:  u16 = 0x198;
const OFFSET_C:  u16 = 0x168;
const OFFSET_CP: u16 = 0x350; // C' (used in type B groups)
const OFFSET_D:  u16 = 0x1B4;

pub const OFFSETS: [(u16, &str); 5] = [
    (OFFSET_A,  "A"),
    (OFFSET_B,  "B"),
    (OFFSET_C,  "C"),
    (OFFSET_CP, "C'"),
    (OFFSET_D,  "D"),
];

/// Block index for each offset in the group sequence
pub const BLOCK_FOR_OFFSET: [usize; 5] = [0, 1, 2, 2, 3]; // A=0, B=1, C=2, C'=2, D=3

/// Expected next block offset indices after each block
pub fn expected_next_offsets(current_offset_idx: usize) -> &'static [usize] {
    match current_offset_idx {
        0 => &[1],       // A → B
        1 => &[2, 3],    // B → C or C'
        2 | 3 => &[4],   // C or C' → D
        4 => &[0],       // D → A
        _ => &[],
    }
}

pub const SYNC_THRESHOLD: usize = 3;  // consecutive good blocks to lock

/// Compute CRC-10 syndrome for a 26-bit block using polynomial long division.
/// For an error-free block, syndrome == offset word for that block type.
pub fn syndrome(block: u32) -> u16 {
    let mut reg = block;
    for i in (10..26).rev() {
        if reg & (1 << i) != 0 {
            reg ^= CRC_POLY << (i - 10);
        }
    }
    (reg & 0x3FF) as u16
}

/// Try to match a block's syndrome against expected offsets, with optional error correction.
/// Returns Some((offset_index, corrected_data)) on success.
pub fn check_block(block: u32, expected_offsets: &[usize], max_correction_bits: u32) -> Option<(usize, u16)> {
    let syn = syndrome(block);

    // First try exact match (no errors)
    for &idx in expected_offsets {
        if syn == OFFSETS[idx].0 {
            return Some((idx, (block >> 10) as u16));
        }
    }

    if max_correction_bits == 0 {
        return None;
    }

    // Try error correction
    for &idx in expected_offsets {
        let error_syndrome = (syn ^ OFFSETS[idx].0) as usize;
        let error_pattern = ERROR_CORRECTION_TABLE[error_syndrome];
        if error_pattern != 0 && error_pattern.count_ones() <= max_correction_bits {
            let corrected = block ^ error_pattern;
            return Some((idx, (corrected >> 10) as u16));
        }
    }

    None
}

/// Error correction lookup table: syndrome → 26-bit error pattern.
/// Covers single-bit errors and burst errors up to 5 bits (367 patterns).
/// Zero entry means uncorrectable.
#[rustfmt::skip]
static ERROR_CORRECTION_TABLE: [u32; 1024] = include!("rds_error_table.inc");

#[derive(Debug, Clone)]
pub struct RdsGroup {
    pub blocks: [u16; 4],
    pub group_type: u8,
    pub version: bool,
    pub pi_code: u16,
    pub rolling_bler: f64,
}

/// Events yielded by the block sync iterator.
pub enum SyncEvent {
    /// A complete group was decoded.
    Group(RdsGroup),
    /// Sync was just lost.
    LostSync,
    /// Sync was just acquired.
    Locked,
    /// Still searching (emitted periodically).
    Searching,
}

// ── Biphase + block synchronization ────────────────────────────────────────

/// A single biphase + differential decoding path.
struct BiphasePath {
    prev_chip: Complex32,
    prev_encoded_bit: bool,
    has_prev_chip: bool,
    /// Which chips produce output: phase 0 outputs on even chips, phase 1 on odd.
    phase: usize,
    chip_count: usize,
    shift_reg: u32,
    bits_pushed: usize,
}

impl BiphasePath {
    fn new(phase: usize) -> Self {
        BiphasePath {
            prev_chip: Complex32::new(0.0, 0.0),
            prev_encoded_bit: false,
            has_prev_chip: false,
            phase,
            chip_count: 0,
            shift_reg: 0,
            bits_pushed: 0,
        }
    }

    /// Push a chip. Returns Some(data_bit) when a complete biphase symbol is decoded.
    fn push_chip(&mut self, chip: Complex32) -> Option<bool> {
        let result = if self.has_prev_chip {
            // Biphase: difference of consecutive chips
            let biphase = (chip - self.prev_chip) * 0.5;
            let encoded_bit = biphase.re >= 0.0;
            let has_value = (self.chip_count % 2) == self.phase;

            if has_value {
                // Differential decode
                let data_bit = encoded_bit != self.prev_encoded_bit;
                self.prev_encoded_bit = encoded_bit;
                Some(data_bit)
            } else {
                // Update encoded_bit tracking on the non-output chip too
                // (the "other half" of the Manchester pair)
                None
            }
        } else {
            None
        };

        self.prev_chip = chip;
        self.has_prev_chip = true;
        self.chip_count += 1;

        // If we got a data bit, push into shift register
        if let Some(bit) = result {
            self.shift_reg = ((self.shift_reg << 1) | (bit as u32)) & 0x03FF_FFFF;
            self.bits_pushed += 1;
        }

        result
    }

    fn reset_shift_reg(&mut self) {
        self.shift_reg = 0;
        self.bits_pushed = 0;
    }
}

enum ChipSyncState {
    /// Both phases active, looking for CRC match
    Searching,
    /// One phase matched, confirming with consecutive blocks
    Tentative {
        phase: usize,
        expected_offsets: &'static [usize],
        bits_remaining: usize,
        consecutive_good: usize,
    },
    /// Locked on a single phase
    Locked {
        phase: usize,
    },
}

pub struct ChipSync {
    paths: [BiphasePath; 2],
    state: ChipSyncState,

    // Block sync state (shared, used in Tentative/Locked)
    blocks: [u16; 4],
    block_idx: usize,
    current_offset_idx: usize,
    bits_in_block: usize,
    consecutive_bad: usize,
    synced_groups: usize,
    bit_count: usize,
    debug: bool,
    loss_threshold: usize,
    crc_max_bits: u32,

    // BLER tracking
    total_blocks_checked: u64,
    total_blocks_passed: u64,
    bler_history: Vec<bool>,
    bler_idx: usize,
    bler_passed: usize,

    // Search rate limiting
    searching_counter: usize,

    // Diagnostics
    pub polarity_at_lock: usize,
    pub total_chips: u64,
}

const BLER_WINDOW: usize = 200;
const SEARCHING_REPORT_INTERVAL: usize = 1187;

impl ChipSync {
    pub fn new(crc_correction_max_bits: u32, loss_threshold: usize, debug: bool) -> Self {
        ChipSync {
            paths: [BiphasePath::new(0), BiphasePath::new(1)],
            state: ChipSyncState::Searching,
            blocks: [0; 4],
            block_idx: 0,
            current_offset_idx: 0,
            bits_in_block: 0,
            consecutive_bad: 0,
            synced_groups: 0,
            bit_count: 0,
            debug,
            loss_threshold,
            crc_max_bits: crc_correction_max_bits,
            total_blocks_checked: 0,
            total_blocks_passed: 0,
            bler_history: vec![false; BLER_WINDOW],
            bler_idx: 0,
            bler_passed: 0,
            searching_counter: 0,
            polarity_at_lock: 0,
            total_chips: 0,
        }
    }

    fn record_block(&mut self, passed: bool) {
        self.total_blocks_checked += 1;
        if passed { self.total_blocks_passed += 1; }
        let old = self.bler_history[self.bler_idx];
        if old { self.bler_passed -= 1; }
        self.bler_history[self.bler_idx] = passed;
        if passed { self.bler_passed += 1; }
        self.bler_idx = (self.bler_idx + 1) % BLER_WINDOW;
    }

    fn rolling_bler(&self) -> f64 {
        let window = self.total_blocks_checked.min(BLER_WINDOW as u64) as usize;
        if window == 0 { return 0.0; }
        1.0 - self.bler_passed as f64 / window as f64
    }

    fn emit_group(&mut self) -> RdsGroup {
        let b = self.blocks[1];
        let group = RdsGroup {
            blocks: self.blocks,
            group_type: ((b >> 12) & 0x0F) as u8,
            version: ((b >> 11) & 1) == 1,
            pi_code: self.blocks[0],
            rolling_bler: self.rolling_bler(),
        };
        self.synced_groups += 1;
        self.blocks = [0; 4];
        group
    }

    /// Push one chip (Complex32 at 2375 Hz). Returns SyncEvent when available.
    pub fn push_chip(&mut self, chip: Complex32) -> Option<SyncEvent> {
        self.total_chips += 1;

        // Feed both paths always (they maintain their own prev_chip state)
        let bit_0 = self.paths[0].push_chip(chip);
        let bit_1 = self.paths[1].push_chip(chip);

        match &self.state {
            ChipSyncState::Searching => {
                // Check both paths for CRC matches
                for phase in 0..2 {
                    let path = &self.paths[phase];
                    if path.bits_pushed < 26 { continue; }

                    let syn = syndrome(path.shift_reg);
                    for (idx, &(offset, name)) in OFFSETS.iter().enumerate() {
                        if syn == offset {
                            if self.debug {
                                eprintln!("RDS search: FOUND {} at phase={} bit={} (syn=0x{:03X})",
                                    name, phase, path.bits_pushed, syn);
                            }
                            let data = (path.shift_reg >> 10) as u16;
                            let block_pos = BLOCK_FOR_OFFSET[idx];
                            self.blocks = [0; 4];
                            self.blocks[block_pos] = data;
                            self.block_idx = block_pos;
                            self.current_offset_idx = idx;
                            self.bits_in_block = 0;
                            self.bit_count = path.bits_pushed;

                            self.state = ChipSyncState::Tentative {
                                phase,
                                expected_offsets: expected_next_offsets(idx),
                                bits_remaining: 26,
                                consecutive_good: 1,
                            };
                            return None;
                        }
                    }
                }

                // Periodic searching notification
                self.searching_counter += 1;
                if self.searching_counter >= SEARCHING_REPORT_INTERVAL {
                    self.searching_counter = 0;
                    return Some(SyncEvent::Searching);
                }
                None
            }

            ChipSyncState::Tentative { phase, expected_offsets, bits_remaining, consecutive_good } => {
                let phase = *phase;
                let expected = *expected_offsets;
                let consecutive_good = *consecutive_good;

                // Only advance on bits from the tentative phase
                let got_bit = if phase == 0 { bit_0 } else { bit_1 };
                if got_bit.is_none() { return None; }

                self.bit_count += 1;
                self.bits_in_block += 1;

                let remaining = *bits_remaining - 1;
                if remaining > 0 {
                    self.state = ChipSyncState::Tentative {
                        phase, expected_offsets: expected, bits_remaining: remaining,
                        consecutive_good,
                    };
                    return None;
                }

                // 26 bits collected — check CRC
                let sr = self.paths[phase].shift_reg;
                if let Some((idx, data)) = check_block(sr, expected, 0) {
                    let new_good = consecutive_good + 1;
                    if self.debug {
                        eprintln!("RDS tentative: GOOD block {} ({}/{}) at phase={}",
                            OFFSETS[idx].1, new_good, SYNC_THRESHOLD, phase);
                    }
                    let block_pos = BLOCK_FOR_OFFSET[idx];
                    self.blocks[block_pos] = data;
                    self.block_idx = block_pos;
                    self.current_offset_idx = idx;
                    self.bits_in_block = 0;

                    let group = if block_pos == 3 {
                        Some(SyncEvent::Group(self.emit_group()))
                    } else {
                        None
                    };

                    if new_good >= SYNC_THRESHOLD {
                        if self.debug {
                            eprintln!("RDS sync: LOCKED at phase={} after {} consecutive blocks",
                                phase, new_good);
                        }
                        self.polarity_at_lock = phase;
                        self.consecutive_bad = 0;
                        self.state = ChipSyncState::Locked { phase };
                        return group.or(Some(SyncEvent::Locked));
                    } else {
                        self.state = ChipSyncState::Tentative {
                            phase,
                            expected_offsets: expected_next_offsets(idx),
                            bits_remaining: 26,
                            consecutive_good: new_good,
                        };
                    }
                    return group;
                } else {
                    if self.debug {
                        let syn = syndrome(self.paths[phase].shift_reg);
                        eprintln!("RDS tentative: FAIL at phase={} (syn=0x{:03X}, expected {:?}, had {} good)",
                            phase, syn,
                            expected.iter().map(|&i| OFFSETS[i].1).collect::<Vec<_>>(),
                            consecutive_good);
                    }
                    // Reset both paths and go back to searching
                    self.paths[0].reset_shift_reg();
                    self.paths[1].reset_shift_reg();
                    self.state = ChipSyncState::Searching;
                    None
                }
            }

            ChipSyncState::Locked { phase } => {
                let phase = *phase;

                // Only process bits from the locked phase
                let got_bit = if phase == 0 { bit_0 } else { bit_1 };
                if got_bit.is_none() { return None; }

                self.bit_count += 1;
                self.bits_in_block += 1;

                if self.bits_in_block < 26 { return None; }

                self.bits_in_block = 0;
                let expected = expected_next_offsets(self.current_offset_idx);
                let sr = self.paths[phase].shift_reg;

                if let Some((idx, data)) = check_block(sr, expected, self.crc_max_bits) {
                    self.record_block(true);
                    let block_pos = BLOCK_FOR_OFFSET[idx];
                    self.blocks[block_pos] = data;
                    self.block_idx = block_pos;
                    self.current_offset_idx = idx;

                    // Exact CRC for health tracking
                    if check_block(sr, expected, 0).is_some() {
                        self.consecutive_bad = 0;
                    }

                    if block_pos == 3 {
                        return Some(SyncEvent::Group(self.emit_group()));
                    }
                } else {
                    self.record_block(false);
                    self.consecutive_bad += 1;

                    if let Some(&next) = expected.first() {
                        self.current_offset_idx = next;
                        self.block_idx = BLOCK_FOR_OFFSET[next];
                    }

                    if self.consecutive_bad >= self.loss_threshold {
                        if self.debug {
                            eprintln!("RDS sync: LOST after {} consecutive bad at chip {} (groups: {})",
                                self.consecutive_bad, self.total_chips, self.synced_groups);
                        }
                        self.consecutive_bad = 0;
                        self.paths[0].reset_shift_reg();
                        self.paths[1].reset_shift_reg();
                        self.state = ChipSyncState::Searching;
                        return Some(SyncEvent::LostSync);
                    }
                }
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Compute the 10-bit check word for 16 data bits + offset word.
    fn encode_block(data: u16, offset: u16) -> u32 {
        let message = (data as u32) << 10;
        let crc = syndrome(message);
        let check = crc ^ offset;
        ((data as u32) << 10) | (check as u32)
    }

    #[test]
    fn test_syndrome_round_trip() {
        let data: u16 = 0xABCD;
        let block = encode_block(data, OFFSET_A);
        assert_eq!(syndrome(block), OFFSET_A,
            "syndrome of a valid block-A should equal OFFSET_A");
    }

    #[test]
    fn test_check_block_exact() {
        let data: u16 = 0xABCD;
        let block = encode_block(data, OFFSET_A);
        let result = check_block(block, &[0], 0);
        assert_eq!(result, Some((0, data)));
    }

    #[test]
    fn test_check_block_with_1bit_error() {
        let data: u16 = 0xABCD;
        let block = encode_block(data, OFFSET_A);
        // Flip one bit (bit 20)
        let corrupted = block ^ (1 << 20);
        // Should fail with 0-bit correction
        assert!(check_block(corrupted, &[0], 0).is_none());
        // Should succeed with 1-bit correction and recover original data
        let result = check_block(corrupted, &[0], 1);
        assert_eq!(result, Some((0, data)));
    }

    #[test]
    fn test_end_to_end_chip_sync() {
        let pi: u16 = 0x1234;
        let block_b: u16 = 0x2400; // group type 2, version A
        let block_c: u16 = 0x4142;
        let block_d: u16 = 0x4344;

        // Encode blocks
        let raw_blocks = [
            encode_block(pi, OFFSET_A),
            encode_block(block_b, OFFSET_B),
            encode_block(block_c, OFFSET_C),
            encode_block(block_d, OFFSET_D),
        ];

        // Convert to bit stream (26 bits per block, MSB first)
        let mut bits: Vec<bool> = Vec::new();
        for &blk in &raw_blocks {
            for i in (0..26).rev() {
                bits.push((blk >> i) & 1 == 1);
            }
        }

        // Repeat groups so sync can lock (need SYNC_THRESHOLD + extra)
        let one_group = bits.clone();
        for _ in 0..7 {
            bits.extend_from_slice(&one_group);
        }

        // Differential encode: encoded[n] = data[n] XOR encoded[n-1]
        let mut encoded_bits: Vec<bool> = Vec::with_capacity(bits.len());
        let mut prev = false;
        for &b in &bits {
            let enc = b ^ prev;
            encoded_bits.push(enc);
            prev = enc;
        }

        // Manchester encode: true → [-1, +1], false → [+1, -1]
        // (Convention matches biphase decoder: chip - prev_chip > 0 → true)
        let mut chips: Vec<Complex32> = Vec::with_capacity(encoded_bits.len() * 2);
        for &enc in &encoded_bits {
            if enc {
                chips.push(Complex32::new(-1.0, 0.0));
                chips.push(Complex32::new(1.0, 0.0));
            } else {
                chips.push(Complex32::new(1.0, 0.0));
                chips.push(Complex32::new(-1.0, 0.0));
            }
        }

        // Feed chips into ChipSync
        let mut sync = ChipSync::new(2, 12, false);
        let mut groups: Vec<RdsGroup> = Vec::new();
        for &chip in &chips {
            if let Some(SyncEvent::Group(g)) = sync.push_chip(chip) {
                groups.push(g);
            }
        }

        assert!(!groups.is_empty(), "Should decode at least one group");
        let g = groups.iter().find(|g| g.pi_code != 0)
            .expect("Should decode at least one group with non-zero PI");
        assert_eq!(g.pi_code, pi, "PI code mismatch");
        assert_eq!(g.group_type, 2, "Group type mismatch");
        assert!(!g.version, "Should be version A");
    }
}