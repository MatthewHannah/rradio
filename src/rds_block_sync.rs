/// Developed using Claude Opus 4.6
/// RDS block synchronization, CRC checking, and group assembly.
///
/// Consumes a bitstream (Iterator<Item = u8>) and yields decoded RDS groups.
///
/// The RDS bitstream is structured as:
///   - Blocks of 26 bits: 16 data + 10 check (CRC + offset word)
///   - Groups of 4 blocks: A, B, C (or C'), D = 104 bits
///
/// Block sync is achieved by computing the CRC-10 syndrome on a sliding
/// 26-bit window and matching against known offset words.

/// CRC-10 generator polynomial: x^10 + x^8 + x^7 + x^5 + x^4 + x^3 + 1
/// Binary (degree 10 to 0): 10110111001 = 0x5B9
const CRC_POLY: u32 = 0x5B9; // 11-bit polynomial including leading x^10 term

/// Offset words for each block type (10 bits each)
const OFFSET_A:  u16 = 0x0FC;
const OFFSET_B:  u16 = 0x198;
const OFFSET_C:  u16 = 0x168;
const OFFSET_CP: u16 = 0x350; // C' (used in type B groups)
const OFFSET_D:  u16 = 0x1B4;

const OFFSETS: [(u16, &str); 5] = [
    (OFFSET_A,  "A"),
    (OFFSET_B,  "B"),
    (OFFSET_C,  "C"),
    (OFFSET_CP, "C'"),
    (OFFSET_D,  "D"),
];

/// Block index for each offset in the group sequence
const BLOCK_FOR_OFFSET: [usize; 5] = [0, 1, 2, 2, 3]; // A=0, B=1, C=2, C'=2, D=3

/// Expected next block offset indices after each block
/// After A(0) → B(1), after B(1) → C(2) or C'(3), after C/C'(2,3) → D(4), after D(4) → A(0)
fn expected_next_offsets(current_offset_idx: usize) -> &'static [usize] {
    match current_offset_idx {
        0 => &[1],       // A → B
        1 => &[2, 3],    // B → C or C'
        2 | 3 => &[4],   // C or C' → D
        4 => &[0],       // D → A
        _ => &[],
    }
}

const SYNC_THRESHOLD: usize = 3;  // consecutive good blocks to lock
const LOSS_THRESHOLD: usize = 5;  // consecutive bad blocks to lose sync

/// Compute CRC-10 syndrome for a 26-bit block using polynomial long division.
/// For an error-free block, syndrome == offset word for that block type.
fn syndrome(block: u32) -> u16 {
    // Treat the 26-bit block as a polynomial and divide by the generator.
    // The remainder after division is the syndrome.
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
fn check_block(block: u32, expected_offsets: &[usize], correct_errors: bool) -> Option<(usize, u16)> {
    let syn = syndrome(block);

    // First try exact match (no errors)
    for &idx in expected_offsets {
        if syn == OFFSETS[idx].0 {
            return Some((idx, (block >> 10) as u16));
        }
    }

    if !correct_errors {
        return None;
    }

    // Try error correction: XOR syndrome with each expected offset to get
    // the error syndrome, then look up in the correction table.
    // Only correct if the error pattern affects <= 2 bits (conservative).
    for &idx in expected_offsets {
        let error_syndrome = (syn ^ OFFSETS[idx].0) as usize;
        let error_pattern = ERROR_CORRECTION_TABLE[error_syndrome];
        if error_pattern != 0 && error_pattern.count_ones() <= 2 {
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
    pub blocks: [u16; 4],  // 16 data bits from blocks A, B, C/C', D
    pub group_type: u8,    // 0–15, from block B bits 12–15
    pub version: bool,     // false = type A, true = type B (block B bit 11)
    pub pi_code: u16,      // block A always carries the PI code
}

enum SyncState {
    Searching,
    Tentative { expected_offsets: &'static [usize], bits_remaining: usize },
    Locked,
}

pub struct RdsBlockSync<I: Iterator<Item = u8>> {
    iter: I,
    shift_reg: u32,        // 26-bit sliding window
    bit_count: usize,      // total bits shifted in
    state: SyncState,
    blocks: [u16; 4],      // accumulated block data
    block_idx: usize,      // current position in group (0=A, 1=B, 2=C, 3=D)
    current_offset_idx: usize, // which offset we last matched
    bits_in_block: usize,  // bits collected toward current block
    consecutive_good: usize,
    consecutive_bad: usize,
    synced_groups: usize,
}

impl<I: Iterator<Item = u8>> RdsBlockSync<I> {
    pub fn new(iter: I) -> Self {
        RdsBlockSync {
            iter,
            shift_reg: 0,
            bit_count: 0,
            state: SyncState::Searching,
            blocks: [0; 4],
            block_idx: 0,
            current_offset_idx: 0,
            bits_in_block: 0,
            consecutive_good: 0,
            consecutive_bad: 0,
            synced_groups: 0,
        }
    }

    fn push_bit(&mut self, bit: u8) -> Option<RdsGroup> {
        self.shift_reg = ((self.shift_reg << 1) | (bit as u32)) & 0x03FF_FFFF; // 26 bits
        self.bit_count += 1;
        self.bits_in_block += 1;

        match &self.state {
            SyncState::Searching => {
                if self.bit_count < 26 {
                    return None;
                }

                // Check syndrome against all offset words
                let syn = syndrome(self.shift_reg);
                for (idx, &(offset, _name)) in OFFSETS.iter().enumerate() {
                    if syn == offset {
                        let data = (self.shift_reg >> 10) as u16;
                        let block_pos = BLOCK_FOR_OFFSET[idx];
                        self.blocks = [0; 4];
                        self.blocks[block_pos] = data;
                        self.block_idx = block_pos;
                        self.current_offset_idx = idx;
                        self.bits_in_block = 0;
                        self.consecutive_good = 1;
                        self.consecutive_bad = 0;

                        self.state = SyncState::Tentative {
                            expected_offsets: expected_next_offsets(idx),
                            bits_remaining: 26,
                        };
                        return None;
                    }
                }
                None
            }

            SyncState::Tentative { expected_offsets, bits_remaining } => {
                let remaining = *bits_remaining - 1;
                if remaining > 0 {
                    self.state = SyncState::Tentative {
                        expected_offsets,
                        bits_remaining: remaining,
                    };
                    return None;
                }

                // 26 bits collected — check syndrome (no error correction during tentative)
                let expected = *expected_offsets;

                if let Some((idx, data)) = check_block(self.shift_reg, expected, false) {
                    let block_pos = BLOCK_FOR_OFFSET[idx];
                    self.blocks[block_pos] = data;
                    self.block_idx = block_pos;
                    self.current_offset_idx = idx;
                    self.bits_in_block = 0;
                    self.consecutive_good += 1;

                    let group = if block_pos == 3 {
                        // Completed group D — emit group
                        Some(self.emit_group())
                    } else {
                        None
                    };

                    if self.consecutive_good >= SYNC_THRESHOLD {
                        eprintln!("RDS sync: LOCKED after {} consecutive blocks", self.consecutive_good);
                        self.state = SyncState::Locked;
                    } else {
                        self.state = SyncState::Tentative {
                            expected_offsets: expected_next_offsets(idx),
                            bits_remaining: 26,
                        };
                    }
                    return group;
                } else {
                    self.state = SyncState::Searching;
                    self.consecutive_good = 0;
                    None
                }
            }

            SyncState::Locked => {
                if self.bits_in_block < 26 {
                    return None;
                }

                self.bits_in_block = 0;
                let expected = expected_next_offsets(self.current_offset_idx);

                if let Some((idx, data)) = check_block(self.shift_reg, expected, true) {
                    let block_pos = BLOCK_FOR_OFFSET[idx];
                    self.blocks[block_pos] = data;
                    self.block_idx = block_pos;
                    self.current_offset_idx = idx;
                    self.consecutive_good += 1;
                    self.consecutive_bad = 0;

                    if block_pos == 3 {
                        return Some(self.emit_group());
                    }
                } else {
                    self.consecutive_bad += 1;

                    // Advance expected position even on error
                    if let Some(&next) = expected.first() {
                        self.current_offset_idx = next;
                        self.block_idx = BLOCK_FOR_OFFSET[next];
                    }

                    if self.consecutive_bad >= LOSS_THRESHOLD {
                        eprintln!("RDS sync: LOST after {} consecutive bad blocks (total groups: {})",
                            self.consecutive_bad, self.synced_groups);
                        self.state = SyncState::Searching;
                        self.consecutive_good = 0;
                        self.consecutive_bad = 0;
                    }
                }
                None
            }
        }
    }

    fn emit_group(&mut self) -> RdsGroup {
        let b = self.blocks[1];
        let group = RdsGroup {
            blocks: self.blocks,
            group_type: ((b >> 12) & 0x0F) as u8,
            version: ((b >> 11) & 1) == 1,
            pi_code: self.blocks[0],
        };
        self.synced_groups += 1;
        self.blocks = [0; 4];
        group
    }
}

impl<I: Iterator<Item = u8>> Iterator for RdsBlockSync<I> {
    type Item = RdsGroup;

    fn next(&mut self) -> Option<RdsGroup> {
        loop {
            let bit = self.iter.next()?;
            if let Some(group) = self.push_bit(bit) {
                return Some(group);
            }
        }
    }
}

pub trait RdsBlockSyncable {
    fn rds_block_sync(self) -> RdsBlockSync<Self>
    where
        Self: Sized + Iterator<Item = u8>;
}

impl<I: Iterator<Item = u8>> RdsBlockSyncable for I {
    fn rds_block_sync(self) -> RdsBlockSync<Self> {
        RdsBlockSync::new(self)
    }
}

/// Accumulates RDS data across groups and prints complete strings.
pub struct RdsDecoder {
    ps: [u8; 8],          // Program Service name (8 chars)
    ps_filled: u8,        // bitmask of which PS segments (0–3) we've received
    ps_last_addr: usize,  // last PS address seen (detect wrap)

    rt: [u8; 64],         // RadioText (up to 64 chars)
    rt_len: usize,        // length up to \r or end
    rt_filled: u16,       // bitmask of which RT segments (0–15) we've received
    rt_last_addr: usize,  // last RT address seen (detect wrap)

    pi_code: u16,         // most recent valid PI code
}

impl RdsDecoder {
    pub fn new() -> Self {
        RdsDecoder {
            ps: [b' '; 8],
            ps_filled: 0,
            ps_last_addr: 0xFF,
            rt: [b' '; 64],
            rt_len: 64,
            rt_filled: 0,
            rt_last_addr: 0xFF,
            pi_code: 0,
        }
    }

    pub fn process(&mut self, group: &RdsGroup) {
        // Track PI code (ignore 0x0000 — likely corrupt)
        if group.pi_code != 0 {
            self.pi_code = group.pi_code;
        }

        match group.group_type {
            0 => self.decode_group_0(group),
            2 => self.decode_group_2(group),
            _ => {}
        }
    }

    fn decode_group_0(&mut self, group: &RdsGroup) {
        let addr = (group.blocks[1] & 0x03) as usize;

        // New cycle: address wrapped back to 0 — print what we have
        if addr == 0 && self.ps_filled != 0 && self.ps_last_addr != 0 {
            self.print_ps();
            self.ps = [b' '; 8];
            self.ps_filled = 0;
        }

        let c0 = ((group.blocks[3] >> 8) & 0xFF) as u8;
        let c1 = (group.blocks[3] & 0xFF) as u8;
        self.ps[addr * 2] = c0;
        self.ps[addr * 2 + 1] = c1;
        self.ps_filled |= 1 << addr;
        self.ps_last_addr = addr;

        // Print current accumulated PS after each segment
        self.print_ps();
    }

    fn decode_group_2(&mut self, group: &RdsGroup) {
        let addr = (group.blocks[1] & 0x0F) as usize;

        // New cycle: address wrapped back — print what we have
        if addr == 0 && self.rt_filled != 0 && self.rt_last_addr != 0 {
            self.print_rt();
            self.rt = [b' '; 64];
            self.rt_len = 64;
            self.rt_filled = 0;
        }

        if !group.version {
            // 2A: 4 chars from blocks C and D
            let chars = [
                ((group.blocks[2] >> 8) & 0xFF) as u8,
                (group.blocks[2] & 0xFF) as u8,
                ((group.blocks[3] >> 8) & 0xFF) as u8,
                (group.blocks[3] & 0xFF) as u8,
            ];
            let base = addr * 4;
            for (i, &c) in chars.iter().enumerate() {
                if base + i < 64 {
                    if c == 0x0D {
                        // \r marks end of RadioText
                        self.rt_len = base + i;
                    } else {
                        self.rt[base + i] = c;
                    }
                }
            }
        } else {
            // 2B: 2 chars from block D
            let chars = [
                ((group.blocks[3] >> 8) & 0xFF) as u8,
                (group.blocks[3] & 0xFF) as u8,
            ];
            let base = addr * 2;
            for (i, &c) in chars.iter().enumerate() {
                if base + i < 64 {
                    if c == 0x0D {
                        self.rt_len = base + i;
                    } else {
                        self.rt[base + i] = c;
                    }
                }
            }
        }

        self.rt_filled |= 1 << addr;
        self.rt_last_addr = addr;

        // Print current accumulated RT after each segment
        self.print_rt();
    }

    fn print_ps(&self) {
        let ps_str: String = self.ps.iter()
            .map(|&c| if c >= 0x20 && c < 0x7F { c as char } else { '?' })
            .collect::<String>()
            .trim_end()
            .to_string();
        if !ps_str.is_empty() {
            eprintln!("RDS PS:   \"{}\"  (PI=0x{:04X})", ps_str, self.pi_code);
        }
    }

    fn print_rt(&self) {
        let len = self.rt_len.min(64);
        let rt_str: String = self.rt[..len].iter()
            .map(|&c| if c >= 0x20 && c < 0x7F { c as char } else { '?' })
            .collect::<String>()
            .trim_end()
            .to_string();
        if !rt_str.is_empty() {
            eprintln!("RDS RT:   \"{}\"  (PI=0x{:04X})", rt_str, self.pi_code);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Compute the 10-bit check word for 16 data bits + offset word.
    fn encode_block(data: u16, offset: u16) -> u32 {
        // Place data bits in positions 25..10, zeros in 9..0
        let message = (data as u32) << 10;
        // CRC = remainder of polynomial division
        let crc = syndrome(message);
        // Check word = CRC XOR offset
        let check = crc ^ offset;
        ((data as u32) << 10) | (check as u32)
    }

    /// Build a bitstream from a sequence of 26-bit blocks
    fn blocks_to_bits(blocks: &[u32]) -> Vec<u8> {
        let mut bits = Vec::new();
        for &block in blocks {
            for i in (0..26).rev() {
                bits.push(((block >> i) & 1) as u8);
            }
        }
        bits
    }

    #[test]
    fn test_syndrome_matches_offset() {
        // An error-free block should have syndrome == its offset word
        let data: u16 = 0xABCD;
        for &(offset, name) in &OFFSETS {
            let block = encode_block(data, offset);
            let syn = syndrome(block);
            assert_eq!(syn, offset,
                "Block {} syndrome 0x{:03X} != offset 0x{:03X}",
                name, syn, offset);
        }
    }

    #[test]
    fn test_group_decode() {
        // Construct a valid group: A, B, C, D
        let pi: u16 = 0x1234;
        let block_b_data: u16 = 0x2400; // group type 2, version A, plus other bits
        let block_c_data: u16 = 0x4142; // "AB"
        let block_d_data: u16 = 0x4344; // "CD"

        let encoded_a = encode_block(pi, OFFSET_A);
        let encoded_b = encode_block(block_b_data, OFFSET_B);
        let encoded_c = encode_block(block_c_data, OFFSET_C);
        let encoded_d = encode_block(block_d_data, OFFSET_D);

        // Prepend some noise then the group, then repeat to help sync
        let mut blocks = vec![encoded_a, encoded_b, encoded_c, encoded_d];
        // Repeat for sync to lock
        for _ in 0..4 {
            blocks.push(encode_block(pi, OFFSET_A));
            blocks.push(encode_block(block_b_data, OFFSET_B));
            blocks.push(encode_block(block_c_data, OFFSET_C));
            blocks.push(encode_block(block_d_data, OFFSET_D));
        }

        let bits = blocks_to_bits(&blocks);
        let groups: Vec<RdsGroup> = bits.into_iter().rds_block_sync().collect();

        assert!(!groups.is_empty(), "Should decode at least one group");

        let g = &groups[0];
        assert_eq!(g.pi_code, pi);
        assert_eq!(g.group_type, 2);
        assert_eq!(g.version, false); // type A
        assert_eq!(g.blocks[2], block_c_data);
        assert_eq!(g.blocks[3], block_d_data);
    }

    #[test]
    fn test_sync_recovers_after_noise() {
        let pi: u16 = 0x5678;
        let b_data: u16 = 0x0800; // group type 0, version B
        let c_data: u16 = 0x0000;
        let d_data: u16 = 0x4849; // "HI"

        // 50 random noise bits
        let mut bits: Vec<u8> = (0..50).map(|i| ((i * 7 + 3) % 2) as u8).collect();

        // Then valid groups
        for _ in 0..6 {
            let blocks = vec![
                encode_block(pi, OFFSET_A),
                encode_block(b_data, OFFSET_B),
                encode_block(c_data, OFFSET_C),
                encode_block(d_data, OFFSET_D),
            ];
            bits.extend(blocks_to_bits(&blocks));
        }

        let groups: Vec<RdsGroup> = bits.into_iter().rds_block_sync().collect();
        assert!(groups.len() >= 2, "Should recover sync and decode groups after noise, got {}", groups.len());
    }
}
