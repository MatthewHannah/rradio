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
fn check_block(block: u32, expected_offsets: &[usize], max_correction_bits: u32) -> Option<(usize, u16)> {
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
    debug: bool,
    loss_threshold: usize,
    crc_max_bits: u32,
    total_blocks_checked: u64,
    total_blocks_passed: u64,
    bler_history: Vec<bool>,  // ring buffer: true = passed
    bler_idx: usize,
    bler_passed: usize,
    searching_counter: usize,
}

const BLER_WINDOW: usize = 200;
const SEARCHING_REPORT_INTERVAL: usize = 1187; // ~1 second of bits

impl<I: Iterator<Item = u8>> RdsBlockSync<I> {
    pub fn new(iter: I, crc_correction_max_bits: u32, loss_threshold: usize, debug: bool) -> Self {
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
            debug,
            loss_threshold,
            crc_max_bits: crc_correction_max_bits,
            total_blocks_checked: 0,
            total_blocks_passed: 0,
            bler_history: vec![false; BLER_WINDOW],
            bler_idx: 0,
            bler_passed: 0,
            searching_counter: 0,
        }
    }

    fn record_block(&mut self, passed: bool) {
        self.total_blocks_checked += 1;
        if passed {
            self.total_blocks_passed += 1;
        }
        // Rolling window
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

    pub fn push_bit(&mut self, bit: u8) -> Option<SyncEvent> {
        self.shift_reg = ((self.shift_reg << 1) | (bit as u32)) & 0x03FF_FFFF; // 26 bits
        self.bit_count += 1;
        self.bits_in_block += 1;

        match &self.state {
            SyncState::Searching => {
                if self.bit_count < 26 {
                    return None;
                }

                // Periodic searching notification
                self.searching_counter += 1;
                let emit_searching = if self.searching_counter >= SEARCHING_REPORT_INTERVAL {
                    self.searching_counter = 0;
                    true
                } else {
                    false
                };

                // Check syndrome against all offset words
                let syn = syndrome(self.shift_reg);
                for (idx, &(offset, name)) in OFFSETS.iter().enumerate() {
                    if syn == offset {
                        if self.debug {
                            eprintln!("RDS search: FOUND {} at bit {} (syn=0x{:03X})",
                                name, self.bit_count, syn);
                        }
                        let data = (self.shift_reg >> 10) as u16;
                        let block_pos = BLOCK_FOR_OFFSET[idx];
                        self.blocks = [0; 4];
                        self.blocks[block_pos] = data;
                        self.block_idx = block_pos;
                        self.current_offset_idx = idx;
                        self.bits_in_block = 0;
                        self.consecutive_good = 1;
                        self.consecutive_bad = 0;
                        self.searching_counter = 0;

                        self.state = SyncState::Tentative {
                            expected_offsets: expected_next_offsets(idx),
                            bits_remaining: 26,
                        };
                        return None;
                    }
                }
                if emit_searching { Some(SyncEvent::Searching) } else { None }
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

                if let Some((idx, data)) = check_block(self.shift_reg, expected, 0) {
                    if self.debug {
                        eprintln!("RDS tentative: GOOD block {} ({}/{} consecutive) at bit {}",
                            OFFSETS[idx].1, self.consecutive_good + 1, SYNC_THRESHOLD, self.bit_count);
                    }
                    let block_pos = BLOCK_FOR_OFFSET[idx];
                    self.blocks[block_pos] = data;
                    self.block_idx = block_pos;
                    self.current_offset_idx = idx;
                    self.bits_in_block = 0;
                    self.consecutive_good += 1;

                    let group = if block_pos == 3 {
                        Some(SyncEvent::Group(self.emit_group()))
                    } else {
                        None
                    };

                    if self.consecutive_good >= SYNC_THRESHOLD {
                        if self.debug {
                            eprintln!("RDS sync: LOCKED after {} consecutive blocks", self.consecutive_good);
                        }
                        self.state = SyncState::Locked;
                        // If we also have a group, return that; otherwise return Locked
                        return group.or(Some(SyncEvent::Locked));
                    } else {
                        self.state = SyncState::Tentative {
                            expected_offsets: expected_next_offsets(idx),
                            bits_remaining: 26,
                        };
                    }
                    return group;
                } else {
                    if self.debug {
                        let syn = syndrome(self.shift_reg);
                        eprintln!("RDS tentative: FAIL at bit {} (syn=0x{:03X}, expected {:?}, had {} good)",
                            self.bit_count, syn,
                            expected.iter().map(|&i| OFFSETS[i].1).collect::<Vec<_>>(),
                            self.consecutive_good);
                    }
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

                // Check with error correction for data recovery
                if let Some((idx, data)) = check_block(self.shift_reg, expected, self.crc_max_bits) {
                    self.record_block(true);
                    let block_pos = BLOCK_FOR_OFFSET[idx];
                    self.blocks[block_pos] = data;
                    self.block_idx = block_pos;
                    self.current_offset_idx = idx;
                    self.consecutive_good += 1;

                    // Track sync health using EXACT CRC only (no error correction).
                    // Error-corrected blocks that aren't exact matches indicate
                    // degraded signal; only exact matches reset the bad counter.
                    if check_block(self.shift_reg, expected, 0).is_some() {
                        self.consecutive_bad = 0;
                    }

                    if block_pos == 3 {
                        return Some(SyncEvent::Group(self.emit_group()));
                    }
                } else {
                    self.record_block(false);
                    self.consecutive_bad += 1;

                    // Advance expected position even on error
                    if let Some(&next) = expected.first() {
                        self.current_offset_idx = next;
                        self.block_idx = BLOCK_FOR_OFFSET[next];
                    }

                    if self.consecutive_bad >= self.loss_threshold {
                        if self.debug {
                            eprintln!("RDS sync: LOST after {} consecutive bad blocks at bit {} (total groups: {})",
                                self.consecutive_bad, self.bit_count, self.synced_groups);
                        }
                        self.state = SyncState::Searching;
                        self.consecutive_good = 0;
                        self.consecutive_bad = 0;
                        return Some(SyncEvent::LostSync);
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
            rolling_bler: self.rolling_bler(),
        };
        self.synced_groups += 1;
        self.blocks = [0; 4];
        group
    }
}

impl<I: Iterator<Item = u8>> Iterator for RdsBlockSync<I> {
    type Item = SyncEvent;

    fn next(&mut self) -> Option<SyncEvent> {
        loop {
            let bit = self.iter.next()?;
            if let Some(event) = self.push_bit(bit) {
                return Some(event);
            }
        }
    }
}

pub trait RdsBlockSyncable {
    fn rds_block_sync(self, crc_correction_max_bits: u32, loss_threshold: usize, debug: bool) -> RdsBlockSync<Self>
    where
        Self: Sized + Iterator<Item = u8>;
}

impl<I: Iterator<Item = u8>> RdsBlockSyncable for I {
    fn rds_block_sync(self, crc_correction_max_bits: u32, loss_threshold: usize, debug: bool) -> RdsBlockSync<Self> {
        RdsBlockSync::new(self, crc_correction_max_bits, loss_threshold, debug)
    }
}


/// Current RDS display state returned by the decoder.
pub struct RdsDisplayState {
    pub ps: String,
    pub rt: String,
    pub pi_code: u16,
    pub groups_decoded: u64,
}

/// Accumulates RDS data across groups.
pub struct RdsDecoder {
    ps: [u8; 8],
    ps_filled: u8,
    ps_last_addr: usize,

    rt: [u8; 64],
    rt_len: usize,
    rt_filled: u16,
    rt_last_addr: usize,

    pi_code: u16,
    groups_decoded: u64,
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
            groups_decoded: 0,
        }
    }

    pub fn process(&mut self, group: &RdsGroup) -> RdsDisplayState {
        if group.pi_code != 0 {
            self.pi_code = group.pi_code;
        }
        self.groups_decoded += 1;

        match group.group_type {
            0 => self.decode_group_0(group),
            2 => self.decode_group_2(group),
            _ => {}
        }

        self.display_state()
    }

    pub fn display_state(&self) -> RdsDisplayState {
        let ps = self.ps.iter()
            .map(|&c| Self::sanitize(c))
            .collect::<String>()
            .trim_end()
            .to_string();
        let len = self.rt_len.min(64);
        let rt = self.rt[..len].iter()
            .map(|&c| Self::sanitize(c))
            .collect::<String>()
            .trim_end()
            .to_string();
        RdsDisplayState {
            ps,
            rt,
            pi_code: self.pi_code,
            groups_decoded: self.groups_decoded,
        }
    }

    fn decode_group_0(&mut self, group: &RdsGroup) {
        let addr = (group.blocks[1] & 0x03) as usize;
        let c0 = ((group.blocks[3] >> 8) & 0xFF) as u8;
        let c1 = (group.blocks[3] & 0xFF) as u8;
        self.ps[addr * 2] = c0;
        self.ps[addr * 2 + 1] = c1;
        self.ps_filled |= 1 << addr;
        self.ps_last_addr = addr;
    }

    fn decode_group_2(&mut self, group: &RdsGroup) {
        let addr = (group.blocks[1] & 0x0F) as usize;
        if !group.version {
            let chars = [
                ((group.blocks[2] >> 8) & 0xFF) as u8,
                (group.blocks[2] & 0xFF) as u8,
                ((group.blocks[3] >> 8) & 0xFF) as u8,
                (group.blocks[3] & 0xFF) as u8,
            ];
            let base = addr * 4;
            for (i, &c) in chars.iter().enumerate() {
                if base + i < 64 {
                    if c == 0x0D { self.rt_len = base + i; }
                    else { self.rt[base + i] = c; }
                }
            }
        } else {
            let chars = [
                ((group.blocks[3] >> 8) & 0xFF) as u8,
                (group.blocks[3] & 0xFF) as u8,
            ];
            let base = addr * 2;
            for (i, &c) in chars.iter().enumerate() {
                if base + i < 64 {
                    if c == 0x0D { self.rt_len = base + i; }
                    else { self.rt[base + i] = c; }
                }
            }
        }
        self.rt_filled |= 1 << addr;
        self.rt_last_addr = addr;
    }

    fn sanitize(c: u8) -> char {
        if c >= 0x20 && c < 0x7F { c as char } else { ' ' }
    }
}

/// Renders the RDS display to stderr using ANSI cursor control.
pub struct RdsDisplay {
    drawn: bool,
    synced: bool,
    last_ps: String,
    last_rt: String,
    last_pi: u16,
    last_groups: u64,
    last_synced: bool,
}

impl RdsDisplay {
    pub fn new() -> Self {
        RdsDisplay {
            drawn: false,
            synced: false,
            last_ps: String::new(),
            last_rt: String::new(),
            last_pi: 0,
            last_groups: 0,
            last_synced: false,
        }
    }

    pub fn set_synced(&mut self, synced: bool) {
        self.synced = synced;
    }

    pub fn render(&mut self, state: &RdsDisplayState) {
        // Skip redraw if nothing changed
        if self.drawn
            && state.ps == self.last_ps
            && state.rt == self.last_rt
            && state.pi_code == self.last_pi
            && state.groups_decoded == self.last_groups
            && self.synced == self.last_synced
        {
            return;
        }

        self.last_ps = state.ps.clone();
        self.last_rt = state.rt.clone();
        self.last_pi = state.pi_code;
        self.last_groups = state.groups_decoded;
        self.last_synced = self.synced;
        let pi = if state.pi_code != 0 {
            format!("{:04X}", state.pi_code)
        } else {
            "----".to_string()
        };

        if self.drawn {
            eprint!("\x1b[4A");
        }
        self.drawn = true;

        let width = 56;
        let sync_icon = if self.synced { "✓" } else { "✗" };
        let top_line = format!("  {}    PI: {}  Sync: {}  Groups: {}",
            state.ps, pi, sync_icon, state.groups_decoded);
        let rt_display = if state.rt.len() > width - 4 {
            &state.rt[..width - 4]
        } else {
            &state.rt
        };

        eprintln!("\u{250c}{}\u{2510}", "\u{2500}".repeat(width));
        eprintln!("\u{2502}{:<width$}\u{2502}", top_line, width = width);
        eprintln!("\u{2502}  {:<w$}\u{2502}", rt_display, w = width - 2);
        eprintln!("\u{2514}{}\u{2518}", "\u{2500}".repeat(width));
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
        let groups: Vec<RdsGroup> = bits.into_iter().rds_block_sync(2, 12, false).filter_map(|e| match e { SyncEvent::Group(g) => Some(g), _ => None }).collect();

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

        let groups: Vec<RdsGroup> = bits.into_iter().rds_block_sync(2, 12, false).filter_map(|e| match e { SyncEvent::Group(g) => Some(g), _ => None }).collect();
        assert!(groups.len() >= 2, "Should recover sync and decode groups after noise, got {}", groups.len());
    }
}
