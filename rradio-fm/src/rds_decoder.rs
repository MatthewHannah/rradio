/// RDS data decoding and display.
///
/// Accumulates PS (programme service name), RT (radio text), and PI code
/// from decoded RDS groups, and renders them to stderr.

use crate::chip_sync::RdsGroup;

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
