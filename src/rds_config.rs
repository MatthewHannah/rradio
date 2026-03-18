/// RDS pipeline configuration.
///
/// All parameters for tuning the RDS decoding pipeline. Can be loaded
/// from a JSON file via `--rds-config <path>`, or uses defaults.
/// Each component owns its own config struct.

use serde::Deserialize;

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct RdsConfig {
    /// LPF cutoff in wfm_audio after 57 kHz mix (Hz)
    pub baseband_lpf_hz: f32,
    /// Anti-alias filter cutoff before decimation (Hz)
    pub anti_alias_cutoff_hz: f32,
    /// Number of cascaded anti-alias biquad stages
    pub anti_alias_order: usize,
    /// Chip periods to span for the matched filter
    pub matched_filter_spans: usize,
    /// Window function for the matched filter
    pub matched_filter_window: String,
    /// AGC config
    pub agc: AgcConfig,
    /// Costas loop config
    pub costas: CostasConfig,
    /// Gardner clock recovery config
    pub gardner: GardnerConfig,
    /// Block sync config
    pub sync: SyncConfig,
}

impl Default for RdsConfig {
    fn default() -> Self {
        RdsConfig {
            baseband_lpf_hz: 4000.0,
            anti_alias_cutoff_hz: 9000.0,
            anti_alias_order: 3,
            matched_filter_spans: 8,
            matched_filter_window: "blackman".to_string(),
            agc: AgcConfig::default(),
            costas: CostasConfig::default(),
            gardner: GardnerConfig::default(),
            sync: SyncConfig::default(),
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct AgcConfig {
    /// AGC tracking bandwidth (Hz)
    pub bandwidth_hz: f32,
    /// Target output RMS amplitude
    pub target_rms: f32,
}

impl Default for AgcConfig {
    fn default() -> Self {
        AgcConfig {
            bandwidth_hz: 10.0,
            target_rms: 0.001,
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct CostasConfig {
    /// Normalized loop bandwidth
    pub loop_bw: f32,
}

impl Default for CostasConfig {
    fn default() -> Self {
        CostasConfig { loop_bw: 0.05 }
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct GardnerConfig {
    /// Normalized loop bandwidth (~1% of symbol rate)
    pub loop_bw: f32,
}

impl Default for GardnerConfig {
    fn default() -> Self {
        GardnerConfig { loop_bw: 0.01 }
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct SyncConfig {
    /// Max error bits to attempt CRC correction
    pub crc_correction_max_bits: u32,
    /// Consecutive bad blocks before declaring sync loss
    pub loss_threshold: usize,
}

impl Default for SyncConfig {
    fn default() -> Self {
        SyncConfig {
            crc_correction_max_bits: 2,
            loss_threshold: 12,
        }
    }
}

impl RdsConfig {
    /// Load config from a JSON file. Missing fields use defaults.
    pub fn from_file(path: &str) -> Result<RdsConfig, String> {
        let file = std::fs::File::open(path)
            .map_err(|e| format!("Failed to open config {}: {}", path, e))?;
        let config: RdsConfig = serde_json::from_reader(std::io::BufReader::new(file))
            .map_err(|e| format!("Failed to parse config {}: {}", path, e))?;
        Ok(config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = RdsConfig::default();
        assert_eq!(config.baseband_lpf_hz, 4000.0);
        assert_eq!(config.anti_alias_cutoff_hz, 9000.0);
        assert_eq!(config.anti_alias_order, 3);
        assert_eq!(config.matched_filter_spans, 8);
        assert_eq!(config.matched_filter_window, "blackman");
        assert_eq!(config.agc.bandwidth_hz, 10.0);
        assert_eq!(config.costas.loop_bw, 0.05);
        assert_eq!(config.gardner.loop_bw, 0.01);
        assert_eq!(config.sync.loss_threshold, 12);
        assert_eq!(config.sync.crc_correction_max_bits, 2);
    }

    #[test]
    fn test_partial_json() {
        let json = r#"{ "costas": { "loop_bw": 0.08 } }"#;
        let config: RdsConfig = serde_json::from_str(json).unwrap();
        assert_eq!(config.costas.loop_bw, 0.08);
        // Everything else should be default
        assert_eq!(config.baseband_lpf_hz, 4000.0);
        assert_eq!(config.gardner.loop_bw, 0.01);
    }

    #[test]
    fn test_empty_json() {
        let json = r#"{}"#;
        let config: RdsConfig = serde_json::from_str(json).unwrap();
        assert_eq!(config.costas.loop_bw, 0.05);
        assert_eq!(config.agc.target_rms, 0.001);
    }

    #[test]
    fn test_full_json() {
        let json = r#"{
            "baseband_lpf_hz": 5000,
            "anti_alias_cutoff_hz": 8000,
            "anti_alias_order": 4,
            "matched_filter_spans": 12,
            "matched_filter_window": "hamming",
            "agc": { "bandwidth_hz": 20, "target_rms": 0.01 },
            "costas": { "loop_bw": 0.1 },
            "gardner": { "loop_bw": 0.02 },
            "sync": { "crc_correction_max_bits": 3, "loss_threshold": 8 }
        }"#;
        let config: RdsConfig = serde_json::from_str(json).unwrap();
        assert_eq!(config.baseband_lpf_hz, 5000.0);
        assert_eq!(config.anti_alias_order, 4);
        assert_eq!(config.matched_filter_window, "hamming");
        assert_eq!(config.agc.bandwidth_hz, 20.0);
        assert_eq!(config.costas.loop_bw, 0.1);
        assert_eq!(config.gardner.loop_bw, 0.02);
        assert_eq!(config.sync.crc_correction_max_bits, 3);
        assert_eq!(config.sync.loss_threshold, 8);
    }
}
