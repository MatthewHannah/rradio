# RDS BLER Optimization Experiment Plan

## Goal

Minimize block error rate (BLER) across a diverse set of FM recordings by
systematically tuning the RDS DSP pipeline parameters.

## Data Collection

### SigMF Recording Capability

Add a `--record <path>` flag to the binary that writes raw IQ samples to
a SigMF file pair (`<path>.sigmf-data` + `<path>.sigmf-meta`) while
simultaneously running the normal audio/RDS pipeline. This lets us
capture test data from live signals for reproducible experiments.

Record from multiple stations at varying signal strengths:
- Strong station (e.g., 94.9 — good RDS)
- Moderate station (e.g., 96.1 — marginal RDS)
- Weak station (find one with poor but nonzero RDS)
- Different SDR hardware if available (RTL-SDR, PlutoSDR)
- Different times of day (signal propagation varies)

Each recording should be 30–60 seconds for meaningful BLER statistics.

## Metrics Output

Each decoded RDS group emits a tagged JSON line to stderr:

```
RDSMETRIC {"t":0.523,"bler":0.05,"groups":4,"state":"locked","pi":"199D"}
```

Fields:
- `t`: elapsed time in seconds since pipeline start
- `bler`: rolling BLER (over last 200 blocks)
- `groups`: total groups decoded so far
- `state`: "searching", "tentative", or "locked"
- `pi`: most recent PI code (hex string)

This format allows:
- `grep RDSMETRIC` to extract only metrics from mixed stderr output
- JSON parsing per line for analysis
- Plotting BLER over time for each experiment run

At exit, emit a summary line:

```
RDSSUMMARY {"total_blocks":878,"passed":741,"lifetime_bler":0.156,"final_rolling_bler":0.12,"groups":195,"duration":21.5}
```

## Tunable Parameters

These are exposed via a JSON configuration file (`--rds-config <path>`):

```json
{
  "anti_alias_cutoff_hz": 9000,
  "anti_alias_order": 3,
  "baseband_lpf_hz": 4000,
  "matched_filter_spans": 8,
  "matched_filter_window": "blackman",
  "costas_loop_bw": 0.05,
  "agc_bandwidth_hz": 10,
  "agc_target_rms": 0.001,
  "gardner_loop_bw": 0.01,
  "crc_correction_max_bits": 2,
  "sync_loss_threshold": 5
}
```

If no config file is provided, use current hardcoded defaults.

### Parameter Descriptions

| Parameter | Default | Range | Effect |
|-----------|---------|-------|--------|
| `anti_alias_cutoff_hz` | 9000 | 4000–9500 | Pre-decimation anti-alias filter. Lower = less noise but attenuates signal edges |
| `anti_alias_order` | 3 | 1–6 | Number of cascaded biquad stages. More = steeper roll-off |
| `baseband_lpf_hz` | 4000 | 2500–6000 | LPF in wfm_audio after 57 kHz mix. Affects both signal and noise bandwidth |
| `matched_filter_spans` | 8 | 4–16 | Chip periods spanned by the matched filter. More = better frequency response but more latency |
| `matched_filter_window` | "blackman" | blackman/hamming/hann | Window function for the filter. Affects sidelobe vs. main lobe tradeoff |
| `costas_loop_bw` | 0.05 | 0.01–0.2 | Costas loop tracking bandwidth. Higher = faster lock, more noise |
| `agc_bandwidth_hz` | 10 | 1–100 | AGC tracking speed. Too fast distorts per-chip amplitude, too slow doesn't track |
| `agc_target_rms` | 0.001 | 0.0001–0.1 | Target output amplitude. Must match Gardner loop gain assumptions |
| `gardner_loop_bw` | 0.01 | 0.005–0.1 | Gardner timing recovery bandwidth. Higher = faster lock, more jitter |
| `crc_correction_max_bits` | 2 | 0–5 | Max error bits to attempt correction. Higher recovers more but risks false corrections |
| `sync_loss_threshold` | 5 | 3–15 | Consecutive bad blocks before declaring sync loss |

## Implementation Order

### Step 1: SigMF Recording ✅

Implemented. `--record <path>` writes IQ to SigMF v1.2.0 files.
`--duration <seconds>` auto-stops after the specified time.
Three test recordings captured (94.9, 96.1, 100.5 MHz).

### Step 2: BLER Metrics Output ✅

Implemented. `--rds-metrics` emits `RDSMETRIC` tagged JSON per group
and `RDSSUMMARY` at exit. Grep-friendly for experiment scripting.

### Step 3: JSON Pipeline Config

Add `--rds-config <path>` flag. Parse the JSON at startup into an
`RdsConfig` struct, clone it to both the signal pipeline thread and
the RDS pipeline thread.

#### Config Structure

A top-level `RdsConfig` struct with nested per-component configs:

```rust
// src/rds_config.rs

#[derive(Clone, Deserialize)]
#[serde(default)]
pub struct RdsConfig {
    pub baseband_lpf_hz: f32,       // used by signal_pipeline (wfm_audio)
    pub anti_alias_cutoff_hz: f32,  // used by rds_pipeline
    pub anti_alias_order: usize,
    pub matched_filter_spans: usize,
    pub matched_filter_window: String,
    pub agc: AgcConfig,
    pub costas: CostasConfig,
    pub gardner: GardnerConfig,
    pub sync: SyncConfig,
}

#[derive(Clone, Deserialize)]
#[serde(default)]
pub struct AgcConfig {
    pub bandwidth_hz: f32,  // default 10
    pub target_rms: f32,    // default 0.001
}

#[derive(Clone, Deserialize)]
#[serde(default)]
pub struct CostasConfig {
    pub loop_bw: f32,  // default 0.05
}

#[derive(Clone, Deserialize)]
#[serde(default)]
pub struct GardnerConfig {
    pub loop_bw: f32,  // default 0.01
}

#[derive(Clone, Deserialize)]
#[serde(default)]
pub struct SyncConfig {
    pub crc_correction_max_bits: u32,  // default 2
    pub loss_threshold: usize,         // default 5
}
```

Each component module (agc.rs, costas.rs, gardner_clock_recovery.rs,
rds_block_sync.rs) accepts its own config in its constructor. The
module owns its config type — RdsConfig just nests them.

`#[serde(default)]` on every struct means a partial or empty JSON
uses defaults for missing fields. An empty `{}` gives all defaults.

#### JSON File Format

```json
{
  "baseband_lpf_hz": 4000,
  "anti_alias_cutoff_hz": 9000,
  "anti_alias_order": 3,
  "matched_filter_spans": 8,
  "matched_filter_window": "blackman",
  "agc": { "bandwidth_hz": 10, "target_rms": 0.001 },
  "costas": { "loop_bw": 0.05 },
  "gardner": { "loop_bw": 0.01 },
  "sync": { "crc_correction_max_bits": 2, "loss_threshold": 5 }
}
```

Any field can be omitted to use its default. A minimal config tuning
just the Costas loop:

```json
{ "costas": { "loop_bw": 0.08 } }
```

#### Threading

`RdsConfig` is `Clone`. In `run()`:

```rust
let rds_config = load_config(config_path); // or RdsConfig::default()
let rds_config2 = rds_config.clone();

// Signal pipeline thread gets the config for baseband_lpf_hz
signal_pipeline(..., rds_config);

// RDS pipeline thread gets the config for everything else
rds_pipeline(..., rds_config2);
```

Both threads own their own copy. No shared state, no synchronization.

#### Matched Filter Tap Generation (Runtime)

Port the tap generation from `py/rds_manchester_rrc.py` to Rust so
taps can be computed at startup from `matched_filter_spans` and
`matched_filter_window`. This is ~30 lines of math:

1. Compute RRC: `h[n] = sinc(u+0.5) + sinc(u-0.5)` for each tap
2. Apply window (Blackman/Hamming/Hann)
3. Convolve with Manchester waveform `[+1...+1, -1...-1]`
4. Normalize to unit energy

The hardcoded `RDS_MANCHESTER_RRC_TAPS` const stays as a fallback
for the default config, but runtime generation is used when config
specifies non-default filter parameters.

#### Chain Construction

Each component's constructor takes its config directly (no dual
`new()` / `with_config()` paths). Tests use `Config::default()`:

```rust
.dsp_filter(Agc::new(rds_fs, &config.agc))
.clock_recover(bit_rate, rds_fs, &config.gardner)
.costas_demod(&config.costas)
.rds_block_sync(&config.sync, debug)
```

### Step 4: Experiment Runner (Python)

A script (`py/rds_experiment.py`) that:

1. Takes a list of SigMF recordings and a parameter search space
2. For each parameter combination:
   a. Writes a JSON config file
   b. Runs `cargo run --release -- sigmf <recording> --rds-metrics --rds-config <config>`
   c. Captures `RDSMETRIC` lines from stderr
   d. Computes objective metric (e.g., mean rolling BLER after 5s settling)
3. Reports best configuration and Pareto front (if multi-objective)

Search strategies:
- **Grid search**: exhaustive over a coarse grid, good for 2–3 parameters
- **Random search**: sample uniformly, surprisingly effective for >3 parameters
- **Bayesian optimization**: `scikit-optimize` or similar, best for expensive evaluations

Given that each run takes ~20s per recording and we have ~11 parameters,
a full grid search is impractical. Random search with ~200 evaluations
per recording is a reasonable starting point.

## Expected Outcomes

- Identify which parameters have the most impact on BLER
- Find a configuration that works well across all recordings
- Quantify the tradeoff between parameters (e.g., Costas BW vs. Gardner BW)
- Potentially discover that some parameters are insensitive and can be removed from the search
