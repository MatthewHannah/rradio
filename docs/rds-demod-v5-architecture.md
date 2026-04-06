# RDS Demodulator Architecture (v5)

## Overview

The RDS demodulator extracts Radio Data System (RDS) digital data from the 57 kHz subcarrier of an FM broadcast signal. It takes real-valued MPX (multiplex) audio at 171 kHz and outputs decoded RDS groups containing station identification, program type, radiotext, and other metadata.

The architecture decouples carrier recovery from timing recovery, running the Costas loop at a higher sample rate (14250 Hz) than the chip rate (2375 Hz). This provides 6× more carrier corrections per chip, improving lock stability on weak signals and during transients.

## Signal Flow

```
FM MPX (171 kHz, real)
  │
  ├─ Coarse NCO at 57 kHz → complex baseband
  │
  ├─ LPF (1001-tap Blackman, 2500 Hz cutoff) + decimate 12×
  │                                    ↓
  │                          14250 Hz complex baseband (6 SPS)
  │
  ├─ Pre-MF AGC (exponential power averaging)
  │
  ├─ Fine Costas Loop ──────────────────────────┐
  │    ├─ RRC Matched Filter (37 taps, α=0.8)   │
  │    ├─ Phase detector: tanh(I) × Q            │
  │    ├─ PI loop (Bn=30 Hz, ζ=0.707)           │
  │    └─ Fine NCO correction ←──────────────────┘
  │                                    ↓
  │                     Carrier-corrected, MF-shaped (14250 Hz)
  │
  ├─ Post-MF AGC (unit power normalization)
  │
  ├─ Polyphase Gardner Timing Recovery
  │    ├─ 16-arm sinc+Blackman interpolator (8 taps/arm)
  │    ├─ Phase accumulator (fire threshold = 96)
  │    ├─ Gardner TED: e = -Re{y_mid × (y_now* - y_prev*)}
  │    ├─ PI loop (Bn=25 Hz, ζ=1.0, K_TED=0.617)
  │    └─ Outputs on-time chip at ~2375 Hz
  │                                    ↓
  │                          Complex chip (±1 BPSK)
  │
  ├─ Biphase (Manchester) Decoder
  │    ├─ Transition detection: (chip[k] - chip[k-1]) / 2
  │    ├─ Clock alignment via 128-chip energy metric
  │    └─ Outputs data bit at ~1187.5 Hz
  │
  ├─ Delta (Differential) Decoder
  │    └─ XOR with previous bit
  │
  └─ Block Synchronizer
       ├─ CRC-10 syndrome matching (5 offset words)
       ├─ 2-bit error correction in Locked state
       ├─ Exact CRC for sync health (loss detection)
       └─ Outputs decoded RDS groups
```

## Stage Details

### Stage 1: Coarse NCO + Decimation

- **Input**: 171 kHz real-valued MPX
- **NCO**: Fixed 57 kHz, mixes to complex baseband
- **LPF**: 1001-tap Blackman FIR, cutoff 2500 Hz at 171 kHz rate
- **Decimation**: 12× → 14250 Hz output (6 samples per chip)
- **Purpose**: Extract the RDS subcarrier band with sufficient oversampling for the Costas and timing loops

### Stage 2: Fine Carrier Recovery (Costas Loop)

- **Rate**: 14250 Hz (every decimated sample)
- **Matched filter**: Root Raised Cosine, α=0.8, span=±3 chips, 37 taps
- **Phase detector**: `tanh(I) × Q` — soft-decision BPSK detector
  - K_det = tanh(1) ≈ 0.762 at unit signal power
  - tanh saturates large I values, reducing sensitivity to amplitude noise
- **Loop filter**: PI with exact s→z pole mapping gains
  - Bn = 30 Hz, ζ = 0.707 (underdamped)
  - Max frequency tracking: ±100 Hz
- **Fine NCO**: Phase and frequency correction applied to every input sample
- **Why tanh(I)×Q**: More robust than I×Q on weak signals (92.9 MHz). The tanh compresses large I values, preventing noise spikes from corrupting the phase error estimate.

### Stage 3: Polyphase Gardner Timing Recovery

- **Input**: Carrier-corrected, MF-shaped complex samples at 14250 Hz
- **AGC**: Post-MF power normalization (rate=0.01, target RMS=1.0)
- **Interpolator**: 16-arm polyphase lowpass (sinc + Blackman window)
  - 8 taps per arm, 128-tap prototype
  - Provides sub-sample timing precision (1/16 of input sample = 1/96 of chip)
  - Standard lowpass, NOT the RRC — signal is already MF-shaped

**Gardner TED**:
```
y_now  = interp(overshoot)                    — on-time chip sample (peak)
y_mid  = interp(overshoot + period/2)         — mid-chip sample (zero crossing)
y_prev = interp(overshoot + period)           — previous chip sample (peak)

error = -Re{y_mid × (conj(y_now) - conj(y_prev))}
```

- K_TED = 3.70 / SPS ≈ 0.617 (from S-curve analysis)
- The Gardner locks `y_mid` to the zero crossing between consecutive chips
- Single stable equilibrium per chip period (no false lock)
- Data-independent: no hard decisions needed

**PI Loop**:
- Bn = 25 Hz, ζ = 1.0 (critically damped)
- Tracks period in input samples (nominal = 6.0)
- Phase accumulator fires when phase ≥ nfilters × SPS = 96

**Output**: One complex chip sample per fire at ~2375 Hz

### Stage 4: Biphase + Block Sync

Standard RDS decoding chain (unchanged from v4):

- **Biphase decoder**: Energy-based even/odd alignment tracking over 128-chip windows. Each Manchester bit consists of two chips with opposite phase.
- **Delta decoder**: Differential decoding (XOR with previous bit)
- **Block sync**: CRC-10 with syndrome matching against 5 offset words (A, B, C, C', D). 2-bit error correction in Locked state. Exact CRC (no correction) used for sync health monitoring — enables fast re-lock (<0.5s) after timing transients.

## Key Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Input rate | 171 kHz | After 240k→171k resample |
| Baseband rate | 14250 Hz | 6 SPS, decimate 12× |
| Chip rate | 2375 Hz | RDS standard |
| Bit rate | 1187.5 Hz | Manchester: 2 chips/bit |
| RRC α | 0.8 | Excess bandwidth |
| RRC span | ±3 chips | 37 taps at SPS=6 |
| Costas Bn | 30 Hz | Carrier tracking bandwidth |
| Costas ζ | 0.707 | Underdamped |
| Timing Bn | 25 Hz | Symbol timing bandwidth |
| Timing ζ | 1.0 | Critically damped |
| Polyphase arms | 16 | Timing interpolation resolution |
| Taps/arm | 8 | Interpolator filter length |
| Gardner K_TED | 0.617 | S-curve slope per input sample |
| Block sync ECC | 2-bit | Error correction in Locked state |
| Loss threshold | 12 blocks | Consecutive exact-CRC failures to declare loss |

## Performance

Tested on PlutoSDR recordings and live signals:

| Station | Recording | Live (60s) | Redsea |
|---------|-----------|-----------|--------|
| 92.9 (weak) | 534 groups | 682 groups | 580 |
| 94.9 | 329 | 680 | ~329 |
| 97.1 | 321 | — | ~326 |
| 99.7 | 328 | 676 | ~329 |
| 100.5 | 327 | 617 | ~328 |
| 5min live | 3413 | — | 3412 |

Processing speed: ~10× real-time on Apple Silicon (release build).

## Files

| File | Purpose |
|------|---------|
| `src/rds_demod_v2.rs` | v5 demodulator: CoarseNco, FineCostas, PolyphaseGardner |
| `src/rds_demod.rs` | v4 demodulator (fallback via USE_V4=1) |
| `src/rds_taps.rs` | RRC and lowpass filter tap generation |
| `src/biphase.rs` | Biphase (Manchester) decoder |
| `src/rds_block_sync.rs` | CRC block synchronizer |
| `src/rds_config.rs` | Configuration (sync thresholds, ECC bits) |
| `src/symbolsync.rs` | v4 polyphase symbol sync (unused in v5) |
| `src/pi_loop.rs` | v4 PI loop filter (unused in v5) |
