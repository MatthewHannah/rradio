# RDS Pipeline Architecture

## Overview

The RDS (Radio Data System) pipeline extracts a 1187.5 bps digital data stream
embedded in the 57 kHz subcarrier of FM broadcast signals. The pipeline design
closely follows [Redsea](https://github.com/windytan/redsea) by Oona Räisänen,
which uses [liquid-dsp](https://github.com/jgaeddert/liquid-dsp) for its signal
processing primitives.

On a 28.9-second recording, the pipeline decodes 329 groups at 0.0% block error
rate, running at 9.6× real-time on an Apple M1. It sustains real-time operation
on live RTL-SDR input.

## Signal Flow

```
IQ samples @ 2.4 MHz (from SDR or recording)
    │
    ├─ ÷2 downsample → 1.2 MHz
    ├─ FM discriminator → baseband MPX @ 240 kHz
    ├─ Stereo decoder (19 kHz pilot PLL) → left/right audio
    │
    └─ MPX tee (real f32 @ 240 kHz) ──→ RDS pipeline thread
                                              │
                                              ▼
                                    ┌─────────────────┐
                                    │  RDS Pipeline    │
                                    └─────────────────┘
```

## RDS Pipeline Stages

### 1. Polyphase Rational Resampler (240 kHz → 171 kHz)

The input MPX signal at 240 kHz is resampled to 171 kHz, matching Redsea's
internal processing rate. This rate was chosen so that an integer decimation by
24 yields exactly 7125 Hz = 3 × 2375 Hz (three samples per chip).

- **Interpolation factor (L):** 57
- **Decimation factor (M):** 80
- **Ratio:** 57/80 = 0.7125 (240 × 0.7125 = 171 kHz)
- **Prototype filter:** 570-tap Blackman-windowed lowpass at 85.5 kHz cutoff,
  designed at L × 240 kHz = 13.68 MHz, scaled by L
- **Implementation:** Polyphase decomposition into 57 arms of 10 taps each

### 2. NCO Downconversion (57 kHz → baseband)

A numerically controlled oscillator locked to 57 kHz mixes the real MPX signal
down to complex baseband. The NCO incorporates a phase-locked loop that receives
phase error feedback from the BPSK modem (stage 7) for fine carrier tracking.

- **Nominal frequency:** 57,000 Hz
- **Mix operation:** `baseband = sample × exp(-j × phase)`
- **PLL bandwidth:** 0.03 Hz (extremely narrow — fine phase correction only)
- **PLL gains (liquid-dsp formula):** α = bw (frequency), β = √bw (phase)
- **PLL feedback multiplier:** ×12

The NCO steps its phase on every sample at 171 kHz. The PLL is updated at the
chip rate (2375 Hz) with phase error from the BPSK modem.

### 3. FIR Lowpass Filter (2400 Hz cutoff)

A lowpass filter isolates the RDS signal bandwidth around DC. The RDS BPSK
signal occupies roughly ±2375 Hz around the 57 kHz subcarrier (now at DC after
mixing).

- **Taps:** 255 (Blackman window)
- **Cutoff:** 2400 Hz at 171 kHz sample rate
- **Scale factor:** 2 × (2400/171000) ≈ 0.028 (matching liquid-dsp's convention)
- **Optimization:** Samples are pushed into the delay line at 171 kHz, but the
  convolution is only computed at decimation points (every 24th sample). This
  reduces computation from 43.6M to 1.8M multiply-accumulates per second.

### 4. Decimation (171 kHz → 7125 Hz)

Every 24th sample is kept, reducing the rate from 171 kHz to 7125 Hz. This gives
exactly 3 samples per chip at the 2375 Hz chip rate. The FIR lowpass from stage 3
provides anti-alias protection.

### 5. Automatic Gain Control

The AGC normalizes signal amplitude using liquid-dsp's algorithm: an exponential
power tracker with logarithmic gain adjustment.

- **Normalized bandwidth:** 500/171000 ≈ 0.003
- **Initial gain:** 0.08
- **Update rule:** `gain *= exp(-0.5 × α × ln(output_power))`
- **Maximum gain:** 10⁶ (120 dB)

### 6. Polyphase Filterbank Symbol Synchronizer (SymSync)

The symbol synchronizer is the heart of the demodulator. It simultaneously
performs matched filtering, fractional timing interpolation, and clock recovery
using a polyphase filterbank architecture.

- **Filter type:** Root Raised Cosine (RRC)
- **Roll-off factor (β):** 0.8
- **Samples per symbol (k):** 3
- **Filter semi-length (m):** 3 symbols
- **Number of polyphase filter banks:** 32
- **Timing loop bandwidth:** 2200/171000 ≈ 0.013
- **Output rate:** 1 symbol per chip period (2375 Hz)

**How it works:**

The RRC prototype filter is decomposed into 32 polyphase arms (sub-filters).
A derivative version of the same filter is also decomposed into 32 arms. As input
samples arrive, a fractional timing phase τ sweeps from 0 to 1, selecting which
polyphase arm to evaluate. When τ crosses 1.0, a symbol boundary has been reached
and one output symbol is produced.

The timing error detector uses the Mengali method:

```
error = Re(conj(matched_filter_output) × derivative_filter_output)
```

This error drives an IIR loop filter that adjusts the rate at which τ advances,
keeping the symbol timing locked to the transmitter's clock.

### 7. BPSK Modem (Phase Error Extraction)

A BPSK (Binary Phase Shift Keying) modem makes a hard decision on each symbol
and extracts the phase error between the received symbol and the nearest
constellation point.

- **Decision:** `bit = (Re(symbol) ≥ 0)`
- **Phase error:** `atan2(±Im, ±Re)` relative to nearest constellation point
  (+1 or −1 on the real axis)

The modem's bit decision is *not used* for data — it only provides the phase
error that feeds back to the NCO PLL (stage 2) with a ×12 gain multiplier.
This decision-directed loop keeps the carrier phase aligned so that the BPSK
constellation sits on the real axis.

### 8. Biphase (Manchester) Decoder (2375 Hz → 1187.5 Hz)

RDS uses biphase coding (differential Manchester) where each data bit is
represented by two chips with opposite phase. The biphase decoder subtracts
consecutive chip symbols and thresholds the real component.

- **Biphase symbol:** `(current_chip − previous_chip) / 2`
- **Bit decision:** `Re(biphase_symbol) ≥ 0`
- **Alignment tracking:** Every 128 chips, the decoder compares even-index and
  odd-index energy sums to determine which phase alignment produces stronger
  symbols. This self-corrects the even/odd ambiguity within ~50 ms.

### 9. Differential Decoder

RDS uses differential encoding: the transmitted bit is the XOR of the current
and previous differentially-encoded bits. The decoder simply XORs consecutive
biphase output bits.

```
data_bit = current_biphase_bit ⊕ previous_biphase_bit
```

### 10. Block Synchronization and Group Assembly

The decoded bit stream is organized into 26-bit blocks (16 data + 10 check bits).
Four consecutive blocks form a group (104 bits total). A CRC-10 syndrome check
identifies block boundaries and corrects up to 2-bit errors.

- **Sync state machine:** Searching → Tentative → Locked
- **Lock criteria:** 5 consecutive good CRC checks
- **Loss threshold:** 12 consecutive bad CRC checks before returning to search
- **Error correction:** Syndrome lookup table with 367 correctable patterns (≤2 bits)

## Performance

### Decoding Quality

Tested across 5 recorded FM stations (~30 seconds each):

| Recording | Groups | BLER |
|-----------|--------|------|
| 94.9 MHz  | 329    | 0.0% |
| 96.1 MHz  | 329    | 0.0% |
| 100.5 MHz | 330    | 0.0% |
| 97.1 MHz  | 320    | 0.0% |
| 99.7 MHz  | 318    | 0.0% |
| **Total** | **1626** | — |

For comparison, Redsea (C++ with liquid-dsp) decodes 1633 groups on the same
recordings — our pipeline achieves 99.6% of Redsea's performance.

### Throughput

| Mode | Time | Speed |
|------|------|-------|
| Recording (28.9s) | 3.0s | 9.6× real-time |
| Live RTL-SDR (60s) | 60.1s | 1.0× real-time |

### Compute Budget

| Stage | Rate | Taps | MAC/s |
|-------|------|------|-------|
| Polyphase resample | 171k out | 10/arm | 1.7M |
| NCO (trig) | 171k | — | 171k evals |
| FIR LPF (push only) | 171k push, 7125 exec | 255 | 1.8M |
| AGC | 7125 | — | negligible |
| SymSync | 7125 | 19/arm | 0.14M |
| **Total** | | | **~4M** |

## Key Design Decisions

**Why 171 kHz?** This rate divides evenly: 171000 / 24 / 3 = 2375 Hz (chip rate).
Integer samples per chip eliminates interpolation error in the symbol synchronizer.

**Why SymSync instead of Gardner?** The polyphase filterbank approach integrates
matched filtering, fractional timing interpolation, and clock recovery into a
single step. The derivative-matched-filter timing error detector (Mengali method)
is more robust than the Gardner TED at low samples-per-symbol (3 samp/chip).

**Why decision-directed PLL instead of Costas loop?** The Costas loop operates on
every sample and needs multiple samples per symbol to track. At 3 samp/chip, it
struggles to acquire. The decision-directed approach updates the NCO PLL once per
chip (2375 Hz) with a clean phase error from the modem, after the SymSync has
already produced well-timed symbols.

**Why push/execute FIR?** The FIR lowpass runs at 171 kHz but only one out of
every 24 outputs is needed. By splitting into push (delay line update) and execute
(convolution), we avoid 23/24 of the convolution work — a 24× reduction in the
dominant compute cost.

## Source Files

| File | Purpose |
|------|---------|
| `main.rs` (`rds_pipeline`) | Pipeline orchestration, sample loop |
| `symsync.rs` | Polyphase filterbank symbol synchronizer |
| `nco.rs` | NCO with PLL for carrier tracking |
| `agc.rs` | Automatic gain control (liquid-dsp compatible) |
| `fir.rs` | FIR filter with push/execute optimization |
| `biphase.rs` | Manchester decoder and differential decoder |
| `psk_modem.rs` | BPSK hard decision and phase error |
| `polyphase.rs` | Polyphase rational resampler |
| `rds_block_sync.rs` | CRC-10, block sync state machine, group decoder |
| `rds_taps.rs` | Filter tap generation (lowpass, RRC, Manchester-RRC) |
| `rds_config.rs` | Configuration structs (JSON-configurable) |
