# RDS Demodulator Pipeline — Implementation Guide

This document describes the from-first-principles RDS (Radio Data System) BPSK
demodulator implemented in `rradio`. It provides enough detail to reimplement
the entire pipeline from scratch.

## Signal Overview

RDS is a low-bitrate digital subcarrier embedded in FM broadcasts:

| Parameter | Value |
|-----------|-------|
| Subcarrier frequency | 57,000 Hz |
| Modulation | BPSK (Binary Phase Shift Keying) |
| Symbol rate | 1,187.5 symbols/s |
| Chip rate (Manchester) | 2,375 chips/s |
| Encoding | Manchester (biphase) + differential |

A "chip" is one half of a Manchester symbol. Data bits are encoded as
transitions between consecutive chips, not absolute chip polarity.

## Pipeline Block Diagram

```
FM MPX input (240 kHz real)
    │
    ▼
[1] Rational Resample  240 kHz → 171 kHz (↑57 ↓80)
    │
    ▼
[2] NCO Mix @ 57 kHz   real → complex baseband
    │                   ← freq_adj (ongoing, from Costas)
    │                   ← phase_adj (one-shot, from Costas)
    ▼
[3] Polyphase Decimate  171 kHz → 7,125 Hz (÷24), LPF fc=2500 Hz
    │
    ▼
[4] Symbol Sync (PFB)   7,125 Hz → 2,375 Hz (chip rate)
    │   ├─ Polyphase matched filter (RRC) + derivative
    │   ├─ Post-MF AGC (unit power normalization)
    │   ├─ ML Timing Error Detector
    │   └─ PI loop filter → phase accumulator → arm select
    │
    ▼
[5] Costas Loop          carrier phase/freq recovery
    │   ├─ Phase error: tanh(I)·Q
    │   └─ PI loop filter → NCO corrections
    │
    ▼
[6] Biphase Decode       chip differencing → data symbols
    │
    ▼
[7] Differential Decode  transition → absolute bits
    │
    ▼
[8] Block Sync + CRC     RDS group assembly
```

## Stage Details

### [1] Rational Resampler: 240 kHz → 171 kHz

Converts the FM MPX sample rate to exactly 3× the subcarrier frequency.

- **Interpolation factor**: 57 (up)
- **Decimation factor**: 80 (down)
- **Anti-alias filter**: Lowpass at 171/2 = 85.5 kHz, ~570 taps at the
  interpolated rate (240k × 57 = 13.68 MHz), Blackman window
- **Why 171 kHz**: F_S = 57,000 × 3, giving exactly integer samples per
  subcarrier cycle. This simplifies the NCO and ensures SPS = 3 at the
  decimated rate.

### [2] NCO (Numerically Controlled Oscillator)

Mixes the real input down to complex baseband centered at 0 Hz.

```
phase[n] = phase[n-1] + base_freq + freq_adj + phase_adj
output[n] = input[n] × (cos(phase[n]) − j·sin(phase[n]))
phase_adj = 0  (cleared after first sample of each chip)
```

- **base_freq**: 2π × 57000 / 171000 = 2π/3 rad/sample
- **freq_adj**: ongoing frequency correction from Costas loop (rad/sample
  at F_S rate)
- **phase_adj**: one-shot phase correction from Costas proportional term,
  applied at the chip boundary then zeroed
- **Phase wrapping**: keep phase in [0, 2π)

### [3] Polyphase Decimator: 171 kHz → 7,125 Hz

A rational resampler configured as a pure decimator (↑1 ↓24) with an
integrated lowpass filter.

- **Decimation ratio**: 24
- **Filter**: 1001-tap FIR, Blackman window, cutoff 2500 Hz (at 171 kHz rate)
- **Output rate**: 171,000 / 24 = 7,125 Hz = 3 samples per chip

### [4] Symbol Synchronizer (Polyphase Filterbank)

This is the core timing recovery block. It consumes samples at 7,125 Hz
and produces one chip-rate output (~2,375 Hz) with optimal sampling phase.

#### Architecture: Phase-Accumulation with Polyphase Interpolation

Rather than consuming a fixed number of samples per chip, the synchronizer
pushes **one sample at a time** and maintains a phase accumulator. When the
accumulator wraps past `nfilters`, a chip boundary has occurred and the
matched filter is evaluated.

This follows GNU Radio's `symbol_sync` architecture and handles non-integer
SPS and clock drift naturally.

#### Polyphase Filter Bank

The matched filter is a Root Raised Cosine (RRC) designed at the prototype
rate of `nfilters × SPS`:

| Parameter | Value |
|-----------|-------|
| Roll-off (α) | 0.8 |
| Span | 3 chip periods |
| SPS | 3 |
| nfilters | 32 |
| Prototype rate | 32 × 3 = 96 samples/chip |
| Prototype length | 2 × 3 × 96 + 1 = 577 taps |
| Taps per arm | ceil(577/32) = 19 |

**RRC impulse response** (standard formula):

```
        ⎧ 1 - α + 4α/π                                     t = 0
        ⎪
h(t) =  ⎨ (α/√2)[(1+2/π)sin(π/4α) + (1-2/π)cos(π/4α)]   |t| = 1/(4α)
        ⎪
        ⎩ [sin(πt(1-α)) + 4αt·cos(πt(1+α))] / [πt(1-(4αt)²)]  otherwise
```

where t is in chip periods. Normalize to unit energy: `h /= √(Σh²)`.

**Derivative matched filter**: central difference of the RRC prototype,
then scaled by `nfilters` to compensate for polyphase partition:

```
dh[n] = (h[n+1] - h[n-1]) / 2    (central difference)
dh *= nfilters                     (scale before partition)
```

**Polyphase partition**: split the prototype into `nfilters` arms.
Arm k gets taps `h[k], h[k + nfilters], h[k + 2·nfilters], ...`.
Reverse each arm for convolution via forward dot product.

Both MF and dMF are partitioned identically.

#### Shared Delay Line

A circular buffer of length `taps_per_arm`, duplicated (length × 2) so that
a contiguous slice always covers the full filter span:

```
push(sample):
    circbuf[head] = sample
    circbuf[head + circ_len] = sample    // mirror copy
    head = (head + 1) % circ_len

evaluate(arm):
    mf_out  = dot(filt[arm],  circbuf[head .. head + taps_per_arm])
    dmf_out = dot(dfilt[arm], circbuf[head .. head + taps_per_arm])
```

#### Phase Accumulator

The accumulator runs at the **input sample rate** (7,125 Hz):

```
inst_rate = nfilters / inst_period      // phase advance per sample (~32/3 ≈ 10.67)

for each input sample:
    push sample into delay line
    phase += inst_rate

    if phase >= nfilters:               // chip boundary
        phase -= nfilters
        arm = floor(phase)              // select polyphase arm [0, nfilters)
        mf_out, dmf_out = evaluate(arm)
        ... AGC, TED, loop filter ...
        inst_rate = nfilters / inst_period   // update rate
        output mf_out
```

This naturally produces ~3 input samples per chip. If the actual chip
period drifts to 3.001 samples, `inst_period` adjusts and some chips
consume 4 samples — no explicit skip/slip needed.

#### Post-MF AGC

Applied to the matched filter output **before** the TED, so both K_TED
and K_DET are calibrated for unit-power signals.

**Acquisition phase** (first 100 chips): accumulates power, computes
initial gain, passes samples through unscaled.

**Tracking phase**: exponential power averaging:

```
power = |mf_out|²
avg_power = (1 - rate) × avg_power + rate × power
gain = √(target_power / avg_power)
gain = min(gain, max_gain)

mf_out  *= gain
dmf_out *= gain     // same gain applied to derivative
```

Parameters: `rate = 0.01`, `target_rms = 1.0`, `max_gain = 100,000`.

#### ML Timing Error Detector

Maximum-likelihood TED using the Mengali method:

```
error = Re{conj(mf_out) × dmf_out} / 2
```

This produces an S-curve centered at zero when sampling at the MF peak.
The slope at the zero crossing is K_TED.

**K_TED** depends on the pulse shape, nfilters, and post-AGC normalization.
For our filter: **K_TED = 0.040125 per accumulator unit** (= 1.284 per
input sample). Computed via `py/compute_detector_gains.py`.

#### Timing PI Loop Filter

Second-order PI loop with gains from exact s→z pole mapping (GNU Radio
`clock_tracking_loop` formulation).

**Gain computation** (computed once at initialization):

Given noise bandwidth `Bn` (Hz), damping `ζ`, and detector gain `K_d`:

```
ωn      = 2·Bn / (ζ + 1/(4ζ))          // natural frequency (rad/s)
ωn_norm = ωn / R_CHIP                    // normalized (dimensionless)

ζω = ζ × ωn_norm
k0 = 2 / K_d
k1 = exp(-ζω)

α = k0 × k1 × sinh(ζω)                 // proportional gain
β = k0 × (1 - k1 × (sinh(ζω) + cosₓ(ωd)))   // integral gain
```

where `cosₓ` depends on damping regime:
- ζ < 1: `cos(ωn_norm × √(1 - ζ²))`
- ζ = 1: `1.0`
- ζ > 1: `cosh(ωn_norm × √(ζ² - 1))`

**CRITICAL**: compute in f64, cast result to f32. The subtraction
`1 - k1×(sinh + cos)` suffers catastrophic cancellation in f32 for
small ωn_norm.

**Loop update** (once per chip):

```
avg_period  += β × error                  // integral arm (input samples)
avg_period   = clamp(avg_period, min, max)
inst_period  = avg_period + α × error     // proportional + integral
inst_rate    = nfilters / inst_period      // update phase advance rate
```

**Parameters**:

| | Value |
|---|---|
| Bn | 25 Hz (was 50, tunable) |
| ζ | 1.0 (critically damped) |
| K_TED (per sample) | 1.284 |
| Period limits | 3.0 ± 1.0 input samples |

### [5] Costas Carrier Recovery Loop

Runs once per chip (2,375 Hz) on the AGC-normalized MF output. Corrects
the NCO which runs at 171 kHz — a multi-rate configuration.

#### Phase Error Detector

Soft-decision BPSK Costas discriminator:

```
error = tanh(I) × Q
```

where `I = Re(mf_out)`, `Q = Im(mf_out)`.

Linearized gain: `K_det = A × tanh(A)`. For unit amplitude: **K_det = tanh(1)
≈ 0.761594**.

The `tanh` limiter reduces noise sensitivity at low SNR compared to `I × Q`.

#### PI Loop Filter

Same gain computation as the timing loop, but with different parameters:

| | Value |
|---|---|
| Bn | 30 Hz |
| ζ | 0.707 (Butterworth) |
| K_det | 0.761594 |
| Integrator limits | ±2π × 10 / 2375 (±10 Hz tracking) |

**Loop update** (once per chip, in chip-rate units):

```
freq     += β × error                     // rad/chip (integrator)
freq      = clamp(freq, ±max_freq)
pi_output = freq + α × error              // rad/chip (full PI)
```

#### Multi-Rate NCO Feedback

The Costas PI produces corrections in rad/chip. The NCO runs at 171 kHz
(72 samples per chip). Two feedback paths:

```
nco.freq_adj  = freq / 72                // ongoing (zero-order hold between chips)
nco.phase_adj = (pi_output - freq)       // = α × error (one-shot at chip boundary)
```

The frequency correction is applied every NCO sample. The phase correction
is applied once then cleared to zero.

### [6] Biphase (Manchester) Decode

RDS uses Manchester encoding: data is in transitions, not absolute polarity.

```
every 2 chips:
    symbol = chip[k] - chip[k-1]
```

Output: one symbol at 1,187.5 Hz. Magnitude ≈ ±2 for correct decoding.

Manchester naturally resolves the 180° BPSK phase ambiguity — the Costas
loop can lock at 0° or 180° and the bit stream is identical.

### [7] Differential Decode

Converts transition-based encoding to absolute bits:

```
bit = current_transition XOR previous_transition
```

### [8] Block Synchronization + CRC

RDS groups are 104 bits: 4 blocks × 26 bits (16 data + 10 check).
Block sync searches for the syndrome pattern across all 4 block positions.
CRC-10 provides error detection (and limited correction).

## Loop Filter Gain Derivation

Both the timing and carrier loops use identical gain math, derived from
matching the discrete closed-loop transfer function to desired s-plane
poles mapped via z = e^{sT}.

The closed-loop transfer function (from `clock_tracking_loop.h`):

```
H(z) = K_d(α+β)z⁻¹ × [1 - α/(α+β) z⁻¹]
       ────────────────────────────────────────────────────
       1 - 2(1 - K_d(α+β)/2) z⁻¹ + (1 - K_d·α) z⁻²
```

Matching the denominator to z-plane poles from `s = -ζωn ± jωd` via `z = e^{sT}`:

```
K_d·α     = 1 - e^{-2ζωnT}
K_d·(α+β) = 2(1 - cosₓ(ωdT) × e^{-ζωnT})
```

Solving for α and β gives the formulas in the implementation section above.

The key property: if K_d is correct, the gains produce **exact** pole
placement at the specified ωn and ζ, regardless of the detector type.

## Detector Gain Calibration

Two detector gains must be determined for the loop filters to work correctly.

### K_TED (Timing Error Detector)

The slope of the TED's S-curve at the zero crossing, measured per accumulator
unit, at unit signal power (post-AGC).

**Computed from the polyphase filter** using `py/compute_detector_gains.py`:

1. Evaluate each arm's center tap to get MF(arm) and dMF(arm)
2. Normalize so MF peak = 1.0 (post-AGC condition)
3. S-curve at each arm: `S[k] = dMF[k] × MF[k]`
4. K_TED = |S[peak+1] - S[peak-1]| / 2 (central difference at zero crossing)

Current value: **0.040125** per acc unit (× 32 nfilters = **1.284** per sample).

### K_DET (Costas Phase Error Detector)

Analytical for the `tanh(I)·Q` detector.

For BPSK at amplitude A with phase error θ:

```
e = tanh(A·cos θ) × A·sin θ
```

Linearized around θ = 0: `e ≈ K_det × θ` where `K_det = A × tanh(A)`.

At unit amplitude: **K_det = tanh(1) ≈ 0.761594**.

## Numerical Considerations

### f32 Catastrophic Cancellation

The gain formula computes `β = (2/K) × (1 - k1 × (sinh + cos))`. When
ωn_norm < 0.001, `k1 × (sinh + cos)` rounds to exactly 1.0 in f32,
making β = 0 (no integral gain). **Always compute gains in f64.**

### Phase Wrapping

The NCO phase must be kept in [0, 2π). Rust's `%` operator on negative
floats preserves sign, so use explicit while loops:

```
while phase >= 2π: phase -= 2π
while phase < 0:   phase += 2π
```

## Parameter Summary

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Input sample rate | F_S | 171,000 | Hz |
| Subcarrier | F_C | 57,000 | Hz |
| Chip rate | R_CHIP | 2,375 | Hz |
| Symbol rate | R_SYM | 1,187.5 | Hz |
| Pre-decimation | | 24× | |
| Decimated rate | F_DEC | 7,125 | Hz |
| Samples per chip | SPS | 3 | |
| F_S per chip | TOTAL_DEC | 72 | |
| RRC roll-off | α | 0.8 | |
| RRC span | | 3 | chips |
| Polyphase arms | nfilters | 32 | |
| Costas Bn | | 30 | Hz |
| Costas ζ | | 0.707 | |
| Costas K_det | | 0.761594 | |
| Costas max freq | | ±10 | Hz |
| Timing Bn | | 25 | Hz |
| Timing ζ | | 1.0 | |
| Timing K_TED/acc | | 0.040125 | |
| Timing K_TED/sample | | 1.284 | |
| AGC rate | | 0.01 | |
| AGC acquisition | | 100 | chips |

## Performance

On 94.9 FM (clean signal, 30s recording):
- **329 groups decoded**, 0% BLER
- First group at ~115 ms
- Matches reference liquid-dsp/redsea pipeline (330 groups)

## File Map

| File | Role |
|------|------|
| `src/rds_demod.rs` | Top-level pipeline: NCO + downsample + SymbolSync + Costas |
| `src/symbolsync.rs` | PFB symbol sync: filter bank, AGC, TED, timing loop |
| `src/pi_loop.rs` | PI loop filter with exact s→z gain computation |
| `src/rds_taps.rs` | RRC and lowpass filter generation |
| `src/resample.rs` | Rational resampler (polyphase) |
| `src/biphase.rs` | Manchester + differential decode |
| `src/rds_block_sync.rs` | Block sync, CRC, group decode |
| `py/compute_detector_gains.py` | K_TED and K_DET computation tool |
