# RDS/RBDS Bitstream Extraction — DSP Chain Discussion

## What the Transmitter Does

Understanding the TX chain tells us exactly what we need to undo:

```
Data bits (1187.5 bps)
    ↓ Differential encode: d[n] = b[n] ⊕ d[n-1]
Differentially-encoded bits (1187.5 bps)
    ↓ Biphase/Manchester encode: '0' → [+1, −1], '1' → [−1, +1]
Biphase chip stream (2375 chips/sec, NRZ ±1)
    ↓ Pulse shaping filter (≈raised cosine, α ≈ 1.0)
Shaped baseband waveform (bandwidth ≈ ±2.4 kHz)
    ↓ DSBSC modulate onto 57 kHz subcarrier
RDS subcarrier signal (54.6–59.4 kHz)
    ↓ Sum into FM composite baseband
    ↓ FM modulate onto RF carrier
```

Key numbers:
- **Bit rate**: 1187.5 bps = 19 kHz / 16 = 57 kHz / 48
- **Biphase chip rate**: 2375 Hz = 19 kHz / 8 = 57 kHz / 24
- **Subcarrier**: 57 kHz = 3 × 19 kHz pilot (phase-coherent)
- **Spectral occupancy**: ≈ ±2.4 kHz around 57 kHz

## The Receiver DSP Chain

### Stage 1: Mix to Baseband  ✅ EXISTS

Multiply the FM-demodulated signal by the phase-locked 57 kHz carrier:

```
baseband(t) = FM_demod(t) × cos(2π·57000·t + φ)
            = RDS_baseband(t) + [component at 114 kHz]
```

**Current code** (wideband_fm_audio.rs:113): `rds_out.out * s`

The PLL (pll.rs) locks to the 19 kHz pilot at ÷3, outputs at 57 kHz.
The phase-coherent relationship means we get the correct carrier phase
(up to a possible 180° ambiguity — handled later by differential coding).

**Sample rate at this point: 240 kHz** → 101 samples per biphase chip.

### Stage 2: Anti-Image Lowpass  ✅ EXISTS (needs modification)

The mixing creates a component at 114 kHz. At 240 kHz sample rate the
Nyquist is 120 kHz, so the 114 kHz image is representable but close to
the edge. The existing cascade of two 3 kHz biquad lowpass filters
removes this and limits noise bandwidth.

**Current code** (lines 113-114): `rds_filt` and `rds_filt2` (both 3 kHz LP)

**Problem**: Line 115 hard-thresholds to binary: `if rds > 0.0 { 1.0 } else { 0.0 }`
This destroys all amplitude/timing information. For proper decoding, we need
the **continuous-valued** signal from here onward.

The 3 kHz cutoff is reasonable (RDS occupies ±2.4 kHz), though the biquad
roll-off is gentle. This is adequate as a preliminary filter — the matched
filter in the next stage provides the precise spectral shaping.

### Stage 3: Matched Filter (IEC 62106 RRC)  ❌ NEEDS IMPLEMENTATION

**What it does**: The matched filter maximizes SNR at the optimal sampling
instants. It correlates the incoming signal with the expected pulse shape,
producing clean peaks at each chip center whose polarity encodes the chip
value.

**What the "impulse shapes" are**: After baseband mixing and lowpass
filtering, what you see are overlapping shaped pulses — each biphase chip
(+1 or −1) convolved with the TX pulse shaping filter. The matched filter
turns these into clean, well-separated peaks that can be sampled and
thresholded.

**The standard's filter** (IEC 62106):

The RDS standard specifies a shaping filter in the frequency domain:

    H(f) = cos(π·f·t_d / 4)    for 0 ≤ f ≤ 2/t_d
    H(f) = 0                    otherwise

where t_d = 1/1187.5 s (the bit period).

This is a smooth lowpass: unity at DC, monotonically rolling off to zero
at f = 2/t_d = 2375 Hz (the chip rate). It is **literally the root raised
cosine with α = 1.0** at the biphase chip rate, written in closed form:

    Textbook RRC at α=1:  H(f) = cos(πf / (2·R_s))
    With R_s = 2375 Hz:   H(f) = cos(πf / 4750) = cos(πf·t_d / 4)  ✓

**Combined TX × RX response**:

    H_total(f) = cos²(π·f·t_d/4) = [1 + cos(π·f·t_d/2)] / 2

This is a **raised cosine** with α = 1.0 at symbol rate 2375 Hz:
- Nyquist bandwidth: 1187.5 Hz
- Total bandwidth: 2375 Hz
- **Zero ISI at every chip boundary** (< 0.3% residual in simulation)

**Impulse response** (time domain):

    h(t) = (2/t_d) · [sinc(4t/t_d + ½) + sinc(4t/t_d − ½)]

This is a sinc-like pulse peaking at t = 0 with decaying sidelobes:
- h(0) = 1.0 (peak)
- h(±T_chip) ≈ −0.067 (small — the individual RRC has mild ISI,
  but combined TX×RX cancels it to zero)
- h(±2·T_chip) ≈ −0.016
- Decays roughly as 1/t beyond the main lobe

**Implementation as an FIR filter**:

Since the impulse response has infinite extent (sinc-like), truncate
and window for a practical FIR:

    For each tap index n = −N_half .. +N_half:
        t_n = n / fs
        u_n = 4 · t_n / t_d
        h[n] = sinc(u_n + 0.5) + sinc(u_n − 0.5)

    Apply window (Blackman recommended) and normalize.

At fs = 240 kHz with ±8 chip spans:
- N_half = round(8 × T_chip × fs) = 809
- Total taps: **1619**
- Decimate ÷25 to 9600 Hz first (→ **65 taps** at 4.04 samp/chip)

**A simpler alternative**: Some decoders skip the formal matched filter
entirely and just use the lowpass-filtered signal for clock recovery and
thresholding. This works for strong signals but loses ~2–3 dB of noise
margin compared to proper matched filtering.

**Output**: A continuous-valued signal where peaks at the correct sampling
instants represent biphase chip values. Not yet "decoded" — just optimally
prepared for sampling.

### Stage 4: Clock Recovery  ❌ NEEDS IMPLEMENTATION

**Goal**: Find the optimal sampling instants within the matched filter
output — i.e., recover the 2375 Hz chip clock phase.

**Critical insight**: The chip clock is coherent with the subcarrier:
**57000 / 24 = 2375 Hz exactly**. Since we have a locked 57 kHz PLL,
we could derive the chip clock by dividing by 24. BUT we still need to
find the correct phase offset (1 of 24 possibilities). Options:

#### Option A: Carrier-Derived Clock (divide 57 kHz by 24)
- Divide the PLL's phase accumulator by 24 to get 2375 Hz
- Try all 24 phase offsets; pick the one that produces valid Manchester
  patterns and eventually valid RDS sync words
- **Pro**: Simplest, no loop dynamics, immediately frequency-locked
- **Con**: Must search 24 phases; no automatic tracking of multipath/drift

#### Option B: Squaring / Spectral Line Method
- Square (or take |x|) of the baseband signal → creates spectral line at 2375 Hz
  (guaranteed by Manchester transitions at every bit midpoint)
- Bandpass filter or PLL to lock onto the 2375 Hz line
- Extract phase from the locked oscillator → sampling instants
- **Pro**: Robust, self-tracking, exploits Manchester structure
- **Con**: More complex; squaring doubles noise

#### Option C: Gardner Timing Error Detector
- Classic feedback approach at 2 samples/symbol
- Interpolate the MF output to exactly 2× chip rate (4750 Hz)
- Error signal: e[n] = (x[n] − x[n−2]) · x[n−1]
  (difference of adjacent symbol samples × midpoint sample)
- Drive a timing NCO that adjusts sampling phase
- **Pro**: Well-understood, data-independent, standard DSP
- **Con**: Requires interpolator; loop bandwidth tuning

#### Option D: Mueller & Müller TED
- Works at 1 sample/symbol
- Error: e[n] = d̂[n−1]·x[n] − d̂[n]·x[n−1]  (decision-directed)
- **Pro**: Simple, works at symbol rate
- **Con**: Needs reliable decisions (chicken-and-egg at startup)

**Recommendation**: Start with **Option A** (carrier-derived) for simplicity
since we already have the locked 57 kHz carrier. The frequency is exact;
we only need to search for phase. Graduate to **Option C** (Gardner) if
robustness or tracking performance requires it.

### Stage 5: Symbol Sampling & Decision  ❌ NEEDS IMPLEMENTATION

At each clock recovery instant, sample the matched filter output.
Threshold at zero:
- sample > 0 → chip = +1
- sample ≤ 0 → chip = −1

Output: 2375 chips/sec stream of ±1 values.

(Soft decisions could be preserved for later error correction if desired,
but hard decisions are standard for RDS.)

### Stage 6: Manchester/Biphase Decoding  ❌ NEEDS IMPLEMENTATION

Group chips into pairs at the **bit** rate (1187.5 Hz). Each pair maps to
one differentially-encoded bit:

| Chip pair | Diff bit |
|-----------|----------|
| [+1, −1]  | 0        |
| [−1, +1]  | 1        |
| [+1, +1]  | INVALID  |
| [−1, −1]  | INVALID  |

**Bit boundary ambiguity**: With a 2375 Hz chip clock, there are **two**
possible pairings (even/odd alignment). To resolve:
1. Try both alignments
2. The correct one produces no (or very few) invalid pairs
3. Confirm with RDS sync word detection (Stage 8)

If using carrier-derived clock (Option A), the 24-phase search naturally
subsumes this — half the phases correspond to each alignment.

### Stage 7: Differential Decoding  ❌ NEEDS IMPLEMENTATION

The differential encoding is: `d[n] = b[n] ⊕ d[n-1]` at the transmitter.

To decode: `b[n] = d[n] ⊕ d[n-1]`

This is a simple XOR of adjacent bits.

**Key benefit**: Differential coding resolves the 57 kHz carrier **phase
ambiguity**. If our PLL locks 180° off, ALL chip values are negated. This
negation propagates through Manchester decoding as bit inversion. But
differential decoding cancels it: `(d[n]⊕1) ⊕ (d[n-1]⊕1) = d[n] ⊕ d[n-1]`.
So we get the correct data bits regardless of carrier phase.

Output: **1187.5 bps data bit stream** — the raw RDS bitstream.

## Summary of the Chain

```
FM-demodulated signal @ 240 kHz
    ↓ × cos(57 kHz) from PLL           [EXISTS]
    ↓ Cascaded 3 kHz LPF               [EXISTS, remove hard threshold]
Baseband RDS (continuous, ~101 samp/chip)
    ↓ Decimate ÷25 → 9600 Hz            [simple integer, 4.04 samp/chip, 0.4% ISI]
    ↓ RRC matched filter (cos(πf·t_d/4), α=1 at 2375 Hz chip rate, 65 taps)
    ↓ Clock recovery (2375 Hz chip clock)
    ↓ Sample at chip instants, threshold → ±1
Biphase chip stream (2375 chips/sec)
    ↓ Pair chips → Manchester decode
Diff-encoded bits (1187.5 bps)
    ↓ XOR adjacent bits
Data bits (1187.5 bps)
    ↓ [future: block sync, CRC, group decode]
```

## Architecture: Separating RDS from Audio

### Previous Architecture (before this work)

The pipeline was a single iterator chain running on one thread at 240 kHz:

```
IQ source (thread 1) → buffer → audio_pipeline (thread 2) → buffer → playback (thread 3)
```

Inside `audio_pipeline`, everything was composed as lazy iterators at a
single sample rate, with the RDS signal hard-thresholded to binary and
only used in a spy callback for visualization.

### The Problem

The RDS decoder needs to run at a different rate than the audio path:
- Audio: 240 kHz → ÷5 → 48 kHz output
- RDS:   240 kHz → ÷25 → 9600 Hz → matched filter → clock recovery → bits

A pure iterator chain can only run at one rate. We need to split.

### Implemented Architecture

The pipeline is split into four threads connected by ring buffers:

```
Thread 1: IQ source → iq_buf
Thread 2: signal_pipeline(iq_buf → FM demod → wfm_audio) → tee
              ├→ audio_buf  Vec<(f32, f32)> @ 240 kHz  (blocking send)
              └→ rds_buf    Vec<f32>        @ 240 kHz  (try_get — non-blocking, drops if full)
Thread 3: rds_pipeline(rds_buf → downsample ÷25 → 9600 Hz → [future: decode])
Main:     audio_buf → downsample(5) → interleave → playback/wav @ 48 kHz
```

The `signal_pipeline` function runs the DSP chain (IQ filtering, FM demod,
stereo/RDS extraction via `wfm_audio`) and tees the output to two typed
ring buffers. The audio buffer receives `(f32, f32)` stereo pairs; the RDS
buffer receives `f32` baseband samples. The audio consumer does lazy
`downsample(5).interleave()` at the edge. The RDS consumer is a separate
function (`rds_pipeline`) running on its own thread.

The RDS buffer uses `SendBuf::try_get()` (non-blocking) so a slow RDS
consumer never stalls the audio path — samples are dropped if the buffer
is full.

### What Changed in the Codebase

1. **`wideband_fm_audio.rs`**: ✅ Removed the hard threshold on line 115.
   RDS field now carries the continuous-valued baseband signal.

2. **`buffer.rs`**: ✅ Added `try_get()` method to `SendBuf` — non-blocking
   buffer acquisition using `mpsc::try_recv()`. Unit tested (5 new tests).

3. **`main.rs`**: ✅ Replaced `audio_pipeline` with `signal_pipeline` that
   tees to two typed buffers. Extracted `rds_pipeline` as a separate
   function. `run()` wires 4 threads: IQ source, signal pipeline, RDS
   consumer, audio consumer (main thread). Audio consumer does lazy
   `downsample(5).interleave()` at the boundary.

4. **New module (future): `rds/`** (or `rds_decoder.rs`): Will contain:
   - FIR matched filter (65 taps at 9600 Hz)
   - Clock recovery state machine
   - Manchester decoder
   - Differential decoder
   - Block sync, CRC, group decode

### Alternative: Rational Resample to 9500 Hz

For optimal chip alignment (ISI < 0.03%), resample to 9500 Hz = 4 × chip
rate. From 240 kHz: ratio = 9500/240000 = 19/480.

This is unnecessary given ÷25 → 9600 Hz already achieves 0.4% ISI with
simple integer decimation. The 9600 Hz rate (4.04 samp/chip) is close
enough to the ideal 9500 Hz (4.00 samp/chip) that the difference is
negligible for practical decoding.

**Recommendation**: Decimate ÷25 to 9600 Hz. Simple, 65-tap matched
filter, 0.4% ISI — well within tolerance for reliable RDS decoding.

## Addressing the Original Questions

**"Matched filter — to decode impulses? or bit levels?"**
Neither exactly. The matched filter doesn't "decode" — it **optimally shapes**
the signal for sampling. It's a correlator that maximizes the signal-to-noise
ratio at each chip center. Think of it as turning the blurry overlapping
impulse shapes you currently see into sharp, well-defined peaks whose
polarities correspond to chip values. The "decoding" happens when you
sample those peaks and threshold them.

**"I don't know where to go from there"**
Your proposed pipeline is correct. The full sequence is:
1. Matched filter → makes peaks clean and separable
2. Clock recovery → finds exactly where to sample those peaks
3. Sample + threshold → ±1 chip stream
4. Manchester decode → differential bit stream
5. Differential decode → data bits

## What Comes After the Bits (Future)

Once you have the 1187.5 bps data stream, the digital layer is:

1. **Block synchronization**: RDS data is organized into blocks of 26 bits
   (16 data + 10 check). Each block has an offset word (A/B/C/C'/D) XORed
   into the check bits. Slide a 26-bit window, compute CRC syndrome, match
   against known offsets to find block boundaries.

2. **Error detection/correction**: CRC-10 with generator polynomial
   g(x) = x¹⁰ + x⁸ + x⁷ + x⁵ + x⁴ + x³ + 1. Can detect multiple-bit
   errors and correct single-bit burst errors.

3. **Group assembly**: 4 consecutive blocks (A+B+C+D = 104 bits) form a
   group. Group type (0A–15B) is encoded in block B.

4. **Application decode**: Extract PI code, program name, radiotext,
   clock/time, AF lists, etc. from group contents.

## Key Design Decisions

1. **Matched filter**: Standard's RRC (cos(πf·t_d/4)) — correct per spec.
   At 9600 Hz after ÷25 decimation: **65 taps**, 0.4% ISI.

2. **Clock recovery**: Carrier-derived (simplest) vs. Gardner (most robust)?
   Recommend starting carrier-derived.

3. **Sample rate for RDS path**: ✅ Decimate ÷25 from 240 kHz → 9600 Hz.
   Simple integer decimation, 4.04 samp/chip.

4. **Architecture**: ✅ Four-thread tee architecture. `signal_pipeline`
   tees wfm_audio output to separate audio and RDS ring buffers. Audio
   consumer does lazy `downsample(5).interleave()`. RDS consumer runs
   on its own thread (`rds_pipeline`). Non-blocking `try_get` on RDS
   buffer prevents stalling audio.

5. **RDS baseband signal**: ✅ Hard threshold removed. Continuous-valued
   RDS baseband flows through the tee to the RDS consumer.

## Next Steps

1. Implement FIR matched filter (65 taps at 9600 Hz) in `rds_pipeline`
2. Implement clock recovery (start with carrier-derived ÷24 from 57 kHz PLL)
3. Implement Manchester decoder and differential decoder
4. Implement block synchronization and CRC
5. Implement RDS group decode (PI code, program name, radiotext, etc.)
