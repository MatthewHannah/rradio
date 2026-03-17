Output from Claude Opus 4.6

# Deriving an FIR Filter from a Frequency Response

This document walks through the general DSP procedure for converting a
continuous frequency response specification into a discrete FIR filter,
using the RDS matched filter as a concrete example.

## The Problem

We have a desired frequency response (from IEC 62106):

    H(f) = cos(π·f·t_d / 4)    for |f| ≤ 2/t_d
    H(f) = 0                    otherwise

where t_d = 1/1187.5 s. We need a finite set of discrete filter taps
that approximate this response at a given sample rate.

## Step 1: Inverse Fourier Transform

The ideal (continuous, infinite-length) impulse response is the inverse
Fourier transform of H(f):

    h(t) = ∫ H(f) · e^(j2πft) df

Expand the cosine using Euler's formula:

    cos(πf·t_d/4) = [e^(jπf·t_d/4) + e^(-jπf·t_d/4)] / 2

Substituting and grouping the exponentials:

    cos(πf·t_d/4) · e^(j2πft) = ½[e^(j2πf(t + t_d/8)) + e^(j2πf(t − t_d/8))]

Each term is now a complex exponential integrated over the band [−B, B]
where B = 2/t_d. We apply the standard Fourier pair:

    ∫_{−B}^{B} e^(j2πfτ) df = 2B · sinc(2Bτ)

where sinc(x) = sin(πx)/(πx) is the normalized sinc function.

This gives:

    h(t) = B · [sinc(2B·(t + t_d/8)) + sinc(2B·(t − t_d/8))]

Substituting B = 2/t_d and noting that 2B·(t_d/8) = ½:

    h(t) = (2/t_d) · [sinc(4t/t_d + ½) + sinc(4t/t_d − ½)]

This is the ideal impulse response. It is infinite in extent because the
sinc function's tails never fully reach zero.

## Step 2: Sample at the Discrete Sample Rate

Replace continuous time t with discrete samples n/fs:

    h[n] = sinc(4n/(t_d·fs) + ½) + sinc(4n/(t_d·fs) − ½)

The constant scaling factor (2/t_d) is dropped here since we will
normalize in a later step. This is still an infinite sequence — we
cannot implement it directly.

## Step 3: Truncate to Finite Length

Choose a span: how many chip periods (T_chip = t_d/2) of the sinc tails
to keep on each side of center.

    N_half = round(num_spans × T_chip × fs)
    Keep h[n] for n = −N_half .. +N_half
    Total taps = 2 × N_half + 1

More taps means the FIR more closely approximates the ideal response.
Fewer taps means less computation per sample. Typical trade-offs:

| Span (±chips) | Taps @ 80 kHz | Taps @ 10 kHz |
|----------------|---------------|---------------|
| ±4             | 271           | 35            |
| ±8             | 539           | 69            |
| ±16            | 1079          | 135           |

Truncation is mathematically equivalent to multiplying the infinite h[n]
by a rectangular window. This causes **spectral leakage**: the
frequency response develops ripples and overshoots (Gibbs phenomenon)
near sharp transitions.

## Step 4: Apply a Window Function

To suppress spectral leakage, multiply the truncated taps by a window
function w[n] that tapers smoothly to zero at the edges:

    h_final[n] = h[n] · w[n]

The window trades off main-lobe width (transition sharpness) against
sidelobe suppression (out-of-band rejection):

| Window      | Sidelobe suppression | Main lobe width | Notes                |
|-------------|---------------------|-----------------|----------------------|
| Rectangular | −13 dB              | Narrowest       | No window (worst)    |
| Hamming     | −43 dB              | Moderate        | General purpose      |
| Blackman    | −58 dB              | Wider           | Good sidelobe control|
| Kaiser(β)   | Tunable             | Tunable         | Flexible, β adjusts  |

Blackman is a safe default for most applications. The slight widening of
the transition band is acceptable when the filter already has a gentle
roll-off (as our α = 1.0 RRC does).

## Step 5: Normalize

The final step depends on the filter's intended use:

**Matched filter** (maximize SNR at sampling instants):

    h[n] /= sqrt(Σ h[n]²)     — unit energy normalization

**Spectral shaping** (preserve signal amplitude):

    h[n] /= Σ h[n]            — unit DC gain normalization

For the RDS receiver's matched filter, unit energy normalization is
appropriate.

## Summary

The complete procedure is:

    Frequency response H(f)
        → Inverse Fourier Transform → continuous h(t)
        → Sample at fs             → infinite discrete h[n]
        → Truncate to ±N_half      → finite h[n] (with spectral leakage)
        → Apply window w[n]        → h[n] with controlled sidelobes
        → Normalize                → final FIR taps

The engineering choices are:
1. **Sample rate** — must be ≥ 2× the filter bandwidth (Nyquist)
2. **Span** — more taps = better approximation, more computation
3. **Window** — controls the leakage vs. transition-width trade-off
4. **Normalization** — depends on the application (matched filter vs. shaping)

## Application to RDS

For the RDS matched filter specifically:
- H(f) = cos(πf·t_d/4), which is an RRC with α = 1.0 at chip rate 2375 Hz
- Bandwidth: 2375 Hz → minimum sample rate ~4750 Hz
- Practical sample rate: 80 kHz (or decimate to ~10 kHz first)
- ±8 chip spans is a good balance of accuracy and cost
- Blackman window
- Unit energy normalization for matched filtering

See `py/rds_matched_filter.py` for a working implementation that generates
the filter taps.
