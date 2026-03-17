"""
Output from Claude Opus 4.6

RDS Matched Filter — FIR Tap Generator

Generates FIR filter taps for the RDS matched filter specified by IEC 62106:

    H(f) = cos(π·f·t_d / 4)    for |f| ≤ 2/t_d, 0 otherwise

This is an RRC with α = 1.0 at the biphase chip rate (2375 Hz).

Usage:
    python rds_matched_filter.py [--fs SAMPLE_RATE] [--spans NUM_CHIP_SPANS]

Output:
    Prints filter taps to stdout, one per line.
    Also prints a summary to stderr with filter parameters and verification.
"""

import argparse
import sys
import numpy as np


RDS_BIT_RATE = 1187.5     # Hz
RDS_CHIP_RATE = 2375.0    # Hz (biphase symbol rate = 2 × bit rate)
T_D = 1.0 / RDS_BIT_RATE  # bit period ≈ 842.1 µs
T_CHIP = 1.0 / RDS_CHIP_RATE  # chip period ≈ 421.1 µs


def rds_matched_filter_taps(fs: float, num_spans: int = 8,
                            window: str = "blackman") -> np.ndarray:
    """
    Generate FIR taps for the RDS matched filter.

    Parameters
    ----------
    fs : float
        Sample rate in Hz.
    num_spans : int
        Number of chip periods to span on each side of center.
        More spans = more taps = closer to ideal response.
    window : str
        Window function name: "blackman", "hamming", "hann", or "rectangular".

    Returns
    -------
    np.ndarray
        FIR filter taps, normalized to unit energy (matched filter convention).
    """
    # Step 1-2: Compute the analytical impulse response sampled at fs.
    #
    # h(t) = (2/t_d) · [sinc(4t/t_d + 1/2) + sinc(4t/t_d - 1/2)]
    #
    # We drop the (2/t_d) constant since we normalize later.

    n_half = int(round(num_spans * T_CHIP * fs))
    n_taps = 2 * n_half + 1
    n = np.arange(-n_half, n_half + 1)
    u = 4.0 * n / (T_D * fs)

    h = np.sinc(u + 0.5) + np.sinc(u - 0.5)

    # Step 3-4: Apply window function to control spectral leakage.
    if window == "blackman":
        w = np.blackman(n_taps)
    elif window == "hamming":
        w = np.hamming(n_taps)
    elif window == "hann":
        w = np.hanning(n_taps)
    elif window == "rectangular":
        w = np.ones(n_taps)
    else:
        raise ValueError(f"Unknown window: {window}")

    h *= w

    # Step 5: Normalize to unit energy (matched filter).
    energy = np.sqrt(np.sum(h ** 2))
    if energy > 0:
        h /= energy

    return h


def verify_filter(h: np.ndarray, fs: float):
    """Print verification info about the filter to stderr."""
    n_half = (len(h) - 1) // 2
    n_taps = len(h)

    # Compute frequency response via FFT
    n_fft = max(4096, n_taps * 4)
    H = np.fft.rfft(h, n=n_fft)
    freqs = np.fft.rfftfreq(n_fft, 1.0 / fs)
    H_mag = np.abs(H)
    H_mag /= np.max(H_mag)  # normalize to peak = 1

    # Ideal response for comparison
    H_ideal = np.zeros_like(freqs)
    B = 2.0 / T_D
    mask = freqs <= B
    H_ideal[mask] = np.cos(np.pi * freqs[mask] * T_D / 4.0)
    H_ideal = np.maximum(H_ideal, 0)

    # Max deviation in passband
    passband = freqs <= B
    deviation = np.max(np.abs(H_mag[passband] - H_ideal[passband]))

    # ISI check: convolve filter with itself (combined TX×RX)
    h_combined = np.convolve(h, h)
    center = len(h_combined) // 2
    peak = h_combined[center]

    info = []
    info.append(f"RDS Matched Filter Summary")
    info.append(f"  Sample rate:       {fs:.0f} Hz")
    info.append(f"  Taps:              {n_taps}")
    info.append(f"  Span:              ±{n_half} samples "
                f"(±{n_half / (T_CHIP * fs):.1f} chips)")
    info.append(f"  Latency:           {n_half / fs * 1e3:.2f} ms")
    info.append(f"  Passband deviation: {deviation:.4f} "
                f"(max |H_actual - H_ideal|)")

    info.append(f"  Combined TX×RX ISI at chip boundaries:")
    for k in range(-4, 5):
        idx = center + int(round(k * T_CHIP * fs))
        if 0 <= idx < len(h_combined):
            val = h_combined[idx] / peak
            label = " ◄ peak" if k == 0 else ""
            if k != 0 and abs(val) < 0.01:
                label = " ✓"
            info.append(f"    k={k:+2d}: {val:+.6f}{label}")

    print("\n".join(info), file=sys.stderr)


def plot_filter(h: np.ndarray, fs: float):
    """Plot frequency response comparison and combined TX×RX ISI."""
    import matplotlib.pyplot as plt

    n_half = (len(h) - 1) // 2
    n_fft = max(4096, len(h) * 4)
    B = 2.0 / T_D

    # --- Frequency response: actual vs ideal ---
    H = np.fft.rfft(h, n=n_fft)
    freqs = np.fft.rfftfreq(n_fft, 1.0 / fs)
    H_mag = np.abs(H)
    H_mag /= np.max(H_mag)

    H_ideal = np.zeros_like(freqs)
    mask = freqs <= B
    H_ideal[mask] = np.cos(np.pi * freqs[mask] * T_D / 4.0)
    H_ideal = np.maximum(H_ideal, 0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    ax1.plot(freqs, H_ideal, "k--", linewidth=1.5, label="Ideal: cos(πf·t_d/4)")
    ax1.plot(freqs, H_mag, "b-", linewidth=1.0, label=f"FIR ({len(h)} taps)")
    ax1.axvline(RDS_CHIP_RATE, color="r", linestyle=":", alpha=0.5,
                label=f"Chip rate ({RDS_CHIP_RATE:.0f} Hz)")
    ax1.axvline(RDS_BIT_RATE, color="g", linestyle=":", alpha=0.5,
                label=f"Bit rate ({RDS_BIT_RATE:.0f} Hz)")
    ax1.set_xlim(0, min(2 * B, fs / 2))
    ax1.set_ylim(-0.05, 1.1)
    ax1.set_xlabel("Frequency (Hz)")
    ax1.set_ylabel("Magnitude (normalized)")
    ax1.set_title("Frequency Response: FIR vs. Ideal")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # --- Combined TX×RX with chip boundary markers ---
    h_combined = np.convolve(h, h)
    center = len(h_combined) // 2
    peak = h_combined[center]
    h_combined /= peak

    # Time axis in chip periods
    t_samples = np.arange(len(h_combined)) - center
    t_chips = t_samples / (T_CHIP * fs)

    display_chips = 6
    display_samples = int(display_chips * T_CHIP * fs)
    lo = center - display_samples
    hi = center + display_samples
    ax2.plot(t_chips[lo:hi], h_combined[lo:hi], "b-", linewidth=1.0)

    # Mark chip boundaries and their values
    for k in range(-display_chips, display_chips + 1):
        idx = center + int(round(k * T_CHIP * fs))
        if 0 <= idx < len(h_combined):
            val = h_combined[idx]
            color = "green" if abs(val) < 0.01 else "red"
            if k == 0:
                color = "blue"
            ax2.plot(k, val, "o", color=color, markersize=5, zorder=5)
        ax2.axvline(k, color="gray", linestyle=":", alpha=0.3)

    ax2.axhline(0, color="gray", linewidth=0.5)
    ax2.set_xlim(-display_chips, display_chips)
    ax2.set_xlabel("Time (chip periods)")
    ax2.set_ylabel("Amplitude (normalized)")
    ax2.set_title("Combined TX×RX Response — "
                  "green dots = zero ISI at chip boundaries")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Generate RDS matched filter FIR taps")
    parser.add_argument("--fs", type=float, default=80000.0,
                        help="Sample rate in Hz (default: 80000)")
    parser.add_argument("--spans", type=int, default=8,
                        help="Chip periods to span on each side (default: 8)")
    parser.add_argument("--window", type=str, default="blackman",
                        choices=["blackman", "hamming", "hann", "rectangular"],
                        help="Window function (default: blackman)")
    parser.add_argument("--format", type=str, default="text",
                        choices=["text", "rust", "json"],
                        help="Output format (default: text)")
    parser.add_argument("--plot", action="store_true",
                        help="Show frequency response and ISI plots")
    args = parser.parse_args()

    h = rds_matched_filter_taps(args.fs, args.spans, args.window)
    verify_filter(h, args.fs)

    if args.plot:
        plot_filter(h, args.fs)

    if args.format == "text":
        print(f"# RDS matched filter: {len(h)} taps, "
              f"fs={args.fs:.0f} Hz, ±{args.spans} chips, "
              f"window={args.window}")
        for i, tap in enumerate(h):
            print(f"{tap:+.10e}")

    elif args.format == "rust":
        print(f"/// RDS matched filter: {len(h)} taps, "
              f"fs={args.fs:.0f} Hz, ±{args.spans} chips, "
              f"window={args.window}")
        print(f"const RDS_MATCHED_FILTER: [f32; {len(h)}] = [")
        for i in range(0, len(h), 4):
            chunk = h[i:min(i+4, len(h))]
            vals = ", ".join(f"{v:.7e}" for v in chunk)
            print(f"    {vals},")
        print("];")

    elif args.format == "json":
        import json
        out = {
            "sample_rate": args.fs,
            "spans": args.spans,
            "window": args.window,
            "num_taps": len(h),
            "taps": h.tolist(),
        }
        print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
