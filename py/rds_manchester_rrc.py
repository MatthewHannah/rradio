"""
RDS Manchester-shaped RRC Matched Filter — FIR Tap Generator

Generates FIR filter taps for the Manchester-shaped RRC matched filter.
This filter simultaneously matches the RDS pulse shape AND decodes the
Manchester (biphase) encoding, eliminating the need for a separate
Manchester decoder.

The filter is constructed by convolving the standard RRC impulse response
with a Manchester waveform spanning one bit period:
    h_manchester_rrc = h_rrc * manchester_waveform

When applied at 19 kHz (16 samples per bit), the Gardner TED operates at
the bit rate (1187.5 Hz) directly. The output peak polarity encodes the
differential bit value.

Usage:
    python rds_manchester_rrc.py [--fs SAMPLE_RATE] [--spans NUM_CHIP_SPANS]

Output:
    Prints filter taps to stdout in the specified format.
"""

import argparse
import sys
import numpy as np


RDS_BIT_RATE = 1187.5     # Hz
RDS_CHIP_RATE = 2375.0    # Hz
T_D = 1.0 / RDS_BIT_RATE  # bit period
T_CHIP = 1.0 / RDS_CHIP_RATE  # chip period


def rds_rrc_taps(fs: float, num_spans: int = 8,
                 window: str = "blackman") -> np.ndarray:
    """Generate the base RRC impulse response at the given sample rate."""
    n_half = int(round(num_spans * T_CHIP * fs))
    n_taps = 2 * n_half + 1
    n = np.arange(-n_half, n_half + 1)
    u = 4.0 * n / (T_D * fs)

    h = np.sinc(u + 0.5) + np.sinc(u - 0.5)

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
    return h


def rds_manchester_rrc_taps(fs: float, num_spans: int = 8,
                             window: str = "blackman") -> np.ndarray:
    """
    Generate Manchester-shaped RRC matched filter taps.

    The RRC impulse response is convolved with a Manchester waveform
    spanning one bit period (first half +1, second half -1). The result
    is a filter that simultaneously matches the pulse shape and decodes
    Manchester encoding.

    Parameters
    ----------
    fs : float
        Sample rate in Hz. Should be an integer multiple of 1187.5 Hz
        for clean operation (e.g., 19000 Hz = 16 × 1187.5).
    num_spans : int
        Number of chip periods to span on each side for the base RRC.
    window : str
        Window function for the base RRC.

    Returns
    -------
    np.ndarray
        FIR filter taps, normalized to unit energy.
    """
    # Base RRC
    h_rrc = rds_rrc_taps(fs, num_spans, window)

    # Manchester waveform for one bit period: [+1...+1, -1...-1]
    sps = int(round(fs / RDS_BIT_RATE))
    manchester = np.ones(sps)
    manchester[sps // 2:] = -1.0

    # Convolve
    h = np.convolve(h_rrc, manchester)

    # Normalize to unit energy
    energy = np.sqrt(np.sum(h ** 2))
    if energy > 0:
        h /= energy

    return h


def verify_filter(h: np.ndarray, fs: float):
    """Print verification info about the filter to stderr."""
    sps = int(round(fs / RDS_BIT_RATE))
    n_taps = len(h)

    info = []
    info.append(f"Manchester-RRC Matched Filter Summary")
    info.append(f"  Sample rate:       {fs:.0f} Hz")
    info.append(f"  Taps:              {n_taps}")
    info.append(f"  Samples per bit:   {sps}")
    info.append(f"  Samples per chip:  {sps // 2}")
    info.append(f"  Latency:           {(n_taps // 2) / fs * 1e3:.2f} ms")

    # Test: generate a short Manchester-encoded sequence and verify
    # that the filter output peaks align with bit boundaries
    test_diff_bits = [1, 0, 0, 1, 1, 0, 1, 0] * 30
    signal = []
    for d in test_diff_bits:
        if d == 0:
            signal.extend([1.0] * (sps // 2) + [-1.0] * (sps // 2))
        else:
            signal.extend([-1.0] * (sps // 2) + [1.0] * (sps // 2))
    signal = np.array(signal)

    filtered = np.convolve(signal, h, mode='same')

    # Sample at bit centers
    offset = sps // 2
    samples = filtered[offset::sps]

    skip = 10
    check = min(len(samples) - skip, len(test_diff_bits) - skip)
    # positive → diff=1 (determined empirically from filter construction)
    matches = sum(1 for i in range(check)
                  if (samples[i + skip] > 0) == (test_diff_bits[i + skip] == 1))
    accuracy = matches / check * 100.0

    info.append(f"  Bit recovery test: {matches}/{check} = {accuracy:.1f}%")
    info.append(f"  Polarity: positive peak → diff bit 1")

    print("\n".join(info), file=sys.stderr)


def plot_filter(h: np.ndarray, fs: float):
    """Plot frequency response and impulse response with bit-period markers."""
    import matplotlib.pyplot as plt

    sps = int(round(fs / RDS_BIT_RATE))
    n_fft = max(4096, len(h) * 4)

    # --- Frequency response ---
    H = np.fft.rfft(h, n=n_fft)
    freqs = np.fft.rfftfreq(n_fft, 1.0 / fs)
    H_mag = np.abs(H)
    H_mag /= np.max(H_mag)

    # Also show the base RRC for comparison
    h_rrc = rds_rrc_taps(fs, 8, "blackman")
    h_rrc /= np.sqrt(np.sum(h_rrc ** 2))
    H_rrc = np.fft.rfft(h_rrc, n=n_fft)
    H_rrc_mag = np.abs(H_rrc)
    H_rrc_mag /= np.max(H_rrc_mag)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10))

    ax1.plot(freqs, H_rrc_mag, "k--", linewidth=1.0, label="Plain RRC")
    ax1.plot(freqs, H_mag, "b-", linewidth=1.5, label=f"Manchester-RRC ({len(h)} taps)")
    ax1.axvline(RDS_CHIP_RATE, color="r", linestyle=":", alpha=0.5,
                label=f"Chip rate ({RDS_CHIP_RATE:.0f} Hz)")
    ax1.axvline(RDS_BIT_RATE, color="g", linestyle=":", alpha=0.5,
                label=f"Bit rate ({RDS_BIT_RATE:.0f} Hz)")
    ax1.set_xlim(0, min(3 * RDS_CHIP_RATE, fs / 2))
    ax1.set_ylim(-0.05, 1.1)
    ax1.set_xlabel("Frequency (Hz)")
    ax1.set_ylabel("Magnitude (normalized)")
    ax1.set_title("Frequency Response: Manchester-RRC vs. Plain RRC")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # --- Impulse response in bit periods ---
    t_samples = np.arange(len(h)) - len(h) // 2
    t_bits = t_samples / sps

    ax2.plot(t_bits, h, "b-", linewidth=1.0)
    for k in range(-5, 6):
        ax2.axvline(k, color="gray", linestyle=":", alpha=0.3)
    ax2.axhline(0, color="gray", linewidth=0.5)
    ax2.set_xlabel("Time (bit periods)")
    ax2.set_ylabel("Amplitude")
    ax2.set_title("Manchester-RRC Impulse Response")
    ax2.grid(True, alpha=0.3)

    # --- Test signal: filter applied to Manchester-encoded bits ---
    test_diff = [1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1] * 3
    signal = []
    for d in test_diff:
        if d == 0:
            signal.extend([1.0] * (sps // 2) + [-1.0] * (sps // 2))
        else:
            signal.extend([-1.0] * (sps // 2) + [1.0] * (sps // 2))
    signal = np.array(signal)

    filtered = np.convolve(signal, h, mode='same')
    t_sig = np.arange(len(filtered)) / sps

    ax3.plot(t_sig, filtered, "b-", linewidth=0.8, label="Filter output")
    # Mark bit centers
    offset = sps // 2
    for i in range(len(test_diff)):
        idx = i * sps + offset
        if idx < len(filtered):
            color = "green" if test_diff[i] == 1 else "red"
            ax3.plot(t_sig[idx], filtered[idx], "o", color=color,
                     markersize=6, zorder=5)
        ax3.axvline(i, color="gray", linestyle=":", alpha=0.2)

    ax3.axhline(0, color="gray", linewidth=0.5)
    ax3.set_xlabel("Time (bit periods)")
    ax3.set_ylabel("Amplitude")
    ax3.set_title("Filter output on test signal — "
                  "green=diff 1 (positive), red=diff 0 (negative)")
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    fig.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Generate RDS Manchester-shaped RRC matched filter taps")
    parser.add_argument("--fs", type=float, default=19000.0,
                        help="Sample rate in Hz (default: 19000)")
    parser.add_argument("--spans", type=int, default=8,
                        help="Chip periods to span on each side (default: 8)")
    parser.add_argument("--window", type=str, default="blackman",
                        choices=["blackman", "hamming", "hann", "rectangular"],
                        help="Window function (default: blackman)")
    parser.add_argument("--format", type=str, default="text",
                        choices=["text", "rust", "json"],
                        help="Output format (default: text)")
    parser.add_argument("--plot", action="store_true",
                        help="Show frequency response, impulse response, and test signal plots")
    args = parser.parse_args()

    h = rds_manchester_rrc_taps(args.fs, args.spans, args.window)
    verify_filter(h, args.fs)

    if args.plot:
        plot_filter(h, args.fs)

    if args.format == "text":
        print(f"# Manchester-RRC filter: {len(h)} taps, "
              f"fs={args.fs:.0f} Hz, ±{args.spans} chips, "
              f"window={args.window}")
        for tap in h:
            print(f"{tap:.10e}")

    elif args.format == "rust":
        print(f"/// Manchester-shaped RRC matched filter: {len(h)} taps, "
              f"fs={args.fs:.0f} Hz, ±{args.spans} chips, "
              f"window={args.window}")
        print(f"const RDS_MANCHESTER_RRC_TAPS: [f32; {len(h)}] = [")
        for i in range(0, len(h), 4):
            chunk = h[i:min(i + 4, len(h))]
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
            "type": "manchester_rrc",
            "taps": h.tolist(),
        }
        print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
