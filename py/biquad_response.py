#!/usr/bin/env python3
"""Biquad filter frequency and phase response viewer.

Computes the frequency response H(z) from biquad coefficients using the same
design formulas as src/biquad.rs: Robert Bristow-Johnson's Audio EQ Cookbook.

Usage:
    python py/biquad_response.py                      # default: lowpass 4kHz @ 240kHz
    python py/biquad_response.py --type lowpass --fs 240000 --cutoff 4000 --q 0.707
    python py/biquad_response.py --type highpass --fs 48000 --cutoff 1000 --q 0.707
    python py/biquad_response.py --cascade 2 --type lowpass --fs 240000 --cutoff 4000
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqz, sosfreqz


def biquad_lowpass(fs, cutoff, q):
    """Lowpass biquad coefficients (matches src/biquad.rs).
    Returns SOS row: [b0, b1, b2, 1, a1, a2] (a0 normalized to 1)."""
    w0 = 2 * np.pi * cutoff / fs
    alpha = np.sin(w0) / (2 * q)
    b0 = (1 - np.cos(w0)) / 2
    b1 = 1 - np.cos(w0)
    b2 = (1 - np.cos(w0)) / 2
    a0 = 1 + alpha
    a1 = -2 * np.cos(w0)
    a2 = 1 - alpha
    return np.array([b0/a0, b1/a0, b2/a0, 1.0, a1/a0, a2/a0])


def biquad_highpass(fs, cutoff, q):
    """Highpass biquad coefficients (matches src/biquad.rs).
    Returns SOS row: [b0, b1, b2, 1, a1, a2]."""
    w0 = 2 * np.pi * cutoff / fs
    alpha = np.sin(w0) / (2 * q)
    b0 = (1 + np.cos(w0)) / 2
    b1 = -(1 + np.cos(w0))
    b2 = (1 + np.cos(w0)) / 2
    a0 = 1 + alpha
    a1 = -2 * np.cos(w0)
    a2 = 1 - alpha
    return np.array([b0/a0, b1/a0, b2/a0, 1.0, a1/a0, a2/a0])


def compute_group_delay(sos, fs, n_points=4096):
    """Compute group delay for SOS filter sections using scipy's freqz.
    
    Group delay = -d(unwrapped_phase)/d(omega), computed via the ratio
    of the derivative polynomial to the original polynomial (exact method).
    """
    # Build cascaded b, a from SOS for group delay computation
    b_total = np.array([1.0])
    a_total = np.array([1.0])
    for section in sos:
        b_total = np.convolve(b_total, section[:3])
        a_total = np.convolve(a_total, section[3:])

    # Exact group delay: -Re(d/dz[ln H(z)] * z) evaluated on unit circle
    # For H(z) = B(z)/A(z):
    #   gd = Re{ z * [B'(z)/B(z) - A'(z)/A(z)] }
    # where B'(z) = d/dz B(z)
    b_deriv = np.polyder(b_total[::-1])[::-1]  # derivative of B(z)
    a_deriv = np.polyder(a_total[::-1])[::-1]

    freqs = np.linspace(0, fs/2, n_points, endpoint=False)
    w = 2 * np.pi * freqs / fs
    z = np.exp(1j * w)

    def polyval_z(coeffs, z_vals):
        """Evaluate polynomial in z^-1: c0 + c1*z^-1 + c2*z^-2 + ..."""
        result = np.zeros_like(z_vals)
        for i, c in enumerate(coeffs):
            result += c * z_vals ** (-i)
        return result

    B = polyval_z(b_total, z)
    A = polyval_z(a_total, z)
    # Derivative polynomials: d/dz of sum(c_k * z^-k) = sum(-k * c_k * z^(-k-1))
    dB = np.zeros_like(z)
    for k in range(len(b_total)):
        dB += (-k) * b_total[k] * z ** (-k - 1)
    dA = np.zeros_like(z)
    for k in range(len(a_total)):
        dA += (-k) * a_total[k] * z ** (-k - 1)

    # Group delay in samples = -Re(z * (dB/B - dA/A))
    # The negative sign comes from group delay = -d(phase)/d(omega)
    gd_samples = -np.real(z * (dB / B - dA / A))

    return freqs, gd_samples / fs  # return in seconds


def plot_response(freqs, H, gd_freqs, gd_seconds, fs, title):
    mag_db = 20 * np.log10(np.abs(H) + 1e-20)
    phase_deg = np.degrees(np.unwrap(np.angle(H)))
    gd_ms = gd_seconds * 1000

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    fig.suptitle(title, fontsize=13)

    # Magnitude
    ax = axes[0]
    ax.plot(freqs, mag_db, 'b', linewidth=1.2)
    ax.set_ylabel('Magnitude (dB)')
    ax.set_ylim(max(mag_db.min(), -80), mag_db.max() + 3)
    ax.grid(True, alpha=0.3)
    ax.axhline(-3, color='r', linestyle='--', alpha=0.5, label='-3 dB')
    ax.legend(fontsize=9)

    # Phase
    ax = axes[1]
    ax.plot(freqs, phase_deg, 'g', linewidth=1.2)
    ax.set_ylabel('Phase (degrees)')
    ax.grid(True, alpha=0.3)

    # Group delay
    ax = axes[2]
    ax.plot(gd_freqs, gd_ms, 'm', linewidth=1.2)
    ax.set_ylabel('Group Delay (ms)')
    ax.set_xlabel('Frequency (Hz)')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def main():
    parser = argparse.ArgumentParser(description='Biquad filter response viewer')
    parser.add_argument('--type', choices=['lowpass', 'highpass'], default='lowpass',
                        help='Filter type (default: lowpass)')
    parser.add_argument('--fs', type=float, default=240000,
                        help='Sample rate in Hz (default: 240000)')
    parser.add_argument('--cutoff', type=float, default=4000,
                        help='Cutoff frequency in Hz (default: 4000)')
    parser.add_argument('--q', type=float, default=0.707,
                        help='Q factor (default: 0.707 = Butterworth)')
    parser.add_argument('--cascade', type=int, default=1,
                        help='Number of cascaded stages (default: 1)')
    parser.add_argument('--save', type=str, default=None,
                        help='Save plot to file instead of showing')
    args = parser.parse_args()

    if args.save:
        import matplotlib
        matplotlib.use('Agg')

    if args.type == 'lowpass':
        section = biquad_lowpass(args.fs, args.cutoff, args.q)
    else:
        section = biquad_highpass(args.fs, args.cutoff, args.q)

    # Build SOS matrix (n_stages × 6)
    sos = np.tile(section, (args.cascade, 1))

    # Frequency response via scipy
    freqs, H = sosfreqz(sos, worN=4096, fs=args.fs)

    # Group delay via exact polynomial method
    gd_freqs, gd_seconds = compute_group_delay(sos, args.fs)

    cascade_str = f' × {args.cascade} cascade' if args.cascade > 1 else ''
    title = (f'{args.type.capitalize()} biquad{cascade_str}\n'
             f'fs={args.fs:.0f} Hz, fc={args.cutoff:.0f} Hz, Q={args.q:.3f}')

    # Print key metrics
    mag_db = 20 * np.log10(np.abs(H) + 1e-20)
    phase_rad = np.unwrap(np.angle(H))
    b, a = section[:3], section[3:]
    print(f'Filter: {args.type}, fs={args.fs} Hz, cutoff={args.cutoff} Hz, Q={args.q}, stages={args.cascade}')
    print(f'Coefficients: b={b.tolist()}, a={a.tolist()}')
    print()

    for f_check in [100, 500, 1000, 2375, 2400, 4000, 5000, 10000, 19000]:
        if f_check <= args.fs / 2:
            idx = np.argmin(np.abs(freqs - f_check))
            gd_idx = np.argmin(np.abs(gd_freqs - f_check))
            gd_ms = gd_seconds[gd_idx] * 1000
            gd_samp = gd_seconds[gd_idx] * args.fs
            print(f'  {f_check:>6} Hz: {mag_db[idx]:7.1f} dB, phase = {np.degrees(phase_rad[idx]):8.1f}°, group delay = {gd_ms:.4f} ms ({gd_samp:.1f} samp)')

    plot_response(freqs, H, gd_freqs, gd_seconds, args.fs, title)

    if args.save:
        plt.savefig(args.save, dpi=150, bbox_inches='tight')
        print(f'\nSaved to {args.save}')
    else:
        plt.show()


if __name__ == '__main__':
    main()


if __name__ == '__main__':
    main()
