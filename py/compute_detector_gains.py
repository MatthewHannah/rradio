#!/usr/bin/env python3
"""Compute K_TED and K_DET detector gains for the RDS demodulator.

These values are used as constants in src/rds_demod.rs. Rerun this script
whenever you change the pulse shape (RRC alpha, span), nfilters, or the
phase error detector type.

Usage:
    python py/compute_detector_gains.py [--alpha 0.8] [--span 8] [--sps 3] [--nfilters 12]

Output:
    Prints Rust constants ready to paste into rds_demod.rs.
"""

import argparse
import numpy as np


# ── RRC filter design ────────────────────────────────────────

def design_rrc(span, sps, alpha):
    """Root raised cosine filter (unit energy normalized)."""
    N = span * sps
    n = np.arange(-N // 2, N // 2 + 1)
    t = n / sps

    h = np.zeros_like(t, dtype=float)
    for i, ti in enumerate(t):
        if ti == 0.0:
            h[i] = 1.0 - alpha + 4.0 * alpha / np.pi
        elif abs(abs(ti) - 1.0 / (4.0 * alpha)) < 1e-10:
            h[i] = (alpha / np.sqrt(2.0)) * (
                (1.0 + 2.0 / np.pi) * np.sin(np.pi / (4.0 * alpha))
                + (1.0 - 2.0 / np.pi) * np.cos(np.pi / (4.0 * alpha))
            )
        else:
            num = (np.sin(np.pi * ti * (1 - alpha))
                   + 4 * alpha * ti * np.cos(np.pi * ti * (1 + alpha)))
            den = np.pi * ti * (1 - (4 * alpha * ti) ** 2)
            h[i] = num / den

    h /= np.sqrt(np.sum(h**2))
    return h


def design_drrc(span, sps, alpha):
    """Derivative RRC via central difference."""
    h = design_rrc(span, sps, alpha)
    dh = np.zeros_like(h)
    dh[1:-1] = (h[2:] - h[:-2]) / 2.0
    dh[0] = h[1] - h[0]
    dh[-1] = h[-1] - h[-2]
    return dh


def polyphase_partition(h, M):
    """Split prototype h into M polyphase arms. Shape: (M, taps_per_arm)."""
    pad = (M - len(h) % M) % M
    h_padded = np.concatenate([h, np.zeros(pad)])
    n_taps = len(h_padded) // M
    return h_padded.reshape(n_taps, M).T


# ── K_TED: timing error detector gain ───────────────────────

def compute_k_ted(h_arms, dh_arms, nfilters):
    """Compute K_TED per accumulator unit from polyphase filter arms.

    For the ML TED  e = Re{dMF · MF*}, the S-curve is:
        S(τ) = p'(τ) · p(τ)    where p = pulse autocorrelation

    K_TED = |dS/dτ| at the zero crossing, per accumulator unit,
    evaluated at unit signal power (post-AGC normalization).
    """
    arm_len = h_arms.shape[1]
    center = arm_len // 2

    mf_vals = np.array([h_arms[k, center] for k in range(nfilters)])
    dmf_vals = np.array([dh_arms[k, center] for k in range(nfilters)])

    peak_arm = np.argmax(np.abs(mf_vals))

    # Normalize to unit MF peak (post-AGC)
    peak_val = abs(mf_vals[peak_arm])
    if peak_val > 1e-20:
        mf_vals = mf_vals / peak_val
        dmf_vals = dmf_vals / peak_val

    s_curve = dmf_vals * mf_vals

    k_p = (peak_arm + 1) % nfilters
    k_m = (peak_arm - 1) % nfilters
    k_ted = abs(s_curve[k_p] - s_curve[k_m]) / 2.0

    return k_ted, s_curve, mf_vals, peak_arm


# ── K_DET: Costas phase error detector gain ──────────────────

def compute_k_det_tanh(amplitude=1.0):
    """K_det for tanh(I)·Q detector.

    For BPSK at amplitude A:  K_det = A · tanh(A)
    At unit amplitude (A=1):  K_det = tanh(1) ≈ 0.7616
    """
    return amplitude * np.tanh(amplitude)


def compute_k_det_iq(amplitude=1.0):
    """K_det for I·Q detector.

    For BPSK at amplitude A:  K_det = A²
    At unit amplitude (A=1):  K_det = 1.0
    """
    return amplitude ** 2


# ── Main ─────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--alpha', type=float, default=0.8, help='RRC roll-off (default: 0.8)')
    parser.add_argument('--span', type=int, default=8, help='RRC span in chips (default: 8)')
    parser.add_argument('--sps', type=int, default=3, help='Samples per chip at decimated rate (default: 3)')
    parser.add_argument('--nfilters', type=int, default=12, help='Polyphase arms (default: 12)')
    parser.add_argument('--detector', choices=['tanh', 'iq'], default='tanh',
                        help='Costas phase error detector type (default: tanh)')
    args = parser.parse_args()

    proto_sps = args.nfilters * args.sps

    # Build polyphase filter bank (same as Rust symbolsync.rs)
    h_proto = design_rrc(args.span, proto_sps, args.alpha)
    dh_proto = design_drrc(args.span, proto_sps, args.alpha)
    dh_proto *= args.nfilters  # compensate for polyphase derivative scaling

    h_arms = polyphase_partition(h_proto, args.nfilters)
    dh_arms = polyphase_partition(dh_proto, args.nfilters)

    # Compute K_TED
    k_ted, s_curve, mf_vals, peak_arm = compute_k_ted(h_arms, dh_arms, args.nfilters)

    # Compute K_DET
    if args.detector == 'tanh':
        k_det = compute_k_det_tanh(1.0)
        det_desc = "tanh(I)·Q at unit power"
    else:
        k_det = compute_k_det_iq(1.0)
        det_desc = "I·Q at unit power"

    k_ted_per_sample = k_ted * args.nfilters

    # Print results
    print(f"Filter: RRC α={args.alpha}, span={args.span}, SPS={args.sps}, nfilters={args.nfilters}")
    print(f"  Prototype: {len(h_proto)} taps at {proto_sps} sps → {h_arms.shape[1]} taps/arm")
    print(f"  MF peak at arm {peak_arm}, value={mf_vals[peak_arm]:.4f} (normalized to 1.0)")
    print()
    print(f"K_TED = {k_ted:.6f}  (per acc unit, post-AGC)")
    print(f"K_TED = {k_ted_per_sample:.6f}  (per input sample = K_TED_acc × nfilters)")
    print(f"K_DET = {k_det:.6f}  ({det_desc})")
    print()
    print(f"S-curve at each arm (normalized):")
    for k in range(args.nfilters):
        marker = " <-- peak" if k == peak_arm else ""
        print(f"  arm {k:2d}: MF={mf_vals[k]:+.4f}  S={s_curve[k]:+.6f}{marker}")
    print()
    print("// Paste into src/rds_demod.rs:")
    print(f"const COSTAS_K_DET: f32 = {k_det:.6f};  // {det_desc}")
    print(f"const TIMING_K_TED: f32 = {k_ted:.6f};  // per acc unit, RRC α={args.alpha} span={args.span} nfilters={args.nfilters}")
    print(f"const TIMING_NFILTERS: usize = {args.nfilters};")
    print(f"// K_TED per input sample = TIMING_K_TED * TIMING_NFILTERS = {k_ted_per_sample:.6f}")


if __name__ == "__main__":
    main()
