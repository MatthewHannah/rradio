#!/usr/bin/env python3
"""
Gardner TED S-curve and K_TED derivation for the v5 RDS polyphase receiver.

Approach:
  1. Generate 1000 random bits → 2000 Manchester chips
  2. TX RRC pulse-shape at massive oversample (SPS × nfilters × 100)
  3. Sweep tau: pick different starting offsets, decimate to SPS × nfilters
  4. RX RRC matched filter at the decimated rate
  5. Evaluate Gardner error (negated, matching Rust) at perfect period → S(tau)
  6. K_TED = slope of S(tau) / SPS

Parameters match rds_demod.rs:
  SPS = 6, nfilters = 16, RRC α = 0.8, span = 3
"""

import numpy as np
import matplotlib.pyplot as plt

# ── Receiver constants (match rds_demod.rs) ──
SPS = 6
NFILTERS = 16
RRC_ALPHA = 0.8
RRC_SPAN = 3

# Oversample factor — sets tau resolution
OVERSAMPLE = 100
PROTO_SPS = SPS * NFILTERS          # 96 samples/chip at polyphase rate
TX_SPS = PROTO_SPS * OVERSAMPLE     # 9600 samples/chip for "continuous time"


def design_rrc(span, sps, alpha):
    """Root raised cosine filter, energy-normalized."""
    N = span * sps
    n = np.arange(-N // 2, N // 2 + 1, dtype=float)
    t = n / sps
    h = np.zeros_like(t)
    for i, ti in enumerate(t):
        if ti == 0:
            h[i] = 1.0 - alpha + 4.0 * alpha / np.pi
        elif abs(abs(ti) - 1.0 / (4.0 * alpha)) < 1e-8:
            h[i] = (alpha / np.sqrt(2.0)) * (
                (1.0 + 2.0 / np.pi) * np.sin(np.pi / (4.0 * alpha))
                + (1.0 - 2.0 / np.pi) * np.cos(np.pi / (4.0 * alpha))
            )
        else:
            num = np.sin(np.pi * ti * (1 - alpha)) + 4 * alpha * ti * np.cos(np.pi * ti * (1 + alpha))
            den = np.pi * ti * (1 - (4 * alpha * ti) ** 2)
            h[i] = num / den
    h /= np.sqrt(np.sum(h ** 2))
    return h


def manchester_encode(bits):
    """Expand bits to Manchester chips: bit 1 → [+1, -1], bit 0 → [-1, +1]."""
    chips = np.zeros(len(bits) * 2)
    for i, b in enumerate(bits):
        if b:
            chips[2 * i] = +1.0
            chips[2 * i + 1] = -1.0
        else:
            chips[2 * i] = -1.0
            chips[2 * i + 1] = +1.0
    return chips


def main():
    print(f"Gardner S-curve analysis")
    print(f"  SPS={SPS}, nfilters={NFILTERS}, RRC α={RRC_ALPHA}, span={RRC_SPAN}")
    print(f"  Oversample={OVERSAMPLE}, TX_SPS={TX_SPS}, PROTO_SPS={PROTO_SPS}")
    print()

    # 1. Generate 50 random bits → 100 Manchester chips
    rng = np.random.default_rng(42)
    bits = rng.integers(0, 2, size=1000)
    chips = manchester_encode(bits)
    n_chips = len(chips)  # 100
    print(f"  {len(bits)} bits → {n_chips} Manchester chips")

    # 2. TX RRC at massive oversample rate
    h_tx = design_rrc(RRC_SPAN, TX_SPS, RRC_ALPHA)
    up = np.zeros(n_chips * TX_SPS)
    up[::TX_SPS] = chips
    tx = np.convolve(up, h_tx, mode="same")

    # 3. RX RRC matched filter at proto rate
    h_rx = design_rrc(RRC_SPAN, PROTO_SPS, RRC_ALPHA)

    # 4. Sweep tau: offset in TX samples, then decimate by OVERSAMPLE
    #    Pad the TX signal so negative offsets index valid data
    pad = TX_SPS  # one chip of padding
    tx_padded = np.concatenate([np.zeros(pad), tx, np.zeros(pad)])

    n_tau = 2 * OVERSAMPLE + 1  # -1.0 to +1.0 chips
    tau_chips = np.linspace(-1.0, 1.0, n_tau)
    s_curve = np.zeros(n_tau)

    skip = 10  # skip chips at edges for filter transient

    for ti, tau in enumerate(tau_chips):
        # Offset in TX samples from the pad origin
        offset = pad + int(round(tau * TX_SPS))

        # Decimate: pick every OVERSAMPLE-th sample starting at offset
        decimated = tx_padded[offset::OVERSAMPLE]

        # RX matched filter
        mf_out = np.convolve(decimated, h_rx, mode="same")

        # Normalize to ±1 on-time amplitude (AGC equivalent)
        ot = mf_out[skip * PROTO_SPS:(n_chips - skip) * PROTO_SPS:PROTO_SPS]
        rms = np.sqrt(np.mean(ot**2))
        if rms > 1e-12:
            mf_out = mf_out / rms

        # Gardner error at perfect period (PROTO_SPS samples/chip)
        err_sum = 0.0
        count = 0
        for k in range(skip, n_chips - skip):
            idx_now = k * PROTO_SPS
            idx_mid = k * PROTO_SPS - PROTO_SPS // 2
            idx_prev = (k - 1) * PROTO_SPS
            if idx_prev < 0 or idx_now >= len(mf_out):
                continue
            y_now = mf_out[idx_now]
            y_mid = mf_out[idx_mid]
            y_prev = mf_out[idx_prev]
            err_sum += y_mid * (y_now - y_prev)
            count += 1

        s_curve[ti] = -(err_sum / count) if count > 0 else 0  # negate to match Rust/Python receiver

    # 5. Find zero crossings and K_TED
    print("\nZero crossings:")
    zero_crossings = []
    for i in range(1, len(tau_chips)):
        if s_curve[i - 1] * s_curve[i] < 0:
            tau_z = tau_chips[i - 1] - s_curve[i - 1] * (tau_chips[i] - tau_chips[i - 1]) / (s_curve[i] - s_curve[i - 1])
            slope = (s_curve[i] - s_curve[i - 1]) / (tau_chips[i] - tau_chips[i - 1])
            zero_crossings.append((tau_z, slope))
            label = "STABLE" if slope < 0 else "unstable"
            print(f"  τ = {tau_z:+.4f} chips, slope = {slope:.4f}/chip  [{label}]")

    stable = [(t, s) for t, s in zero_crossings if s < 0]
    if stable:
        lock_tau, lock_slope = stable[0]
        k_chip = abs(lock_slope)
        # Arms are unity-DC-gain normalized, so TED output is not scaled by nfilters.
        # PI loop updates once per chip (at R_CHIP rate), K_TED is per input sample.
        k_per_sample = k_chip / SPS
        print(f"\n── K_TED ──")
        print(f"  Lock point: τ = {lock_tau:.4f} chips")
        print(f"  S-curve slope: {lock_slope:.4f} per chip")
        print(f"  K_TED = |slope| / SPS = {k_chip:.4f} / {SPS} = {k_per_sample:.4f}")
        print(f"  Rust TIMING_K_TED = 0.617  (= 3.70 / {SPS})")

    # 6. Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(tau_chips, s_curve, "b-", linewidth=2, label="Gardner S-curve")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5, linestyle=":")

    for tau_z, slope in zero_crossings:
        if slope < 0:
            ax.plot(tau_z, 0, "go", markersize=10, zorder=5,
                    label=f"Lock @ τ={tau_z:.3f}, K_TED={abs(slope)/SPS:.4f}/sample")
            tang_x = np.array([tau_z - 0.15, tau_z + 0.15])
            ax.plot(tang_x, slope * (tang_x - tau_z), "g--", linewidth=1, alpha=0.7)
        else:
            ax.plot(tau_z, 0, "rx", markersize=10, zorder=5)

    ax.set_xlabel("Timing offset τ (chip periods)")
    ax.set_ylabel("Average Gardner TED error")
    ax.set_title(f"Gardner S-curve — SPS={SPS}, RRC α={RRC_ALPHA}, span={RRC_SPAN}, "
                 f"{NFILTERS}-arm polyphase")
    ax.set_xlim(-1.05, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig("gardner_s_curve.png", dpi=150)
    print(f"\nPlot saved to gardner_s_curve.png")
    plt.show()


if __name__ == "__main__":
    main()
