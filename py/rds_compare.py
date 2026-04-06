#!/usr/bin/env python3
"""
RDS decoder comparison: PySDR reference implementation vs. rradio.

Runs the PySDR RDS decoder (adapted from pysdr.org/content/rds.html)
against our SigMF recordings and compares group counts with rradio.

Usage:
    python py/rds_compare.py res/recordings/*.sigmf-meta
"""

import json
import os
import subprocess
import sys
import numpy as np
from scipy.signal import resample_poly, firwin

RRADIO_BINARY = "target/release/rradio"


def load_sigmf(meta_path):
    """Load IQ samples from a SigMF recording."""
    with open(meta_path) as f:
        meta = json.load(f)
    sample_rate = meta["global"]["core:sample_rate"]
    freq = meta["captures"][0].get("core:frequency", 0)
    data_path = meta_path.replace(".sigmf-meta", ".sigmf-data")
    samples = np.fromfile(data_path, dtype=np.complex64)
    return samples, sample_rate, freq


def pysdr_rds_decode(samples, sample_rate):
    """
    PySDR RDS decoder, adapted from pysdr.org/content/rds.html.
    Returns the number of decoded groups.
    """
    x = samples

    # If sample rate is much higher than 250 kHz, decimate IQ first
    # PySDR's pipeline expects ~250 kHz input
    target_iq_rate = 250e3
    if sample_rate > target_iq_rate * 1.5:
        iq_dec = int(round(sample_rate / target_iq_rate))
        # Simple LP filter before decimation
        from scipy.signal import decimate as scipy_decimate
        x = scipy_decimate(x, iq_dec, ftype='fir', n=31)
        sample_rate = sample_rate / iq_dec

    # FM demodulation (quadrature demod)
    x = 0.5 * np.angle(x[0:-1] * np.conj(x[1:]))

    # Frequency shift to center RDS at 0 Hz
    N = len(x)
    f_o = -57e3
    t = np.arange(N) / sample_rate
    x = x * np.exp(2j * np.pi * f_o * t)

    # Low-pass filter
    taps = firwin(numtaps=101, cutoff=7.5e3, fs=sample_rate)
    x = np.convolve(x, taps, 'valid')

    # Decimate to ~25 kHz
    dec = max(1, int(round(sample_rate / 25e3)))
    x = x[::dec]
    sr = sample_rate / dec

    # Resample to 19 kHz
    up = 19
    down = int(round(sr / 1000))
    if down == 0:
        return 0
    x = resample_poly(x, up, down)
    sample_rate = 19e3

    # Mueller and Muller clock recovery
    samples_in = x
    samples_interpolated = resample_poly(samples_in, 32, 1)
    sps = 16
    mu = 0.01
    out = np.zeros(len(samples_in) + 10, dtype=np.complex64)
    out_rail = np.zeros(len(samples_in) + 10, dtype=np.complex64)
    i_in = 0
    i_out = 2
    while i_out < len(samples_in) and i_in + 32 < len(samples_in):
        out[i_out] = samples_interpolated[i_in * 32 + int(mu * 32)]
        out_rail[i_out] = int(np.real(out[i_out]) > 0) + 1j * int(np.imag(out[i_out]) > 0)
        x_val = (out_rail[i_out] - out_rail[i_out - 2]) * np.conj(out[i_out - 1])
        y_val = (out[i_out] - out[i_out - 2]) * np.conj(out_rail[i_out - 1])
        mm_val = np.real(y_val - x_val)
        mu += sps + 0.01 * mm_val
        i_in += int(np.floor(mu))
        mu = mu - np.floor(mu)
        i_out += 1
    x = out[2:i_out]

    # Costas loop
    N = len(x)
    phase = 0
    freq = 0
    alpha = 8.0
    beta = 0.02
    out = np.zeros(N, dtype=np.complex64)
    for i in range(N):
        out[i] = x[i] * np.exp(-1j * phase)
        error = np.real(out[i]) * np.imag(out[i])
        freq += (beta * error)
        phase += freq + (alpha * error)
        while phase >= 2 * np.pi:
            phase -= 2 * np.pi
        while phase < 0:
            phase += 2 * np.pi
    x = out

    # Demod BPSK
    bits = (np.real(x) > 0).astype(int)

    # Differential decoding
    bits = (bits[1:] - bits[0:-1]) % 2
    bits = bits.astype(np.uint8)

    # Block sync and group decode (from PySDR)
    syndrome = [383, 14, 303, 663, 748]
    offset_pos = [0, 1, 2, 3, 2]
    offset_word = [252, 408, 360, 436, 848]

    def calc_syndrome(x, mlen):
        reg = 0
        plen = 10
        for ii in range(mlen, 0, -1):
            reg = (reg << 1) | ((x >> (ii - 1)) & 0x01)
            if (reg & (1 << plen)):
                reg = reg ^ 0x5B9
        for ii in range(plen, 0, -1):
            reg = reg << 1
            if (reg & (1 << plen)):
                reg = reg ^ 0x5B9
        return reg & ((1 << plen) - 1)

    synced = False
    presync = False
    wrong_blocks_counter = 0
    blocks_counter = 0
    group_good_blocks_counter = 0
    reg = np.uint32(0)
    lastseen_offset_counter = 0
    lastseen_offset = 0
    groups_decoded = 0

    for i in range(len(bits)):
        reg = np.bitwise_or(np.left_shift(reg, 1), bits[i])
        if not synced:
            reg_syndrome = calc_syndrome(reg, 26)
            for j in range(5):
                if reg_syndrome == syndrome[j]:
                    if not presync:
                        lastseen_offset = j
                        lastseen_offset_counter = i
                        presync = True
                    else:
                        if offset_pos[lastseen_offset] >= offset_pos[j]:
                            block_distance = offset_pos[j] + 4 - offset_pos[lastseen_offset]
                        else:
                            block_distance = offset_pos[j] - offset_pos[lastseen_offset]
                        if (block_distance * 26) != (i - lastseen_offset_counter):
                            presync = False
                        else:
                            wrong_blocks_counter = 0
                            blocks_counter = 0
                            block_bit_counter = 0
                            block_number = (j + 1) % 4
                            group_assembly_started = False
                            synced = True
                    break
        else:
            if block_bit_counter < 25:
                block_bit_counter += 1
            else:
                good_block = False
                dataword = (reg >> 10) & 0xffff
                block_calculated_crc = calc_syndrome(dataword, 16)
                checkword = reg & 0x3ff
                if block_number == 2:
                    block_received_crc = checkword ^ offset_word[block_number]
                    if block_received_crc == block_calculated_crc:
                        good_block = True
                    else:
                        block_received_crc = checkword ^ offset_word[4]
                        if block_received_crc == block_calculated_crc:
                            good_block = True
                        else:
                            wrong_blocks_counter += 1
                else:
                    block_received_crc = checkword ^ offset_word[block_number]
                    if block_received_crc == block_calculated_crc:
                        good_block = True
                    else:
                        wrong_blocks_counter += 1

                if block_number == 0 and good_block:
                    group_assembly_started = True
                    group_good_blocks_counter = 1
                if group_assembly_started:
                    if not good_block:
                        group_assembly_started = False
                    else:
                        group_good_blocks_counter += 1
                    if group_good_blocks_counter == 5:
                        groups_decoded += 1
                        group_assembly_started = False

                block_bit_counter = 0
                block_number = (block_number + 1) % 4
                blocks_counter += 1
                if blocks_counter == 50:
                    if wrong_blocks_counter > 35:
                        synced = False
                        presync = False
                    blocks_counter = 0
                    wrong_blocks_counter = 0

    return groups_decoded


def rradio_rds_decode(meta_path):
    """Run rradio on a SigMF recording and return group count."""
    proc = subprocess.run(
        [RRADIO_BINARY, "sigmf", meta_path, "--wav", "/dev/null", "--rds-metrics"],
        capture_output=True, text=True, timeout=120,
    )
    for line in proc.stderr.splitlines():
        if line.startswith("RDSSUMMARY "):
            try:
                summary = json.loads(line[len("RDSSUMMARY "):])
                return summary.get("groups", 0)
            except json.JSONDecodeError:
                pass
    return 0


def main():
    recordings = sys.argv[1:]
    if not recordings:
        print("Usage: python py/rds_compare.py res/recordings/*.sigmf-meta")
        sys.exit(1)

    print(f"{'Recording':<20} {'PySDR':>8} {'rradio':>8} {'Winner':>10}")
    print("-" * 50)

    total_pysdr = 0
    total_rradio = 0

    for meta_path in recordings:
        name = os.path.basename(meta_path).replace(".sigmf-meta", "")

        # PySDR decode
        samples, sr, freq = load_sigmf(meta_path)
        pysdr_groups = pysdr_rds_decode(samples, sr)

        # rradio decode
        rradio_groups = rradio_rds_decode(meta_path)

        winner = "PySDR" if pysdr_groups > rradio_groups else (
            "rradio" if rradio_groups > pysdr_groups else "tie")

        print(f"{name:<20} {pysdr_groups:>8} {rradio_groups:>8} {winner:>10}")

        total_pysdr += pysdr_groups
        total_rradio += rradio_groups

    print("-" * 50)
    print(f"{'TOTAL':<20} {total_pysdr:>8} {total_rradio:>8}")


if __name__ == "__main__":
    main()
