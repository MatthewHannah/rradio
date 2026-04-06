# Testing RDS against SigMF recordings

## Quick reference

Run all recordings and check group counts:

```bash
for rec in res/recordings/*.sigmf-meta; do
  name=$(basename "$rec" .sigmf-meta)
  groups=$(cargo run --release -- sigmf "$rec" --wav /dev/null --rds-metrics 2>&1 \
    | grep "RDSSUMMARY" | sed 's/.*"groups":\([0-9]*\).*/\1/')
  echo "$name: ${groups:-0}"
done
```

## Expected baseline (v5 pipeline)

| Recording | Groups | Notes |
|-----------|--------|-------|
| 94_9_30s | 327 | Strong station, gold standard |
| 99_7_30s | 329 | Strong station |
| 100_5_30s | 328 | Strong station |
| 96_1_30s | 330 | Strong station |
| 97_1_30s | 323 | Strong station |
| 92_9_1min | 565 | Weak station (-8 dB SNR), 1 minute |
| 99_7_5min_live | 3414 | 5-minute live capture with dropout at ~148s |

Strong stations should all be within ±5 groups of these values. 92.9 is the
most sensitive to changes — it's the weakest signal and drops first when
something regresses.

## Redsea reference pipeline

Compare against the redsea-based pipeline with `USE_REDSEA=1`:

```bash
USE_REDSEA=1 cargo run --release -- sigmf res/recordings/94_9_30s.sigmf-meta \
  --wav /dev/null --rds-metrics 2>&1 | grep "RDSSUMMARY"
```

Both pipelines should produce nearly identical results on recordings (within
±5 groups). Significant divergence means the shared signal_pipeline() front-end
changed in a way that affects one pipeline more than the other.

## What to watch for

**Strong stations (94.9, 99.7, 100.5, 96.1, 97.1):**
These are near the theoretical maximum group rate (~330 groups/30s = 11 groups/s).
A drop of more than 5 groups indicates a real regression.

**Weak station (92.9):**
This station has ~25% BLER. It's the canary — most regressions show up here
first. At the start of this project it decoded 0 groups; it should now decode
550+ in the 1-minute recording.

**5-minute live (99.7):**
Tests dropout recovery. There's an audio transient at t≈148s that causes a
biphase polarity flip. The block sync should re-lock within <0.5s (previously
took 24s before the CRC loss detection fix). Total groups should be ~3414.

## CLI flags

| Flag | Purpose |
|------|---------|
| `--rds-metrics` | Print RDSSUMMARY JSON at end and per-group RDSMETRIC lines |
| `--rds-debug` | Verbose sync state logging (FOUND/FAIL/LOST/LOCKED events) |
| `--wav /dev/null` | Suppress audio output (required even for metrics-only runs) |
| `--diag <path.csv>` | Per-chip diagnostic CSV (Costas phase, timing period, AGC, etc.) |

## Performance

Processing time for 94_9_30s should be under 5s on a modern machine (6× real-time
or better). The 5-minute recording should complete in ~30s. If processing time
doubles, check whether a change added per-sample overhead to the signal_pipeline()
hot path (1.2 MHz sample rate).
