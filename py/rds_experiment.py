#!/usr/bin/env python3
"""
RDS BLER Optimization Experiment Runner

Sweeps RDS pipeline parameters across SigMF recordings, collecting
BLER metrics to find the best configuration.

Usage:
    python py/rds_experiment.py --recordings res/recordings/*.sigmf-meta
    python py/rds_experiment.py --recordings res/recordings/*.sigmf-meta --trials 50
    python py/rds_experiment.py --recordings res/recordings/*.sigmf-meta --config base_config.json
"""

import argparse
import json
import os
import random
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path


# Parameter search space: (min, max) for each tunable parameter
SEARCH_SPACE = {
    "anti_alias_cutoff_hz": (4000.0, 9500.0),
    "anti_alias_order": (1, 6),
    "matched_filter_spans": (4, 16),
    "matched_filter_window": ["blackman", "hamming", "hann"],
    "baseband_lpf_hz": (2500.0, 6000.0),
    "agc.bandwidth_hz": (1.0, 100.0),
    "agc.target_rms": (0.0001, 0.1),
    "costas.loop_bw": (0.01, 0.2),
    "gardner.loop_bw": (0.005, 0.1),
    "sync.crc_correction_max_bits": (0, 5),
    "sync.loss_threshold": (3, 15),
}

# Default config (matches RdsConfig::default() in Rust)
DEFAULT_CONFIG = {
    "baseband_lpf_hz": 4000.0,
    "anti_alias_cutoff_hz": 9000.0,
    "anti_alias_order": 3,
    "matched_filter_spans": 8,
    "matched_filter_window": "blackman",
    "agc": {"bandwidth_hz": 10.0, "target_rms": 0.001},
    "costas": {"loop_bw": 0.05},
    "gardner": {"loop_bw": 0.01},
    "sync": {"crc_correction_max_bits": 2, "loss_threshold": 5},
}

BINARY = "target/release/rradio"


@dataclass
class TrialResult:
    config: dict
    recording: str
    groups: int = 0
    duration: float = 0.0
    mean_bler: float = 1.0
    metrics: list = field(default_factory=list)


def set_nested(d: dict, key: str, value) -> dict:
    """Set a nested key like 'agc.bandwidth_hz' in a dict."""
    d = json.loads(json.dumps(d))  # deep copy
    parts = key.split(".")
    target = d
    for part in parts[:-1]:
        target = target.setdefault(part, {})
    target[parts[-1]] = value
    return d


def sample_random_config() -> dict:
    """Generate a random config by sampling each parameter."""
    config = json.loads(json.dumps(DEFAULT_CONFIG))
    for key, space in SEARCH_SPACE.items():
        if isinstance(space, list):
            value = random.choice(space)
        elif isinstance(space[0], int):
            value = random.randint(space[0], space[1])
        else:
            # Log-uniform for parameters spanning orders of magnitude
            lo, hi = space
            if hi / lo > 10:
                import math
                value = math.exp(random.uniform(math.log(lo), math.log(hi)))
            else:
                value = random.uniform(lo, hi)
        config = set_nested(config, key, value)
    return config


def run_trial(recording: str, config: dict, wav_output: str = "/dev/null") -> TrialResult:
    """Run the binary with a config and collect BLER metrics."""
    result = TrialResult(config=config, recording=recording)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(config, f)
        config_path = f.name

    try:
        cmd = [
            BINARY, "sigmf", recording,
            "--wav", wav_output,
            "--rds-metrics",
            "--rds-config", config_path,
        ]

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
        )

        # Parse RDSMETRIC lines from stderr
        for line in proc.stderr.splitlines():
            if line.startswith("RDSMETRIC "):
                try:
                    metric = json.loads(line[len("RDSMETRIC "):])
                    result.metrics.append(metric)
                except json.JSONDecodeError:
                    pass
            elif line.startswith("RDSSUMMARY "):
                try:
                    summary = json.loads(line[len("RDSSUMMARY "):])
                    result.groups = summary.get("groups", 0)
                    result.duration = summary.get("duration", 0.0)
                except json.JSONDecodeError:
                    pass

        # Compute mean BLER from metrics (skip first 20% for settling)
        if result.metrics:
            settle = max(1, len(result.metrics) // 5)
            settled = result.metrics[settle:]
            if settled:
                result.mean_bler = sum(m["bler"] for m in settled) / len(settled)
            else:
                result.mean_bler = 1.0
        else:
            result.mean_bler = 1.0

    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT on {recording}", file=sys.stderr)
        result.mean_bler = 1.0
    except Exception as e:
        print(f"  ERROR: {e}", file=sys.stderr)
        result.mean_bler = 1.0
    finally:
        os.unlink(config_path)

    return result


def evaluate_config(recordings: list, config: dict) -> float:
    """Run a config across all recordings and return the aggregate score.
    Lower is better (mean of mean BLERs across recordings)."""
    total_bler = 0.0
    for rec in recordings:
        result = run_trial(rec, config)
        total_bler += result.mean_bler
    return total_bler / len(recordings)


def evaluate_trial(args_tuple):
    """Wrapper for parallel execution. Takes (trial_num, recordings, config)."""
    trial_num, recordings, config = args_tuple
    score = evaluate_config(recordings, config)
    return trial_num, score, config


def main():
    parser = argparse.ArgumentParser(
        description="RDS BLER optimization experiment runner")
    parser.add_argument("--recordings", nargs="+", required=True,
                        help="SigMF recording .sigmf-meta files")
    parser.add_argument("--trials", type=int, default=100,
                        help="Number of random configs to try (default: 100)")
    parser.add_argument("--config", type=str, default=None,
                        help="Base config JSON to evaluate (skips random search)")
    parser.add_argument("--build", action="store_true",
                        help="Build the release binary before running")
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel workers (default: CPU count)")
    args = parser.parse_args()

    # Verify recordings exist
    for rec in args.recordings:
        if not os.path.exists(rec):
            print(f"Recording not found: {rec}", file=sys.stderr)
            sys.exit(1)

    # Optionally build
    if args.build:
        print("Building release binary...", file=sys.stderr)
        result = subprocess.run(
            ["cargo", "build", "--release"],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"Build failed:\n{result.stderr}", file=sys.stderr)
            sys.exit(1)

    # Verify binary exists
    if not os.path.exists(BINARY):
        print(f"Binary not found: {BINARY}", file=sys.stderr)
        print("Run 'cargo build --release' first, or use --build", file=sys.stderr)
        sys.exit(1)

    # Single config evaluation
    if args.config:
        with open(args.config) as f:
            config = json.load(f)
        # Merge with defaults for missing fields
        merged = json.loads(json.dumps(DEFAULT_CONFIG))
        merged.update(config)

        print(f"Evaluating config: {args.config}")
        for rec in args.recordings:
            result = run_trial(rec, merged)
            print(f"  {Path(rec).stem}: {result.groups} groups, "
                  f"mean BLER={result.mean_bler:.4f}")
        sys.exit(0)

    # Random search
    workers = args.workers or os.cpu_count() or 4
    print(f"Running {args.trials} random trials across "
          f"{len(args.recordings)} recordings with {workers} workers...")
    print(f"Recordings: {[Path(r).stem for r in args.recordings]}")
    print()

    best_score = 1.0
    best_config = None
    results_log = []

    # Always evaluate default first (sequential)
    print(f"Trial 0/{args.trials}: DEFAULT CONFIG")
    default_score = evaluate_config(args.recordings, DEFAULT_CONFIG)
    print(f"  Score (mean BLER): {default_score:.4f}")
    results_log.append({"trial": 0, "score": default_score,
                        "config": DEFAULT_CONFIG, "type": "default"})
    best_score = default_score
    best_config = DEFAULT_CONFIG
    print()

    # Generate all random configs upfront
    trial_args = []
    for trial in range(1, args.trials + 1):
        config = sample_random_config()
        trial_args.append((trial, args.recordings, config))

    # Run in parallel
    completed = 0
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(evaluate_trial, ta): ta[0]
                   for ta in trial_args}

        for future in as_completed(futures):
            trial_num, score, config = future.result()
            completed += 1
            improved = score < best_score

            if improved:
                best_score = score
                best_config = config
                print(f"  [{completed}/{args.trials}] Trial {trial_num}: "
                      f"Score={score:.4f} *** NEW BEST ***")
            elif completed % 10 == 0 or completed == args.trials:
                print(f"  [{completed}/{args.trials}] "
                      f"best so far: {best_score:.4f}")

            results_log.append({"trial": trial_num, "score": score,
                                "config": config, "type": "random"})

    # Report results
    print()
    print("=" * 60)
    print(f"BEST SCORE: {best_score:.4f} (mean BLER across recordings)")
    print("BEST CONFIG:")
    print(json.dumps(best_config, indent=2))
    print()

    # Save results
    output_path = "rds_experiment_results.json"
    with open(output_path, "w") as f:
        json.dump({
            "best_score": best_score,
            "best_config": best_config,
            "trials": results_log,
            "recordings": args.recordings,
        }, f, indent=2)
    print(f"Full results saved to {output_path}")


if __name__ == "__main__":
    main()
