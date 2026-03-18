#!/usr/bin/env python3
"""
Focused grid search over 3 RDS parameters:
  - costas.loop_bw
  - agc.bandwidth_hz
  - sync.loss_threshold

All other parameters are fixed at their determined-best values.
"""

import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

BINARY = "target/release/rradio"

# Fixed parameters from experiment analysis
BASE_CONFIG = {
    "baseband_lpf_hz": 4000.0,
    "anti_alias_cutoff_hz": 6266.0,
    "anti_alias_order": 1,
    "matched_filter_spans": 5,
    "matched_filter_window": "hann",
    "agc": {"bandwidth_hz": 10.0, "target_rms": 0.001},
    "costas": {"loop_bw": 0.05},
    "gardner": {"loop_bw": 0.01},
    "sync": {"crc_correction_max_bits": 2, "loss_threshold": 5},
}

# Grid values for the 3 tunable parameters
COSTAS_BW = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.18]
AGC_BW = [2, 5, 10, 20, 30, 50, 80]
SYNC_LOSS = [3, 5, 8, 12]


def make_config(costas_bw, agc_bw, loss_thresh):
    config = json.loads(json.dumps(BASE_CONFIG))
    config["costas"]["loop_bw"] = costas_bw
    config["agc"]["bandwidth_hz"] = agc_bw
    config["sync"]["loss_threshold"] = loss_thresh
    return config


def run_trial(recording, config):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(config, f)
        config_path = f.name
    try:
        proc = subprocess.run(
            [BINARY, "sigmf", recording, "--wav", "/dev/null",
             "--rds-metrics", "--rds-config", config_path],
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
    except Exception:
        return 0
    finally:
        os.unlink(config_path)


def evaluate_grid_point(args):
    costas_bw, agc_bw, loss_thresh, recordings = args
    config = make_config(costas_bw, agc_bw, loss_thresh)
    total_groups = sum(run_trial(rec, config) for rec in recordings)
    return costas_bw, agc_bw, loss_thresh, total_groups


def main():
    recordings = sys.argv[1:]
    if not recordings:
        print("Usage: python py/rds_grid_search.py rec1.sigmf-meta rec2.sigmf-meta ...")
        sys.exit(1)

    grid = [(c, a, s, recordings)
            for c in COSTAS_BW for a in AGC_BW for s in SYNC_LOSS]

    total = len(grid)
    workers = os.cpu_count() or 4
    print(f"Grid search: {len(COSTAS_BW)}×{len(AGC_BW)}×{len(SYNC_LOSS)} = "
          f"{total} configs × {len(recordings)} recordings, {workers} workers")

    results = []
    best_groups = 0
    completed = 0

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(evaluate_grid_point, g): g for g in grid}
        for future in as_completed(futures):
            costas_bw, agc_bw, loss_thresh, total_groups = future.result()
            completed += 1
            results.append({
                "costas_bw": costas_bw, "agc_bw": agc_bw,
                "loss_threshold": loss_thresh, "groups": total_groups
            })
            if total_groups > best_groups:
                best_groups = total_groups
                print(f"  [{completed}/{total}] costas={costas_bw} agc={agc_bw} "
                      f"loss={loss_thresh}: {total_groups} groups *** BEST ***")
            elif completed % 20 == 0:
                print(f"  [{completed}/{total}] best so far: {best_groups} groups")

    results.sort(key=lambda r: -r["groups"])
    print()
    print("=" * 60)
    print("TOP 10 CONFIGS:")
    for i, r in enumerate(results[:10]):
        print(f"  {i+1}. costas={r['costas_bw']:.2f} agc={r['agc_bw']:.0f} "
              f"loss={r['loss_threshold']} → {r['groups']} groups")

    print()
    best = results[0]
    print(f"BEST: costas={best['costas_bw']:.2f} agc={best['agc_bw']:.0f} "
          f"loss={best['loss_threshold']} → {best['groups']} groups")

    # Save full config
    best_config = make_config(best["costas_bw"], best["agc_bw"], best["loss_threshold"])
    with open("res/rds_grid_best.json", "w") as f:
        json.dump(best_config, f, indent=2)
    print("Saved to res/rds_grid_best.json")

    # Save all results
    with open("rds_grid_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"Full results ({len(results)} configs) saved to rds_grid_results.json")


if __name__ == "__main__":
    main()
