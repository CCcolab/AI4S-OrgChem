"""Wait for the stage-A pre-relaxed seed, then run the theta scan (stage B/C).

Keeps the two expensive stages back to back without idle wall time, while
still letting stage A be inspected/reused on its own.
"""
from __future__ import annotations

import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SEED = ROOT / "results" / "P3" / "seed_prerelaxed.xyz"
LOG = ROOT / "results" / "P3" / "logs" / "stageB.log"
WAIT_LIMIT_S = 5 * 3600


def stage_a_alive() -> bool:
    out = subprocess.run(
        ["pgrep", "-af", "run_tight"], capture_output=True, text=True, check=False
    ).stdout
    return any("seed-only" in ln for ln in out.splitlines())


def main() -> int:
    t0 = time.time()
    print(f"[chain] waiting for {SEED}", flush=True)
    while not (SEED.exists() and SEED.stat().st_size > 0):
        if not stage_a_alive():
            print("[chain] stage A gone and no seed: aborting", flush=True)
            return 1
        if time.time() - t0 > WAIT_LIMIT_S:
            print("[chain] timed out waiting for seed", flush=True)
            return 1
        time.sleep(30)

    print(f"[chain] seed ready {datetime.now():%H:%M:%S}", flush=True)
    out = ROOT / "results" / "P3"

    # The B3LYP free optimum is the decisive datum for P3: it locates the true
    # minimum with no constraint, at the same level used for the energies.
    jobs = [
        (
            "free_b3lyp",
            [
                "--seed-xyz", str(SEED), "--refine-seed",
                "--cascade", "6-31g*:60",
                "--geom-method", "B3LYP", "--seed-only",
                "--out", str(out / "free_b3lyp"),
            ],
        ),
        (
            "stageB",
            [
                "--seed-xyz", str(SEED),
                "--geom-method", "RHF", "--geom-basis", "3-21g",
                "--energy-method", "B3LYP", "--energy-basis", "6-31g*",
                "--angles", "0,15,30,45,60,75,90",
                "--scan-maxiter", "40", "--k-penalty", "2.0",
                "--direction", "both",
                "--out", str(out),
            ],
        ),
    ]
    LOG.parent.mkdir(parents=True, exist_ok=True)
    procs = []
    for name, args in jobs:
        fh = (LOG.parent / f"{name}.log").open("w", encoding="utf-8")
        procs.append(
            (
                name,
                fh,
                subprocess.Popen(
                    [sys.executable, "-m", "src.p3_nba.run_tight", *args],
                    cwd=ROOT,
                    stdout=fh,
                    stderr=subprocess.STDOUT,
                ),
            )
        )
        print(f"[chain] started {name} at {datetime.now():%H:%M:%S}", flush=True)
    rc = 0
    for name, fh, p in procs:
        code = p.wait()
        fh.close()
        print(f"[chain] {name} exit={code} at {datetime.now():%H:%M:%S}", flush=True)
        rc = rc or code
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
