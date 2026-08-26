#!/usr/bin/env bash
# Read-only env assumed. Run P1 from repo root on /mnt/d.
set -euo pipefail
cd /mnt/d/AI4S-OrgChem
export PYTHONPATH="/mnt/d/AI4S-OrgChem:${PYTHONPATH:-}"
METHOD="${1:-B3LYP}"
BASIS="${2:-6-31g*}"
python3 -m src.p1_butadiene.run --method "$METHOD" --basis "$BASIS"
