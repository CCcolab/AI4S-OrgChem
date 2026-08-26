#!/usr/bin/env bash
# Read-only environment probe. Does NOT install or modify anything.
echo "=== python ==="
command -v python3 || true
python3 --version || true
echo "=== pyscf ==="
python3 - <<'PY'
try:
    import pyscf
    print("pyscf", pyscf.__version__)
except Exception as e:
    print("pyscf_import_error:", type(e).__name__, e)
PY
echo "=== gpu ==="
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi -L
else
  echo "nvidia-smi not found"
fi
