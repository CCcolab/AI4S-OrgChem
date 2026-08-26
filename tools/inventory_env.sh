#!/usr/bin/env bash
# Read-only environment inventory. Does NOT install or modify anything.
echo "=== WSL OS ==="
uname -a
. /etc/os-release 2>/dev/null && echo "PRETTY_NAME=$PRETTY_NAME"
echo "=== WSL CPU ==="
nproc
lscpu 2>/dev/null | egrep 'Model name|CPU\(s\)|Thread|Core|Architecture|MHz' | head -20
echo "=== WSL MEM ==="
free -h
echo "=== WSL PYTHON ==="
command -v python3
python3 --version
echo "=== WSL PYSCF ==="
python3 - <<'PY'
try:
    import pyscf
    print("pyscf", pyscf.__version__)
except Exception as e:
    print("pyscf_error", type(e).__name__, e)
try:
    import numpy as np
    print("numpy", np.__version__)
except Exception as e:
    print("numpy_error", type(e).__name__, e)
try:
    import scipy
    print("scipy", scipy.__version__)
except Exception as e:
    print("scipy_error", type(e).__name__, e)
try:
    import h5py
    print("h5py", h5py.__version__)
except Exception as e:
    print("h5py_error", type(e).__name__, e)
PY
echo "=== WSL GPU ==="
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
else
  echo "nvidia-smi not found"
fi
echo "=== WSL CUDA ==="
command -v nvcc >/dev/null && nvcc --version | tail -1 || echo "nvcc not found"
echo "CUDA_HOME=${CUDA_HOME:-unset}"
echo "=== WSL DISK /mnt/d ==="
df -h /mnt/d 2>/dev/null | tail -1
echo "=== WSL EXTRA QC ==="
for pkg in psi4 rdkit ase matplotlib pandas; do
  python3 - <<PY
try:
    import $pkg as m
    v = getattr(m, "__version__", "ok")
    print("$pkg", v)
except Exception as e:
    print("$pkg", "not_imported", type(e).__name__)
PY
done
