#!/bin/bash
# Snapshot WSL2 project-related layout into Windows tree: /mnt/d/AI4S-OrgChem/WSL2/
set -euo pipefail
ROOT=/mnt/d/AI4S-OrgChem/WSL2
INV="$ROOT/inventory"
HOME_SNAP="$ROOT/home_cascade"
PUB_USER=wsluser
PUB_HOME=/home/wsluser
PUB_HOST=wsl-host
mkdir -p "$INV" "$HOME_SNAP" "$ROOT/mounts"

redact_inventory() {
  local user host
  user="$(whoami)"
  host="$(uname -n)"
  for f in "$INV"/*.txt; do
    [ -f "$f" ] || continue
    sed -i \
      -e "s|/home/${user}|${PUB_HOME}|g" \
      -e "s|${user} ${user}|${PUB_USER} ${PUB_USER}|g" \
      -e "s/^Linux ${host} /Linux ${PUB_HOST} /" \
      "$f"
  done
}

# --- environment inventory (lightweight; no multi-GB envs) ---
uname -a > "$INV/uname.txt"
cp -f /etc/os-release "$INV/os-release.txt" 2>/dev/null || true
nvidia-smi > "$INV/nvidia-smi.txt" 2>&1 || echo "nvidia-smi unavailable" > "$INV/nvidia-smi.txt"
{
  echo "whoami: ${PUB_USER}"
  echo "HOME: ${PUB_HOME}"
  echo "pwd: $(pwd)"
  echo "date: $(date -Iseconds)"
  echo "WSL_DISTRO: ${WSL_DISTRO_NAME:-unknown}"
} > "$INV/session.txt"

/usr/bin/python3 - <<'PY' > "$INV/python_system_pyscf.txt" 2>&1
import sys, platform
print("executable:", sys.executable)
print("version:", sys.version)
print("platform:", platform.platform())
try:
    import pyscf
    print("pyscf_version:", pyscf.__version__)
    print("pyscf_file:", pyscf.__file__)
except Exception as e:
    print("pyscf_error:", e)
PY

# user-site freeze (where pyscf 2.14.0 lives)
/usr/bin/python3 -m pip freeze --user > "$INV/pip_freeze_user.txt" 2>/dev/null || true
/usr/bin/python3 -m pip show pyscf numpy scipy > "$INV/pip_show_pyscf_stack.txt" 2>&1 || true

# micromamba env metadata only (NOT the 8.9G tree)
if [ -x "$HOME/micromamba/envs/ai4orgchem/bin/python" ]; then
  "$HOME/micromamba/envs/ai4orgchem/bin/python" -V > "$INV/mamba_ai4orgchem_python_version.txt" 2>&1 || true
  if command -v micromamba >/dev/null 2>&1; then
    micromamba list -n ai4orgchem > "$INV/mamba_ai4orgchem_packages.txt" 2>&1 || true
  elif [ -x "$HOME/micromamba/bin/micromamba" ]; then
    "$HOME/micromamba/bin/micromamba" list -n ai4orgchem > "$INV/mamba_ai4orgchem_packages.txt" 2>&1 || true
  fi
  # tiny pointer file instead of copying env
  {
    echo "path: ${PUB_HOME}/micromamba/envs/ai4orgchem"
    du -sh "$HOME/micromamba/envs/ai4orgchem" 2>/dev/null
    echo "NOTE: full env NOT copied into Windows tree (multi-GB); see packages list above."
  } > "$INV/mamba_ai4orgchem_POINTER.txt"
fi

# directory map of home project-related paths
{
  echo "# WSL Linux filesystem paths related to AI4S-OrgChem"
  echo
  du -sh "$HOME/p3run" "$HOME/ai4orgchem_env" "$HOME/ai4orgchem-scratch" \
    "$HOME/micromamba/envs/ai4orgchem" "$HOME/micromamba/envs/ai4orgchem-nequip" \
    "$HOME/.local/lib/python3.12/site-packages/pyscf" 2>/dev/null
  echo
  echo "# Windows mount (same inode tree as D:\\AI4S-OrgChem)"
  ls -la /mnt/d/AI4S-OrgChem | head -30
} > "$INV/path_map.txt"

redact_inventory

# --- copy small project-related home dirs ---
rsync -a --delete \
  --exclude '__pycache__/' \
  --exclude '*.pyc' \
  "$HOME/p3run/" "$HOME_SNAP/p3run/"

rsync -a --delete \
  --exclude '__pycache__/' \
  "$HOME/ai4orgchem_env/" "$HOME_SNAP/ai4orgchem_env/" 2>/dev/null || \
  cp -a "$HOME/ai4orgchem_env" "$HOME_SNAP/" 2>/dev/null || true

rsync -a --delete \
  "$HOME/ai4orgchem-scratch/" "$HOME_SNAP/ai4orgchem-scratch/" 2>/dev/null || \
  cp -a "$HOME/ai4orgchem-scratch" "$HOME_SNAP/" 2>/dev/null || true

# Replace Linux→/usr symlinks (broken on Windows/GitHub) with notes
if [ -d "$HOME_SNAP/ai4orgchem_env" ]; then
  rm -f "$HOME_SNAP/ai4orgchem_env/bin/python" \
        "$HOME_SNAP/ai4orgchem_env/bin/python3" \
        "$HOME_SNAP/ai4orgchem_env/bin/python3.12" \
        "$HOME_SNAP/ai4orgchem_env/lib64" 2>/dev/null || true
  mkdir -p "$HOME_SNAP/ai4orgchem_env/bin"
  printf '%s\n' 'Symlinks to /usr/bin/python3 removed for Windows/GitHub.' \
    'Use WSL: /usr/bin/python3' > "$HOME_SNAP/ai4orgchem_env/bin/README.txt"
  printf '%s\n' 'Empty venv shell; PySCF is user-site (~/.local). See WSL2/inventory/.' \
    > "$HOME_SNAP/ai4orgchem_env/SYMLINKS_NOTE.txt"
  if [ -f "$HOME_SNAP/ai4orgchem_env/pyvenv.cfg" ]; then
    sed -i -e "s|${HOME}|${PUB_HOME}|g" "$HOME_SNAP/ai4orgchem_env/pyvenv.cfg"
  fi
fi

# mounts note
cat > "$ROOT/mounts/D_drive_project.md" <<EOF
# \`/mnt/d/AI4S-OrgChem\` ↔ \`D:\\AI4S-OrgChem\`

WSL2 通过 DrvFs 把 Windows \`D:\` 挂到 \`/mnt/d\`。

| 视角 | 路径 |
|------|------|
| Windows | \`D:\\AI4S-OrgChem\\\` |
| WSL2 | \`/mnt/d/AI4S-OrgChem/\` |

**主计算工作目录就是该挂载点**（与 Windows 目录为同一树，不是第二份拷贝）。  
本目录 \`WSL2/\` 额外归档的是 Linux 根文件系统里、挂载点之外的运行痕迹与环境清单（\`${PUB_HOME}/...\`）。
EOF

# manifest
{
  echo "created: $(date -Iseconds)"
  echo "host: ${PUB_HOST}"
  find "$ROOT" -type f | sed "s|^$ROOT/||" | sort
} > "$ROOT/MANIFEST.txt"

echo "OK -> $ROOT"
