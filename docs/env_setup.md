# 环境记录（只读核查，未改配置）

> 按项目要求：**不安装、不升级、不修改** WSL 中任何已有软件或 GPU 相关设置。本文件只记录核查结果。  
> **完整 WSL2 计算平面归档**见根目录 [`WSL2/`](../WSL2/)（含 `inventory/nvidia-smi.txt`、`python_system_pyscf.txt`、home 侧小目录拷贝）。

## 核查时间
- 初查：2026-08-23  
- WSL2 目录快照：2026-08-26（`WSL2/MANIFEST.txt`）

## 结果
| 项 | 值 |
|----|-----|
| 发行版 | Ubuntu 24.04.4 LTS（WSL2） |
| 工作目录 | `/mnt/d/AI4S-OrgChem` ≡ `D:\AI4S-OrgChem`（同一树） |
| Python | `/usr/bin/python3` · 3.12.3 |
| PySCF | **2.14.0**（`~/.local/lib/python3.12/site-packages/pyscf`） |
| GPU | NVIDIA GeForce RTX 4060 Laptop（`nvidia-smi` 可见；快照见 `WSL2/inventory/nvidia-smi.txt`） |

## 使用约定
- **所有量子化学计算在 WSL2 中**调用现有 `python3` + `pyscf`（含 GPU 可见性）
- Windows 侧：文档、交付、Canvas；不在 Windows conda 中另装 PySCF 作主证据
- 探测 / 快照：`tools/check_wsl_env.sh`、`tools/wsl_snapshot.sh`（只读 / 归档，不改环境）
