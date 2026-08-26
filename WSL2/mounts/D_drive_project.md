# `/mnt/d/AI4S-OrgChem` ↔ `D:\AI4S-OrgChem`

WSL2 通过 DrvFs 把 Windows `D:` 挂到 `/mnt/d`。

| 视角 | 路径 |
|------|------|
| Windows | `D:\AI4S-OrgChem\` |
| WSL2 | `/mnt/d/AI4S-OrgChem/` |

**主计算工作目录就是该挂载点**（与 Windows 目录为同一树，不是第二份拷贝）。  
本目录 `WSL2/` 额外归档的是 Linux 根文件系统里、挂载点之外的运行痕迹与环境清单（`/home/cascade/...`）。
