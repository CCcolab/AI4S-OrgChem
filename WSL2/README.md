# WSL2/ — 计算平面归档（Ubuntu 24.04 + GPU）

本项目的**主要量子化学计算全部在 WSL2 中执行**（Ubuntu 24.04.4 LTS，NVIDIA GPU 对本发行版可见），**不是**在原生 Windows Python / Anaconda 中跑 PySCF。

本目录把「Windows 工程树」与「WSL Linux 本机文件系统中的运行痕迹」显式对齐，并说明**在 WSL2 里算了什么、怎么算、结果写到哪里**。

---

## 1. 计算总览（先读）

| 项 | 内容 |
|----|------|
| **计算平面** | WSL2 · `Ubuntu-24.04` · 内核 `5.15.x-microsoft-standard-WSL2` |
| **工作目录** | `/mnt/d/AI4S-OrgChem` ≡ `D:\AI4S-OrgChem`（同一树，经 DrvFs 挂载） |
| **解释器** | `/usr/bin/python3` · **3.12.3** |
| **主证据引擎** | **PySCF 2.14.0**（user site：`~/.local/lib/python3.12/site-packages/pyscf`） |
| **GPU** | NVIDIA GeForce **RTX 4060 Laptop**（约 8 GB）；`nvidia-smi` 在 WSL2 内可见（驱动/CUDA 快照见 `inventory/nvidia-smi.txt`） |
| **主方法层次** | 多数命题：**B3LYP / 6-31G\***；敏感性含 RHF、MP2、ZPE、双基组、2011-lite |
| **特殊层次** | P2 LFMO-lite：**RHF / STO-3G**（对齐原著早期 π–σ 分解量级） |
| **环境策略** | **不擅自** `apt`/`pip` 安装或升级；缺包只记录缺口（见项目边界 Rule） |
| **产出落盘** | 计算写入 `results/Pn/` → 精选进 `deliverables/`；正式判定只在 `VERDICT.md` |
| **质控** | 每条正式结果须过 G1–G5（拓扑 / 几何 / 收敛 / 尺度 / 路径）；见 [`docs/quality_gates.md`](../docs/quality_gates.md) |

一句话：AI4S Agent 在 Cursor 侧编排，**真正的 SCF / 优化 / 扫描 / 定域能量分解都在 WSL2 + PySCF 上跑完**，再把带 `quality_gate` 的表写回同一工程树。

---

## 2. 双平面关系

| 平面 | 路径 | 角色 |
|------|------|------|
| **Windows 工程根** | `D:\AI4S-OrgChem\` | 文档、交付、自研 `src/`、`results/` 主树、本 `WSL2/` 归档 |
| **WSL 挂载同一树** | `/mnt/d/AI4S-OrgChem/` | **主计算工作目录**（与上表为同一目录，非第二份拷贝） |
| **WSL Linux 本机盘** | `/home/wsluser/…` | PySCF 安装位置、micromamba、早期 `p3run` 切片等 |

```
Windows D:\AI4S-OrgChem\          <──DrvFs──>   WSL /mnt/d/AI4S-OrgChem/
        │                                         │
        │  src/  results/  deliverables/          │  PYTHONPATH=. python3 -m src.…
        │                                         │  → 写回 results/Pn/tables/
        └── WSL2\  (本目录)  <──归档──  /home/wsluser/{p3run,ai4orgchem_*,…}
                                              + inventory/（环境指纹，非 8.9G 全量 env）
```

详见 [`mounts/D_drive_project.md`](mounts/D_drive_project.md)。

---

## 3. 详细计算介绍

### 3.1 为什么必须在 WSL2 算

1. **引擎绑定**：项目边界规定主证据引擎为 WSL 中既有 **PySCF**，不在 Windows 侧另装一套作判定依据。
2. **Linux 工具链**：几何优化（SciPy BFGS + 解析梯度）、定域态 DIIS→Newton、大批量扫描，均在 Linux Python 栈上完成。
3. **GPU 可见性**：WSL2 暴露 NVIDIA 设备（本机 RTX 4060 Laptop）；环境核查与快照均要求 `nvidia-smi` 可读，保证计算节点具备加速硬件，而不是“纯文档机”。
4. **路径统一**：脚本固定以 `/mnt/d/AI4S-OrgChem` 为根，`results/` 与 `deliverables/` 对 Windows / GitHub 读者直接可见，无需二次搬运数值。

### 3.2 硬件与软件栈（快照口径）

| 层 | 快照内容（详见 `inventory/`） |
|----|------------------------------|
| OS | Ubuntu 24.04.4 LTS（`os-release.txt`） |
| 内核 | `Linux …-microsoft-standard-WSL2`（`uname.txt`） |
| GPU | GeForce RTX 4060 Laptop · Driver 566.24 · CUDA 12.7 报告（`nvidia-smi.txt`） |
| Python | `/usr/bin/python3` 3.12.3（`python_system_pyscf.txt`） |
| PySCF | **2.14.0**，安装于 `~/.local/.../site-packages/pyscf`（约 190 MB） |
| 依赖 | `pip_freeze_user.txt` / `pip_show_pyscf_stack.txt`（numpy、scipy 等） |
| 旁路 env | `micromamba/envs/ai4orgchem` ≈ 8.9 GB（**未全量拷贝**；仅包清单 + 指针） |

> 说明：电子结构主路径是 **PySCF SCF/DFT**；GPU 作为 WSL2 计算节点的可用加速资源记录在案。主判定能量不依赖 MACE 等 ML 势，也不依赖未批准的旁路引擎。

### 3.3 在 WSL2 里做了哪些计算（任务类型）

| 类型 | 典型操作 | 出现在 |
|------|----------|--------|
| **单点能量** | RHF / B3LYP / MP2 能量；Ee/EN 分解 | P1、P3 单点、P4 扫描点 |
| **几何优化** | BFGS + `mf.nuc_grad_method()`；约束二面角 / 对称参数优化 | P1 物种优化、P3 θ 约束、P7 D₃ₕ、P8 G* |
| **势能面 / 扫描** | 扭转角 θ、BLA δ、Δr、环尺寸 N | P3、P4、P7、P9 |
| **定域参考态** | 2007-GL：π-AO 识别 + Fock/重叠块删除；GE-m；ESE；杂环 hetero_gl；PLG | P1 L3、P5、P6、P7、P8、P9 |
| **敏感性** | 层次（RHF/MP2/ZPE）、基组（6-31G* vs 6-31G）、2011-lite（`zero_exchange`）、BLA 几何 | P1、P5–P9 |
| **LFMO-lite** | 按 Fig. 5-15 公开态定义的 AO 代理密度能量（非 Kost 复现） | P2 |
| **质控断言** | 拓扑 / 最小间距 / 收敛或滞后代理 / 能量窗 / 路径洁净 | 全命题；失败则 `invalid_*/` |

单位约定：内部 Hartree；报告 **kcal/mol**（1 Ha = 627.5095）；稳定化为负、去稳定为正。

### 3.4 九命题 × WSL2 定判入口

下列命令均在 **WSL2** 内、仓库根执行。定判以表中入口为准；更早的 `run_v1` / `run_pilot` 仅为版本链。

```bash
cd /mnt/d/AI4S-OrgChem
export PYTHONPATH=.
```

| 命题 | 在 WSL2 上算什么 | 定判入口模块 | 主方法 | 主结果（`results/`） |
|------|------------------|--------------|--------|----------------------|
| **P1** | 氢化热 CE₁/CE₂；RHF/B3LYP/MP2/ZPE；自研 2007-GL | `src.p1_butadiene.run` · `run_zpe` · `run_gl` | B3LYP/6-31G* 等 | `ce_*.json` · `gl2007_*.json` |
| **P2** | 汇总 + NBA 扭转 LFMO-lite（EV/Enσσ/Eπσ） | `src.p2_aggregate.run_lfmo_lite` · `summarize` | RHF/STO-3G（lite） | `p2_v2_*.json` |
| **P3** | NBA 双向 θ 约束 PES + 自由 B3LYP 复核 | `src.p3_nba.run_tight` | RHF 几何 + B3LYP 单点 | `pes_tight_*.csv/json` |
| **P4** | 苯 BLA 扫描；Ee/EN/E_tot | `src.p4_benzene.run` | B3LYP/6-31G* | `bla_scan_*.json` |
| **P5** | ΔEAm、Δr；v5c 三异议 | `src.p5_local.run_v5c_objections` | 2007 GL + B3LYP | `p5_v5c_objections_*.json` |
| **P6** | 苯 ESE、环丁二烯 ΔEA；v4 三异议 | `src.p6_ese.run_v4_objections` | 2007 + B3LYP | `p6_v4_objections_*.json` |
| **P7** | C₁₂H₆：G vs PLG 的 Δr | `src.p7_strained.run_v1c` | B3LYP；PLG DIIS→Newton | `p7_v1c_*.json` |
| **P8** | 呋喃/吡咯/噁唑 LDE vs 苯 ESE | `src.p8_furan.run_v2` | B3LYP；双基组 | `p8_v2_*.json` |
| **P9** | [N]annulene VDE/π 与 gap 收敛 | `src.p9_annulene.run_v2b` | B3LYP；2011-lite/BLA | `p9_v2b.json` |

全套命令亦见 [`deliverables/final/FULL_REPORT.md`](../deliverables/final/FULL_REPORT.md) §7。

### 3.5 单次作业的数据流（WSL2 内）

```
implementation/Pn_*（预注册判据）
        ↓
src/p*_…/run_*.py     ← 在 WSL2 中由 python3 -m 启动
        ↓
PySCF SCF / 优化 / 定域构造
        ↓
quality_gate 断言（失败 → SystemExit 或 agree=null）
        ↓
results/Pn/tables/*.json|csv   （含 quality_gate；可入库）
   ├─ 通过 → deliverables/unit/Pn/{VERDICT,report,evidence}
   └─ 失败协议 → results/Pn/invalid_*/（保留审计，禁支撑终裁）
```

**示例（正式加严 P3，曾用于纠正建环 bug 后的定判扫描）**：

```bash
cd /mnt/d/AI4S-OrgChem
PYTHONPATH=. python3 -m src.p3_nba.run_tight
# 产出：results/P3/tables/pes_tight_B3LYP_6-31gs_on_RHF_3-21g_both.{csv,json}
# 作废几何绝不可用：results/P3/invalid_wrong_geometry/
```

**示例（无参数 ESE 基准，项目最强命中之一）**：

```bash
PYTHONPATH=. python3 -m src.p6_ese.run_v4_objections --basis '6-31g*'
# 苯 ESE ≈ −35.44 kcal/mol；环丁二烯 vert@G* ΔEA ≈ +53.98 kcal/mol
```

### 3.6 计算纪律（在 WSL2 上同样强制）

1. **正确分子优先**：拓扑断言失败则不算能量、不定判（P3 事故教训）。
2. **收敛可辩护**：达优化阈，或有独立代理（如 P3 双向 PES 滞后 ≤ 1 kcal/mol）。
3. **能量尺度窗**：ΔE 离谱（如构象差 10² kcal）→ 拒收本轮，不写成“惊人发现”。
4. **路径洁净**：正式表只在 `results/Pn/tables/`；坏数据进 `invalid_*/`。
5. **判定唯一落盘**：聊天结论无效；只认 `deliverables/unit/Pn/VERDICT.md`。

---

## 4. 本目录内容（归档结构）

```
WSL2/
├── README.md                 # 本说明（含详细计算介绍）
├── MANIFEST.txt              # 快照文件清单与时间戳
├── inventory/                # 只读环境指纹
│   ├── uname.txt
│   ├── os-release.txt
│   ├── nvidia-smi.txt        # RTX 4060 Laptop · CUDA 12.7（快照时）
│   ├── python_system_pyscf.txt
│   ├── pip_freeze_user.txt
│   ├── pip_show_pyscf_stack.txt
│   ├── path_map.txt
│   ├── mamba_ai4orgchem_POINTER.txt
│   └── …
├── home_cascade/             # 自 /home/wsluser 拷贝的小体积目录
│   ├── p3run/                # 早期 P3 运行用的 src 切片
│   ├── ai4orgchem_env/       # 小型 venv 壳（指向 /usr 的符号链接已改为说明文件）
│   └── ai4orgchem-scratch/
└── mounts/
    └── D_drive_project.md
```

### 刻意未拷贝的内容

| 路径（WSL） | 约大小 | 原因 |
|-------------|--------|------|
| `/home/wsluser/micromamba/envs/ai4orgchem` | ~8.9 GB | 体积过大；以 `inventory/mamba_*` 代替 |
| `/home/wsluser/micromamba/envs/ai4orgchem-nequip` | ~6.6 GB | 非本项目主证据引擎 |
| `/home/wsluser/.local/lib/.../pyscf` 整树 | ~190 MB | 本机已装；版本见 inventory |

---

## 5. 在 WSL2 中复现（操作摘要）

```bash
# Windows：wsl -d Ubuntu-24.04
cd /mnt/d/AI4S-OrgChem

# 只读自检（不改环境）
bash tools/check_wsl_env.sh
nvidia-smi

# 跑一条定判级作业（示例）
PYTHONPATH=. python3 -m src.p6_ese.run_v4_objections --basis '6-31g*'
```

`quality_gate.passed != true` 时脚本应非零退出或将 `agree` 置空，**不得**暗示可写 VERDICT。

环境策略：**不擅自 pip/apt 升级**；缺包则记录缺口并询问，不得自行安装。

---

## 6. 刷新本快照

在 WSL 内执行：

```bash
bash /mnt/d/AI4S-OrgChem/tools/wsl_snapshot.sh
```

会重写 `inventory/`、同步 `home_cascade/` 小目录，并更新 `MANIFEST.txt`（仍**不会**拷贝多 GB micromamba 环境）。

---

## 7. 相关文档

| 文档 | 内容 |
|------|------|
| [`docs/env_setup.md`](../docs/env_setup.md) | 环境只读核查 |
| [`docs/quality_gates.md`](../docs/quality_gates.md) | G1–G5 |
| [`results/README.md`](../results/README.md) | 各命题定判文件名 |
| [`deliverables/final/FULL_REPORT.md`](../deliverables/final/FULL_REPORT.md) | 方法总述 + 全套复现命令 |
| [`src/README.md`](../src/README.md) | 自研包布局与定判入口 |
