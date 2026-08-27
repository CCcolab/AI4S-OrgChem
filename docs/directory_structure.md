# 目录结构（GitHub 公开口径）

> **AI4S-OrgChem** 仓库目录职责与落盘约定。  
> 原则：**源码可复现 · 结果可追溯 · 判定可审计 · 版权材料不入库**。  
> 本文件同时是 GitHub 读者看到的结构说明书。

---

## 一、GitHub 公开目录树（权威）

下列树即为推送到 GitHub 后应呈现的结构。`source/`、`data/` 下**只有说明 README**；PDF 与全书抽文本由 `.gitignore` 排除，读者自行网上检索原著。

```
AI4S-OrgChem/
│
├── README.md                 # English homepage（GitHub 默认打开）
├── README.zh-CN.md           # 中文首页
├── NOTICE                    # 版权与独立实现声明
├── LICENSE                   # MIT
├── CITATION.cff              # 学术引用元数据
├── .gitignore                # 排除 PDF / 抽文本 / raw / __pycache__ / .cursor 等
│
├── docs/                     # 规范与知识（权威文档）
│   ├── propositions.md       # P1–P9 命题权威规格
│   ├── quality_gates.md      # ★ 证据质控硬闸 G1–G5
│   ├── research_plan.md      # 实施方案
│   ├── directory_structure.md# 本文件
│   ├── book_overview.md      # 原著论断地图（无全书文本）
│   ├── env_setup.md          # 环境只读核查记录
│   └── github_upload_plan.md # GitHub 发布方案
│
├── implementation/           # 预注册计划表与实施跟踪（原「实施/」）
│   ├── README.md
│   ├── 命题总表.md
│   ├── 环境与论证工具清单.md
│   ├── AI4S工作流映射.md
│   ├── P1_计划表.md / P1_计算与论证表.md
│   └── P2…P9_*_计划实施表.md / *_计算与论证表.md
│
├── src/                      # 自研可复现代码（不复制原著程序）
│   ├── README.md
│   ├── common/               # 单位、I/O、能量分量、几何
│   ├── localization/         # GL / GE-m / ESE / PLG / LFMO-lite 等
│   ├── p1_butadiene/
│   ├── p2_aggregate/
│   ├── p3_nba/
│   ├── p4_benzene/
│   ├── p5_local/
│   ├── p6_ese/
│   ├── p7_strained/
│   ├── p8_furan/
│   └── p9_annulene/
│
├── results/                  # 单元计算（可重跑）；大文件 raw 不入库
│   ├── README.md
│   ├── P1/ … P9/
│   │   ├── tables/           # CSV/JSON（含 quality_gate）★ 入库
│   │   ├── figures/          # 静态图（若有）
│   │   ├── summary*.md
│   │   ├── invalid_*/        # 作废数据隔离（可审计，禁支撑 VERDICT）
│   │   └── raw/              # 日志/chk — gitignore，仅本地
│   └── _scratch/             # 临时试算 — gitignore
│
├── deliverables/             # ★ 正式交付（对外结题）
│   ├── README.md
│   ├── unit/                 # 每命题一套
│   │   ├── INDEX.md
│   │   └── P1/ … P9/
│   │       ├── VERDICT.md    # 唯一正式判定：一致 / 非一致
│   │       ├── report.md
│   │       ├── evidence/     # 精选表图
│   │       └── canvas_link.md
│   ├── final/                # 全项目结题五件套（已冻结）
│   │   ├── README.md
│   │   ├── EXECUTIVE_SUMMARY.md
│   │   ├── VERDICT_TABLE.md
│   │   ├── FULL_REPORT.md
│   │   ├── checklist.md
│   │   └── evidence_pack/
│   └── papers/               # 学术产出骨架
│       ├── chem_preprint_outline.md
│       └── ai4s_method_outline.md
│
├── tools/                    # 运维/探测（非科学主逻辑）
│   ├── check_wsl_env.sh
│   ├── inventory_env.sh
│   ├── wsl_snapshot.sh       # 刷新 WSL2/ 归档
│   ├── extract_pdf.py        # 仅本地抽文本；输出不得提交
│   └── …
│
├── WSL2/                     # ★ WSL2 计算平面归档（Ubuntu + GPU）
│   ├── README.md             # 双平面说明：/mnt/d ↔ D:\ 与 /home/cascade
│   ├── inventory/            # OS / nvidia-smi / PySCF 指纹（非全量 env）
│   ├── home_cascade/         # 自 Linux 本机盘拷贝的小目录（p3run 等）
│   └── mounts/               # 挂载关系说明
│
├── source/
│   └── README.md             # 说明：PDF 仅本地、不入库
│
└── data/
    └── README.md             # 说明：抽文本仅本地、不入库
```

### 明确不进 GitHub 的内容

| 路径 / 类型 | 原因 |
|-------------|------|
| `source/*.PDF` | 原著版权材料；读者自行检索下载 |
| `data/book_full.txt`、`data/book_text/` | 全书派生文本，等同再分发 |
| `results/**/raw/**`、`*.chk`、`*.molden` 等 | 大体量中间文件，交付不依赖 |
| `results/_scratch/` | 临时试算 |
| `.cursor/`、`.specstory/`、`*.canvas.tsx` | IDE/Agent 本地；Canvas 由 `canvas_link.md` 指向 |
| `__pycache__/`、`*.pyc` | 编译缓存 |
| `/home/.../micromamba/envs/ai4orgchem` 全量 | ~8.9 GB；仅在 `WSL2/inventory/` 留包清单与路径指针 |

---

## 二、目录职责速查

| 目录 | 入库? | 放什么 | 不放什么 |
|------|-------|--------|----------|
| `docs/` | ✅ | 命题规格、门禁、方案、结构说明 | 大数值结果、全书文本 |
| `implementation/` | ✅ | 预注册计划表、总表、环境/工作流清单 | 正式判定正文（判定在 deliverables） |
| `src/` | ✅ | 自研脚本与定域库 | 原著代码、大输出 |
| `results/Pn/tables/` | ✅ | 可复现表（含 `quality_gate`） | — |
| `results/Pn/invalid_*/` | ✅ | 作废轨迹（禁引用为终裁） | — |
| `results/Pn/raw/` | ❌ | 本机日志/检查点 | 不得 push |
| `deliverables/unit/` | ✅ | VERDICT + report + evidence | 全量 raw |
| `deliverables/final/` | ✅ | 结题五件套 | 半成品判定 |
| `deliverables/papers/` | ✅ | 预印本/方法篇骨架 | — |
| `tools/` | ✅ | 环境探测、WSL2 快照、PDF 抽取 | 科学主流程 |
| `WSL2/` | ✅ | WSL2 计算平面归档（inventory + 小体积 home 拷贝） | 多 GB micromamba 全量 env |
| `source/`、`data/` | 仅 README | 本机 PDF / 抽文本 | **PDF 与全文不得入库** |

---

## 三、阅读顺序（给 GitHub 访客）

```
README.md
  → deliverables/final/EXECUTIVE_SUMMARY.md     # 一页结论
  → deliverables/final/VERDICT_TABLE.md         # 九命题审计表
  → deliverables/unit/Pn/VERDICT.md             # 单命题判定
  → docs/quality_gates.md + src/                # 如何复现
  → implementation/命题总表.md                  # 预注册与跟踪
```

数据流（落盘时序）：

```
docs/propositions.md
  → implementation/Pn_*
  → src/
  → results/Pn/tables/
  → deliverables/unit/Pn/
  → deliverables/final/
```

---

## 四、单元交付（Unit）——P1–P9

### 4.1 完成定义（DoD）

| # | 交付项 | 位置 |
|---|--------|------|
| 1 | 可复现计算 | `results/Pn/` |
| 2 | 自研脚本入口 | `src/…` |
| 3 | G1–G5 通过 | `results/Pn/tables/*` 中 `quality_gate` |
| 4 | **VERDICT.md**（仅一致/非一致；禁引 `invalid_*/`） | `deliverables/unit/Pn/` |
| 5 | `report.md` | 同上 |
| 6 | `evidence/` | 同上 |
| 7 | `canvas_link.md` | 同上（Canvas 本身在 Cursor 本地） |
| 8 | 回填 | `implementation/Pn_*`、`命题总表.md`、`docs/propositions.md` |

### 4.2 命名约定

| 类型 | 约定 |
|------|------|
| 结果 / 交付目录 | `results/P1` … `P9`；`deliverables/unit/P1` …（大写 P） |
| 作废归档 | `results/Pn/invalid_<原因>/` |
| 脚本包 | `src/p1_butadiene/` 等小写+下划线 |
| 实施跟踪目录 | **`implementation/`**（统一英文；不再使用 `实施/`） |

---

## 五、最终交付（Final）

启动条件：P1–P9 均有 `VERDICT.md` 且 G1–G5 全过（已满足，2026-08-25 冻结）。

| 文件 | 内容 |
|------|------|
| `EXECUTIVE_SUMMARY.md` | AI4S Agent 定位 + 一页结论 |
| `VERDICT_TABLE.md` | 纯审计判定总表 |
| `FULL_REPORT.md` | 独立复现总报告 |
| `checklist.md` | 结题核查 |
| `evidence_pack/` | 跨命题主证据 |
| `README.md` | 本区说明 |

```
results/Pn  →  deliverables/unit/Pn  →  deliverables/final
```

数字必须能回溯到 `results/Pn/tables/` 或 `evidence/`；禁止只写在聊天里。

---

## 六、落盘状态（冻结核对）

| 路径 | 状态 |
|------|------|
| 根：`README.md` / `README.zh-CN.md` / `NOTICE` / `LICENSE` / `CITATION.cff` / `.gitignore` | ✅ |
| `docs/`、`implementation/`、`src/`、`tools/` | ✅ |
| `results/P1–P9/tables/` + P3 `invalid_wrong_geometry/` | ✅ |
| `deliverables/unit/P1–P9/` 四件套 | ✅ 9/9 |
| `deliverables/final/` 五件套 | ✅ 已冻结 |
| `deliverables/papers/` | ✅ 两篇骨架 |
| `source/README.md`、`data/README.md` | ✅ 仅说明 |
| `docs/formalism.md` | ⬜ 未建（不阻断结题；公式散见 report 与 `src/localization/`） |

---

## 七、维护

1. **新增顶层目录前先改本文件**，再改 `README.md` / `implementation/README.md` / `github_upload_plan.md`。
2. 大文件进 `results/Pn/raw/` 并由 `.gitignore` 忽略；交付只带精选。
3. 「一致 / 非一致」只以 `deliverables/unit/Pn/VERDICT.md` 与 `deliverables/final/` 为准。
4. 本地可有 PDF 与抽文本；**任何 push 前**确认 `git status` 无 `.PDF`、无 `book_full`、无 `book_text/`。
