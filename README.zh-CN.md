# AI4S-OrgChem

[English](README.md) | **中文**

**本项目是采用最新 AI 技术的 AI for Science（AI4S）Agent 项目**：由智能体在预注册判据与硬证据门禁下，端到端完成命题拆解、独立编码、量子化学计算编排、质控与正式判定落盘——而不是用聊天模型「口头评论」原著。

**科学任务**：对**虞忠衡教授** **《Questioning Fundamental Principles of Organic Chemistry》**中与传统有机结构理论相悖的九条核心命题（P1–P9）做第三方独立复现，每条交付 **一致 / 非一致**；计算主引擎为既有环境中的 **PySCF**。

> **本仓库不提供原著 PDF 或全书抽文本。** 读者请自行在网上依法检索、获取原著。本地工作区中的 `source/`、`data/` **永不上传** GitHub。

---

## 项目背景

共轭稳定化与芳香性能量判据是有机结构理论的基石表述之一。**虞忠衡教授**专著系统质疑其中若干「基本原理」——例如：局部/成对共轭在定域参考下呈去稳定、几何畸变可归因于 π 而非张力、大环 annulene 趋于多烯、呋喃类不宜按苯类 ESE 判芳香等。这些主张高度依赖作者自研的定域化能量分解，**此前缺乏按公开定义独立实现的第三方检验**。

本项目将上述挑战转化为 AI4S Agent 可执行的九条命题（P1–P9）：计算前写死判据与阈值（预注册），仅依据 PDF 中公开的数学定义、轨道图示与能量表达式 **由 Agent 独立编写** PySCF 脚本，经硬证据门禁（拓扑 / 几何 / 收敛 / 尺度 / 路径洁净）后，才允许写入正式判定。

**体裁**：独立复现研究（replication / verification study），嵌在 AI4S Agent 工作流中——不是立场论战，也不是原著推广材料。「一致」表示在原著自身定义下可复现；**不等于**该能量分解是唯一正确的物理图像，**也不等于**传统理论被推翻。

### 三家大模型独立验证计划

整个验证计划设计为 **三大前沿 LLM Agent 各自独立完成全套或对等验证**，互不共享实现代码与中间判定，事后交叉比对，以降低单一模型偏好与「讨好原文」风险——这本身也是 AI4S 可信度实验。

| 验证支路 | 执行环境 | 状态 |
|----------|----------|------|
| **本仓库** | 独立 AI4S Agent 完成 | **已冻结**（2026-08-25）：8 一致 + 1 非一致 |
| 支路 B | 另一家大模型 Agent（独立仓库 / 独立实现） | 计划中 |
| 支路 C | 第三家大模型 Agent（独立仓库 / 独立实现） | 计划中 |

本仓库只代表**本支路**的完整审计轨迹；最终跨模型共识表将在三支齐套后另行发布。

---

---

## 这是什么：AI4S Agent，而不只是「用 AI 写报告」

| 维度 | 本项目做法 |
|------|------------|
| **范式** | **AI for Science（AI4S）**：用前沿大模型 Agent 驱动可审计的科学研究流程 |
| **执行形态** | **Agent 闭环**：读规格 → 写/改代码 → 调 PySCF → 过 G1–G5 门禁 → 写入 `VERDICT`；失败则隔离作废数据并重跑 |
| **技术栈** | 最新一代 LLM Agent（本支路）+ 电子结构计算（PySCF）+ 项目级 Rule/门禁（非一次性提示词） |
| **产出** | 可复现脚本、带 `quality_gate` 的结果表、单元判定书、结题五件套——全部可被他人重跑核验 |
| **刻意不做的** | 不用 ML 势/符号回归代替主证据；不把聊天结论当作正式判定；不复制原著程序 |

科学结论的权威来自 **门禁授权的证据链**；AI4S Agent 的价值在于把这条链跑通、留痕、可复核。

---

---

## 研究价值

1. **AI4S Agent 的可审计科研范式**  
   展示「最新 AI 技术」如何用于**真科研闭环**（编码—计算—质控—判定），而不是生成不可核验的叙述。门禁不过则只能「待定」——给 Agent 合法的「证据不足」出口，压缩自信错误结论的空间。

2. **可复现性优先于对错裁判**  
   先回答「这些主张在公开定义下能否被第三方复现」，再谈理论取舍。对长期依赖单一程序、缺少独立复现的异端框架，这是缺失的一环。

3. **预注册判据 + 硬证据门禁**  
   判据在算前锁定；G1–G5 不过则不得定判。P3 曾因建错分子与欠收敛 PES 几乎写成错误「非一致」，被门禁拦下、坏数据隔离、重跑后翻为「一致」——判定来自证据流程，而非模型口头偏好。

4. **仅凭公开定义由 Agent 独立实现**  
   不接触、不复制原著程序；检验定义层面的可复现性，比「同一代码算出同一数字」更强。

5. **多模型独立验证的实验设计**  
   三家 LLM Agent 分仓独立做同一命题集，为 AI4S 提供可引用的交叉验证案例：比拼的不是文笔，而是门禁下的符号与量级是否一致。

6. **学术呈现的克制边界**  
   明确写出能量分解的参考态依赖、口径限定，以及一条正式 **非一致**（P1）。有异议的复现报告，比全票通过更可信。

---

---

## 最终结论 — 对虞忠衡教授 P1–P9 的独立验证

**在原著公开定义与参考态口径下：**

| | |
|--|--|
| **总判** | **8 一致 · 1 非一致**（2026-08-25 冻结） |
| **非一致** | **P1** — 丁二烯共轭能随参考分子符号翻转：**未复现**（四种理论层次 CE₁ 均为正） |
| **一致** | P2、P3、P4、P5、P6、P7、P8、P9 |

**结论入口（建议先看）：**

| 文档 | 链接 |
|------|------|
| **判定总表（审计）** | [deliverables/final/VERDICT_TABLE.md](deliverables/final/VERDICT_TABLE.md) · [English](deliverables/final/VERDICT_TABLE.en.md) |
| **执行摘要** | [deliverables/final/EXECUTIVE_SUMMARY.md](deliverables/final/EXECUTIVE_SUMMARY.md) · [English](deliverables/final/EXECUTIVE_SUMMARY.en.md) |
| **单元判定 P1–P9** | [deliverables/unit/INDEX.md](deliverables/unit/INDEX.md) · [English](deliverables/unit/INDEX.en.md) |

> **一致 / 非一致 = 在原著口径下能否被第三方复现**，不等于「全书对错」或「教科书理论被推翻」。限定说明见执行摘要。

---

**本项目是采用最新 AI 技术的 AI for Science（AI4S）Agent 项目**：由智能体在预注册判据与硬证据门禁下，端到端完成命题拆解、独立编码、量子化学计算编排、质控与正式判定落盘——而不是用聊天模型「口头评论」原著。

**科学任务**：对**虞忠衡教授** **《Questioning Fundamental Principles of Organic Chemistry》**中与传统有机结构理论相悖的九条核心命题（P1–P9）做第三方独立复现，每条交付 **一致 / 非一致**；计算主引擎为既有环境中的 **PySCF**。

> **本仓库不提供原著 PDF 或全书抽文本。** 读者请自行在网上依法检索、获取原著。本地工作区中的 `source/`、`data/` **永不上传** GitHub。

---

## 这是什么：AI4S Agent，而不只是「用 AI 写报告」

| 维度 | 本项目做法 |
|------|------------|
| **范式** | **AI for Science（AI4S）**：用前沿大模型 Agent 驱动可审计的科学研究流程 |
| **执行形态** | **Agent 闭环**：读规格 → 写/改代码 → 调 PySCF → 过 G1–G5 门禁 → 写入 `VERDICT`；失败则隔离作废数据并重跑 |
| **技术栈** | 最新一代 LLM Agent（本支路）+ 电子结构计算（PySCF）+ 项目级 Rule/门禁（非一次性提示词） |
| **产出** | 可复现脚本、带 `quality_gate` 的结果表、单元判定书、结题五件套——全部可被他人重跑核验 |
| **刻意不做的** | 不用 ML 势/符号回归代替主证据；不把聊天结论当作正式判定；不复制原著程序 |

科学结论的权威来自 **门禁授权的证据链**；AI4S Agent 的价值在于把这条链跑通、留痕、可复核。

---

## 项目背景

共轭稳定化与芳香性能量判据是有机结构理论的基石表述之一。**虞忠衡教授**专著系统质疑其中若干「基本原理」——例如：局部/成对共轭在定域参考下呈去稳定、几何畸变可归因于 π 而非张力、大环 annulene 趋于多烯、呋喃类不宜按苯类 ESE 判芳香等。这些主张高度依赖作者自研的定域化能量分解，**此前缺乏按公开定义独立实现的第三方检验**。

本项目将上述挑战转化为 AI4S Agent 可执行的九条命题（P1–P9）：计算前写死判据与阈值（预注册），仅依据 PDF 中公开的数学定义、轨道图示与能量表达式 **由 Agent 独立编写** PySCF 脚本，经硬证据门禁（拓扑 / 几何 / 收敛 / 尺度 / 路径洁净）后，才允许写入正式判定。

**体裁**：独立复现研究（replication / verification study），嵌在 AI4S Agent 工作流中——不是立场论战，也不是原著推广材料。「一致」表示在原著自身定义下可复现；**不等于**该能量分解是唯一正确的物理图像，**也不等于**传统理论被推翻。

### 三家大模型独立验证计划

整个验证计划设计为 **三大前沿 LLM Agent 各自独立完成全套或对等验证**，互不共享实现代码与中间判定，事后交叉比对，以降低单一模型偏好与「讨好原文」风险——这本身也是 AI4S 可信度实验。

| 验证支路 | 执行环境 | 状态 |
|----------|----------|------|
| **本仓库** | 独立 AI4S Agent 完成 | **已冻结**（2026-08-25）：8 一致 + 1 非一致 |
| 支路 B | 另一家大模型 Agent（独立仓库 / 独立实现） | 计划中 |
| 支路 C | 第三家大模型 Agent（独立仓库 / 独立实现） | 计划中 |

本仓库只代表**本支路**的完整审计轨迹；最终跨模型共识表将在三支齐套后另行发布。

---

## 研究价值

1. **AI4S Agent 的可审计科研范式**  
   展示「最新 AI 技术」如何用于**真科研闭环**（编码—计算—质控—判定），而不是生成不可核验的叙述。门禁不过则只能「待定」——给 Agent 合法的「证据不足」出口，压缩自信错误结论的空间。

2. **可复现性优先于对错裁判**  
   先回答「这些主张在公开定义下能否被第三方复现」，再谈理论取舍。对长期依赖单一程序、缺少独立复现的异端框架，这是缺失的一环。

3. **预注册判据 + 硬证据门禁**  
   判据在算前锁定；G1–G5 不过则不得定判。P3 曾因建错分子与欠收敛 PES 几乎写成错误「非一致」，被门禁拦下、坏数据隔离、重跑后翻为「一致」——判定来自证据流程，而非模型口头偏好。

4. **仅凭公开定义由 Agent 独立实现**  
   不接触、不复制原著程序；检验定义层面的可复现性，比「同一代码算出同一数字」更强。

5. **多模型独立验证的实验设计**  
   三家 LLM Agent 分仓独立做同一命题集，为 AI4S 提供可引用的交叉验证案例：比拼的不是文笔，而是门禁下的符号与量级是否一致。

6. **学术呈现的克制边界**  
   明确写出能量分解的参考态依赖、口径限定，以及一条正式 **非一致**（P1）。有异议的复现报告，比全票通过更可信。

---

## 项目成就（本支路）

### 总判定（2026-08-25 冻结）

**8 一致 + 1 非一致。**

| 判定 | 命题 |
|------|------|
| **一致** | P2, P3, P4, P5, P6, P7, P8, P9 |
| **非一致** | **P1**（丁二烯氢化热参考态符号翻转：四层次下 CE₁ 均为正，未出现翻转） |

一页摘要 · 总表 · 总报告：

- [deliverables/final/EXECUTIVE_SUMMARY.md](deliverables/final/EXECUTIVE_SUMMARY.md)
- [deliverables/final/VERDICT_TABLE.md](deliverables/final/VERDICT_TABLE.md)
- [deliverables/final/FULL_REPORT.md](deliverables/final/FULL_REPORT.md)
- [deliverables/final/checklist.md](deliverables/final/checklist.md)

### 三条最硬的正面证据

| # | 命题 | 结果要点 |
|---|------|----------|
| 1 | **P6** | 无参数 ESE：苯 **−35.44** kcal/mol（原著约 −36.3）；环丁二烯垂直 ΔEA **+53.98**（原著约 +53.6）——**两个反向基准同时命中** |
| 2 | **P7** | 张力芳香：切断中心–外周 π 耦合后，中心环 BLA 由 **+0.207 Å → +0.020 Å**，与 Mills–Nixon 角张力主导预期相悖 |
| 3 | **P3** | NBA 弛豫 PES 最低点 **θ ≈ 44.9°**；无约束自由优化独立落在 **≈34.5°**——冲击「共轭稳定化使共平面最稳」的常见推论 |

### 方法学资产（AI4S Agent 管线）

- **Agent 端到端**：命题规格 → 自研定域代码（`src/localization/`）→ 批量计算 → 门禁 → 判定书，全程可追溯
- 九命题单元齐套：`VERDICT.md` + `report.md` + `evidence/` + 质控闸字段
- 作废数据隔离：`results/P3/invalid_wrong_geometry/`（可审计，禁止支撑终裁）——Agent「算错了却讲得很顺」时的纠错样板
- 学术产出骨架：[deliverables/papers/](deliverables/papers/)（化学篇预印本 + **AI4S 方法篇**：证据门禁与判定反转）

### 必须同时记住的限定

- 判定是 **可复现性** 判定，不是物理唯一性判定。
- 能量分解 **参考态依赖**；口径全部落盘并附敏感性。
- 原著自身分口径：局部共轭去稳定（ΔEAm>0）与苯环额外离域稳定（ESE<0）不可混读成「一切共轭都去稳定」。
- 本仓库 **未** 使用、复制或翻译原著任何程序代码；主证据仅为既有环境中的 **PySCF**。

---

---

## 本仓库包含 / 不包含

| 包含（公开） | 不包含（本地或自行获取） |
|--------------|--------------------------|
| `src/` 自研脚本 | **`source/` 原著 PDF** |
| `docs/`、`implementation/` 规格与预注册表 | **`data/` 全书抽文本** |
| `results/Pn/tables/` 与作废区说明 | 原著附录/程序代码 |
| `deliverables/` 判定、报告、证据包 | Cursor 本机 Canvas 二进制 |
| `tools/` 环境核查等 | 任何需改系统包版本才能跑的依赖假设 |

原著书名可检索：**Questioning Fundamental Principles of Organic Chemistry**（**虞忠衡教授** / Professor Zhongheng Yu）。获取后，公开定义与数值表即可对照本仓库的判据与结果；**无需**本仓库提供电子版。

---

## 关键入口

| 文档 | 说明 |
|------|------|
| [docs/propositions.md](docs/propositions.md) | 命题权威规格 |
| [docs/quality_gates.md](docs/quality_gates.md) | 证据质控硬闸 G1–G5 |
| [**docs/expert_quick_review_guide.md**](docs/expert_quick_review_guide.md) | **量子化学专家快速审阅指南（中英）** |
| [docs/research_plan.md](docs/research_plan.md) | 实施方案 |
| [docs/directory_structure.md](docs/directory_structure.md) | 目录与落盘约定 |
| [WSL2/README.md](WSL2/README.md) | **WSL2 计算平面**（Ubuntu + GPU 归档） |
| [docs/github_upload_plan.md](docs/github_upload_plan.md) | GitHub 发布方案 |
| [docs/RELEASE_NOTES.md](docs/RELEASE_NOTES.md) | **v1.0.0 冻结版 Release 说明** |
| [implementation/](implementation/) | 计划表、环境与工具、AI4S 映射 |
| [NOTICE](NOTICE) | 版权与独立实现声明 |

---

## 目录结构（GitHub）

权威说明：[docs/directory_structure.md](docs/directory_structure.md)

```
AI4S-OrgChem/
├── README.md · README.zh-CN.md · NOTICE · LICENSE · CITATION.cff · .gitignore
├── docs/                 # 命题规格、质控闸、方案、结构说明
├── implementation/       # 预注册计划表与跟踪（英文目录名）
├── src/                  # 自研代码：common / localization / p1…p9
├── results/Pn/           # tables（入库）· invalid_*（审计）· raw（不入库）
├── deliverables/
│   ├── unit/Pn/          # VERDICT · report · evidence · canvas_link
│   ├── final/            # 结题五件套（已冻结）
│   └── papers/           # 化学篇 / AI4S 方法篇骨架
├── WSL2/                 # ★ 计算平面：Ubuntu-24.04 + GPU 归档
├── tools/                # 环境核查 · wsl_snapshot.sh
├── source/README.md      # 仅说明；PDF 不入库
└── data/README.md        # 仅说明；抽文本不入库
```

---

## 复现（环境约定）

**主计算在 WSL2**（Ubuntu 24.04，PySCF 2.14.0，GPU 可见）。Windows 树与 WSL 挂载点为同一目录：`D:\AI4S-OrgChem` ≡ `/mnt/d/AI4S-OrgChem`。详见 [`WSL2/README.md`](WSL2/README.md)。

```bash
# 在 WSL2 内
cd /mnt/d/AI4S-OrgChem
PYTHONPATH=. python3 -m src.p6_ese.run_v4_objections --basis '6-31g*'
```

全套入口命令见 [FULL_REPORT.md §7](deliverables/final/FULL_REPORT.md)。  
`quality_gate.passed != true` 时脚本非零退出，不产出判定建议。

---

## 引用与许可

- **原著作者**：**虞忠衡教授**（Professor Zhongheng Yu）
- **Copyright (c) 2026 Xiao Chen** · e-mail: [chenxiao0101@gmail.com](mailto:chenxiao0101@gmail.com)
- 引用本复现支路：见 [`CITATION.cff`](CITATION.cff)
- 代码许可：见 [`LICENSE`](LICENSE)（MIT）
- 原著版权归**虞忠衡教授**与出版方；本仓库不构成对原著材料的再分发。
