# 量子化学专家快速审阅指南  
# Quick Review Guide for Quantum-Chemistry Experts

> **读者**：熟悉共轭/芳香性/能量分解，但不关心 AI 编程细节的专家审阅者。  
> **Audience**: Experts in conjugation/aromaticity/energy decomposition who do **not** need AI or software-engineering details.  
> **本仓库 GitHub**：`https://github.com/CCcolab/AI4S-OrgChem`（Owner: Xiao Chen / CCcolab）  
> **建议用时**：核心论证 **15–30 分钟**；单命题深读 **+10 分钟/条**。

---

## 0. 先读这三句 / Read these three sentences first

| 中文 | English |
|------|---------|
| 这是对虞忠衡专著中 **9 条可检验化学主张** 的 **第三方独立复现**，计算引擎为 **PySCF**（B3LYP/6-31G* 为主），**不是**对原著程序的复算。 | This is an **independent third-party replication** of **nine testable chemical claims** from Yu’s monograph, using **PySCF** (mainly B3LYP/6-31G*), **not** a rerun of the author’s programs. |
| **「一致 / 非一致」= 在原著公开定义下能否复现符号与量级**，不等于「原著物理上全对」或「传统理论已被推翻」。 | **“Agree / Disagree”** means whether **signs and magnitudes are reproducible under Yu’s published definitions**—not that Yu is universally “correct” or that textbook theory is “overturned.” |
| 本支路结果：**8 一致 + 1 非一致（P1）**。唯一非一致条目反而提高可信度——全票通过反而像背书。 | This replication arm: **8 agree + 1 disagree (P1)**. The single dissent increases credibility; unanimous agreement would look like endorsement. |

**可忽略的部分（审阅论证时可跳过）**  
**Safe to skip when reviewing the science**

- AI Agent、Cursor、三模型计划等 **方法论包装** → 见 [`deliverables/papers/`](../deliverables/papers/)（非化学主证据）  
- `WSL2/inventory/` 环境指纹、`implementation/` 中文计划表（预注册细节）→ 仅当您要审计流程时再打开  
- 仓库 **不含** 原著 PDF；需自行依法获取后对照公开定义  

---

## 1. 15 分钟阅读路线 / 15-minute reading path

按顺序打开下列文件即可把握 **全项目论证**，无需读代码。

| 步骤 | 文件 | 您将得到什么 |
|------|------|----------------|
| **①** | [`deliverables/final/EXECUTIVE_SUMMARY.md`](../deliverables/final/EXECUTIVE_SUMMARY.md) | 一页结论：8+1、三条最硬证据、三条限定 |
| **②** | [`deliverables/final/VERDICT_TABLE.md`](../deliverables/final/VERDICT_TABLE.md) | 九命题 **预注册判据 ↔ 实测值 ↔ 是否达阈**（审计表） |
| **③** | 任选 1–2 条您最关心的 `deliverables/unit/Pn/VERDICT.md` | 单命题正式判定 + 判据对照表 |
| **④** | 对应 `deliverables/unit/Pn/report.md` | 方法、敏感性、**专家异议预案** |
| **⑤** | 对应 `deliverables/unit/Pn/evidence/` 或 [`deliverables/final/evidence_pack/`](../deliverables/final/evidence_pack/) | 关键 JSON/CSV 数值（可下载核对） |

**English — same path**

| Step | File | What you get |
|------|------|----------------|
| **①** | `deliverables/final/EXECUTIVE_SUMMARY.md` | One-page scoreboard, three strongest lines of evidence, caveats |
| **②** | `deliverables/final/VERDICT_TABLE.md` | Pre-registered criteria vs. results for all nine propositions |
| **③** | `deliverables/unit/Pn/VERDICT.md` (pick 1–2) | Formal verdict for one proposition |
| **④** | `deliverables/unit/Pn/report.md` | Methods, sensitivities, anticipated objections |
| **⑤** | `deliverables/unit/Pn/evidence/` or `deliverables/final/evidence_pack/` | Key numbers in JSON/CSV |

---

## 2. 判定用语（化学审阅必知）/ Verdict vocabulary

| 中文 | English | 含义 / Meaning |
|------|---------|----------------|
| **一致** | **Agree** | 在**该命题预先写明的**定义、参考态、方法窗内，独立 PySCF 实现复现了原著主张的**符号与关键量级**。 |
| **非一致** | **Disagree** | 符号、量级或关键结论与命题判据**相悖**，或在该定义下**不可复现**。 |
| **待定 / 论证中** | **Pending** | 证据闸未过（如未收敛、分子拓扑错误）→ **不得**写一致/非一致；不是「验证失败」。 |
| **可复现性判定** | **Replicability verdict** | 不是「哪套理论正确」的哲学裁决。 |
| **quality_gate** | **quality_gate** | 结果 JSON 中的质控字段；`passed: true` 才可用于判定。 |

**符号约定（与原著一致）**  
**Sign convention (same as Yu)**

- 稳定化能：**负**；去稳定化能：**正**  
- Stabilization: **negative**; destabilization: **positive**  
- 报告单位：**kcal/mol**（内部 Hartree）

---

## 3. 九命题一张表（专家速查）/ Nine propositions at a glance

| ID | 化学问题（一句话） | One-line chemical question | 判定 Verdict | 专家优先看什么 Look at first |
|----|-------------------|----------------------------|--------------|------------------------------|
| **P1** | 丁二烯共轭能：换参考分子是否翻转符号？ | Butadiene CE: sign flip when reference molecule changes? | **非一致 Disagree** | `evidence/ce_summary.csv`；四层次均未翻转 |
| **P2** | 共轭是否去稳定且驱动畸变（总命题）？ | Conjugation destabilizes and drives distortion? | **一致 Agree** | 汇总 P3/P5/P6 + LFMO-lite |
| **P3** | NBA 大扭转角是否能量最低？ | NBA: crowded large twist as energy minimum? | **一致 Agree** | `pes_tight_*.csv`；θ_min≈44.9° |
| **P4** | 苯 D₆ₕ 是否核排斥主导？ | Benzene D₆h: nuclear repulsion driven? | **一致 Agree** | BLA 扫描；ΔEN vs ΔEe |
| **P5** | 局部 ΔEAm 是否恒正且单键伸长？ | Local ΔEAm always positive + bond lengthening? | **一致 Agree** | `p5_v5c_*.json`；6/6 为正 |
| **P6** | 无参数 ESE：苯 ≈−36、CBD ≈+54？ | Parameter-free ESE benchmarks? | **一致 Agree** | 苯 −35.44；CBD +53.98 kcal/mol |
| **P7** | 张力芳香 BLA 是否 π 离域而非角张力？ | Strained aromatic BLA: π not angle strain? | **一致 Agree** | Δr(G)→Δr(PLG) 坍塌 |
| **P8** | 呋喃类是否应称 LDE 而非苯类 ESE？ | Furan-like: LDE not benzene-like ESE? | **一致 Agree** | 呋喃/苯 ΔEA **符号相反** |
| **P9** | 大环 annulene 是否趋于多烯？ | Large [N]annulene → polyene-like? | **一致 Agree** | VDE/π gap 随 N 收敛 |

**原著映射**见 [`docs/propositions.md`](propositions.md) §0。  
**Mapping to the book** → `docs/propositions.md` §0.

---

## 4. 三条最值得先核对的硬证据 / Three lines of evidence to check first

### 4.1 P6 — 无参数双基准（最强）

| 中文 | English |
|------|---------|
| 苯 **ESE = −35.44** kcal/mol（原著约 −36.3） | Benzene **ESE = −35.44** kcal/mol (Yu ≈ −36.3) |
| 环丁二烯垂直 **ΔEA = +53.98**（原著约 +53.6） | Cyclobutadiene vertical **ΔEA = +53.98** (Yu ≈ +53.6) |
| **同一套** 2007 定域 + 虚拟参考流程，**两个相反符号同时命中** | **Same** 2007 localization + VR protocol; **both sign directions hit** |

文件：`deliverables/final/evidence_pack/P6_v4.json` · 报告 `unit/P6/report.md`

---

### 4.2 P7 — 对 Mills–Nixon 的直接检验

| 中文 | English |
|------|---------|
| 完整几何：中心环 BLA **Δr ≈ +0.207 Å** | Full geometry: central ring BLA **Δr ≈ +0.207 Å** |
| 切断中心–外周 π 耦合（PLG）后：**Δr ≈ +0.020 Å** | After cutting center–periphery π coupling (PLG): **Δr ≈ +0.020 Å** |
| 若角张力主导，切断 π 不应使 BLA 几乎消失 | If angle strain dominated, cutting π should **not** collapse BLA |

文件：`evidence_pack/P7_v1c.json` · `unit/P7/report.md`

---

### 4.3 P3 — 构象反例（冲击「共平面最稳」）

| 中文 | English |
|------|---------|
| 约束弛豫 PES 最低：**θ ≈ 44.9°** | Constrained-relaxation PES minimum: **θ ≈ 44.9°** |
| 无约束 B3LYP 自由优化：**≈ 34.5°**（与晶体学量级一致） | Unconstrained B3LYP: **≈ 34.5°** (crystallographic range) |
| **勿引用** `results/P3/invalid_wrong_geometry/`（建环错误，已作废） | **Do not cite** `invalid_wrong_geometry/` (wrong topology; voided) |

文件：`evidence_pack/P3_pes_tight.csv` · `unit/P3/report.md`

---

## 5. 唯一非一致：P1 该怎么读 / The one dissent: how to read P1

| 中文 | English |
|------|---------|
| 原著入口叙事：用 **1-丁烯** 作参考得 CE 稳定化，用 **trans-2-丁烯** 得 CE 去稳定化（符号翻转） | Yu’s entry narrative: **1-butene** reference → stabilizing CE; **trans-2-butene** → destabilizing CE (sign flip) |
| 本项目（B3LYP/6-31G* 等）：**CE₁ 与 CE₂ 均为正**（约 +8～+9 与 +1.9 kcal/mol），**无翻转** | This work: **both CE₁ and CE₂ positive**; **no flip** |
| 自研 2007-GL：ΔE **+4.06** kcal/mol，**去稳定同号**，但未落入书中 **+1.4** 容差（2014 法未实现） | Homemade 2007-GL: ΔE **+4.06**; same sign (destabilizing) but not Yu’s **+1.4** (2014 method not implemented) |
| **解读**：否的是**热化学入口叙事**的可复现性；**不**等于否定 P5/P6 定域层面的去稳定结论 | **Reading**: rejects **thermochemical entry narrative** replicability; does **not** negate delocalization destabilization in P5/P6 |

文件：`unit/P1/VERDICT.md` · `evidence/ce_summary.csv`

---

## 6. 方法学最小信息（够审阅用）/ Minimum methodology (enough to review)

| 项 Item | 内容 Content |
|---------|----------------|
| 软件 Software | **PySCF 2.14.0**（WSL2 / Ubuntu 24.04 上运行；见 [`WSL2/README.md`](../WSL2/README.md)） |
| 主层次 Main level | **B3LYP / 6-31G*** |
| 定域 Localization | 2007 法：π-AO 识别 + 片段间 π–π **Fock/重叠块删除**；GE-m 超片段；ESE/VR；PLG；杂环 hetero_gl |
| 未实现 Not implemented | 严格 **2011/2014 交换积分删除**（用「2011-lite」仅作符号检验） |
| 非主证据 Not primary evidence | 机器学习势、符号回归、聊天结论 |
| 敏感性 Sensitivities | P1：RHF/MP2/ZPE；P8：双基组；P9：2011-lite、BLA 几何 |

**能量分解的参考态依赖**是领域共识弱点；本项目把所有口径写入脚本与 JSON，并在各 `report.md` 中声明——审阅时请重点看 **参考态定义** 是否与您接受的口径一致。

**Reference-state dependence** in energy decomposition is a known issue; all protocols are documented in scripts/JSON and `report.md`—please check whether the **reference-state definition** matches what you would accept.

---

## 7. 质控：如何判断数字能否信 / Quality control: can you trust the numbers?

正式判定前须过 **五道闸 G1–G5**（详见 [`docs/quality_gates.md`](quality_gates.md)）：

| 闸 Gate | 化学含义 Chemical meaning |
|---------|---------------------------|
| **G1** | 算的是**对的分子**（键连、环尺寸正确） |
| **G2** | 几何健全（无严重原子重叠、键长合理） |
| **G3** | 优化收敛，或有可辩护代理（如 P3 双向扫描滞后 ≤1 kcal/mol） |
| **G4** | 能量差在化学合理窗（如构象差不应达 10² kcal/mol） |
| **G5** | 作废数据在 `invalid_*/`，**未**混入正式表 |

**审阅技巧**：打开 JSON，查 `"quality_gate": { "passed": true }`。若为 `false` 或缺失，该点**不得**支撑 VERDICT。

**Review tip**: In JSON, check `"quality_gate": { "passed": true }`. If false or missing, do not use for verdicts.

**已知作废区（可审计，禁止引用为终裁）**  
**Known voided data (auditable; do not cite for final verdicts)**

- `results/P3/invalid_wrong_geometry/` — 建环错误，曾几乎误判「非一致」

---

## 8. 常见误读（请避免）/ Common misreadings (please avoid)

| 误读 Misreading | 正确理解 Correct understanding |
|-----------------|--------------------------------|
| 「8 一致 = 原著全对、教科书全错」 | 仅表示 **在 Yu 公开定义下可复现**；传统理论有多条独立证据线未在本项目中检验 |
| 「1 非一致 = 项目失败」 | P1 非一致是 **刻意保留的审计结果**；定域证据仍支持多条「一致」 |
| 「一切共轭都去稳定」 | **错误**。苯 **ESE<0**（P6）与局部 **ΔEAm>0**（P5）是原著**分口径**的两件事 |
| 「AI 算了数，不可信」 | 判定数字来自 **PySCF**；AI 负责编排与写脚本，**聊天不算证据** |
| 引用 `invalid_*` 或旧版 pilot JSON | 以各命题 VERDICT 与 **定判文件** 表为准（如 P5→v5c，P6→v4） |

---

## 9. 若您要反驳或追问 / If you want to challenge or dig deeper

建议按下列顺序提出异议（与 [`unit/Pn/report.md`](../deliverables/unit/) 中「专家异议预案」对应）：

1. **参考态 / 分解定义** — GL vs GE-m vs ESE vs LDE 是否公平？  
2. **基组与泛函** — 6-31G* / B3LYP 是否足够？（见各命题敏感性表）  
3. **几何模型** — 垂直 vs 半绝热；平面 Kekulé 大环（P9）  
4. **收敛与 PES** — P3 滞后代理是否可接受？  
5. **P1 热化学** — 是否应改用实验氢化热（本项目故意绑定 **独立 QC** 判据）

**深读单命题**时打开：

- 判据来源：[`docs/propositions.md`](propositions.md)  
- 预注册阈值：[`implementation/`](../implementation/) 中 `Pn_*计划*.md`  
- 复现命令：[`deliverables/final/FULL_REPORT.md`](../deliverables/final/FULL_REPORT.md) §7（需 WSL + PySCF，**不必**为审阅而重算）

---

## 10. 三条口径限定（写进任何引用里）/ Three caveats for any citation

1. **可复现性 ≠ 物理唯一性** · Replicability ≠ unique physical truth  
2. **能量分解参考态依赖** · Reference-state dependence in decomposition  
3. **本仓库 = Cursor Grok 一支**；三模型交叉验证计划进行中 · This repo = one of three planned independent LLM arms  

---

## 11. 文件地图（只看化学）/ File map (chemistry only)

```
README.md                          ← 项目背景（含 AI 说明，可略读）
deliverables/final/
  EXECUTIVE_SUMMARY.md               ← ★ 先读
  VERDICT_TABLE.md                   ← ★ 审计总表
  FULL_REPORT.md                     ← 方法与异议专章
  evidence_pack/                     ← ★ 九命题主数值
deliverables/unit/P1…P9/
  VERDICT.md                         ← ★ 正式判定
  report.md                          ← 方法 + 异议回应
  evidence/                          ← 精选表
docs/
  propositions.md                    ← 九命题权威规格
  quality_gates.md                   ← 质控闸（可选）
  expert_quick_review_guide.md       ← 本文件
src/                                 ← 复现脚本（仅重算时需要）
results/Pn/tables/                   ← 完整计算表
WSL2/README.md                       ← 计算在 WSL2+PySCF 上完成（可选）
```

---

## 12. 引用 / Citation

见仓库根目录 [`CITATION.cff`](../CITATION.cff)。  
See [`CITATION.cff`](../CITATION.cff) at repo root.

---

*文档版本 Document version：2026-08-27 · 对齐结题冻结 Aligned with freeze 2026-08-25*
