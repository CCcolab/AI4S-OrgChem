# AI4S-OrgChem 独立技术研究实施方案

> **项目定位**：**Cursor-assisted AI4S** 证据支路；以第三方身份用 PySCF 验证虞忠衡教授 P1–P9。  
> **对外状态（L0）**：见 [`FROZEN_VERDICT_AUTHORITY.md`](FROZEN_VERDICT_AUTHORITY.md) · 第三方复审（2026-09-01）：本地 `AI4S_Post_Remediation_Review_2026-09-01.pdf`（§4；**部分通过**）。  
> **L1 审计快照（2026-08-25）**：历史记录；**不是**科学总分。  
> **合规红线**：仅依据 PDF 公开的数学定义、轨道图示、能量表达式与数值表格；不使用、不复制、不反编译、不翻译原著程序代码。  
> **命题权威清单**：[propositions.md](propositions.md)（P1–P9）。  
> **实施跟踪**：[`implementation/`](../implementation/)。  
> **目录与交付落盘**：[`directory_structure.md`](directory_structure.md)。  
> **严格项目边界 Rule**：[`.cursor/rules/ai4s-orgchem-boundaries.mdc`](../.cursor/rules/ai4s-orgchem-boundaries.mdc)（本地 Agent；公开仓可不含 `.cursor/`）。  
> **证据质控硬闸 Rule**：[`.cursor/rules/ai4s-quality-gates.mdc`](../.cursor/rules/ai4s-quality-gates.mdc) · [`quality_gates.md`](quality_gates.md)。  
> **AI4S 工作流参考**：[公众号文章](https://mp.weixin.qq.com/s/v6qcn_LI7ZdcoxB0xIcO3A) → [`implementation/AI4S工作流映射.md`](../implementation/AI4S工作流映射.md)。

---

## 一、已锁定的工程决策

| 项目 | 决策 |
|------|------|
| 引擎 | **PySCF（在 WSL 中）** |
| 环境策略 | **不改动** WSL 中任何已有软件/配置（含 GPU 相关设置） |
| 论证深度 | 交付与原著的 **一致 / 非一致**；须可辩护、可复现、经得起专家追问 |
| 交付形式 | 每条命题 + **Cursor Canvas** 交互可视化报告；另有总看板 |
| 单位约定 | 内部 Hartree；报告 kcal/mol（1 Ha = 627.5095 kcal/mol）；稳定化为负、去稳定为正 |

### 环境只读核查结果（2026-08-23）

- WSL Python：3.12.3（`/usr/bin/python3`）
- PySCF：**2.14.0**（已可用）
- GPU：NVIDIA GeForce RTX 4060 Laptop GPU（`nvidia-smi` 可见）
- **未执行任何安装、升级或配置修改**

---

## 二、命题追踪计划表（权威版见 propositions.md）

判定只允许：**一致** / **非一致** / **待定**（论证中不提前判）。  
**写 VERDICT 前必须过质控闸 G1–G5**（拓扑/正确分子、收敛、能量尺度等）——见 [`quality_gates.md`](quality_gates.md)；严禁未收敛或错误分子上定判。

| ID | 原著命题（摘要） | 论证工具 / 技术 | 状态 | 判定 | 完成度 |
|----|------------------|-----------------|------|------|--------|
| **P1** | 丁二烯共轭能：参考分子选择翻转符号；经典 −3.9 不可作为共轭稳定证据 | 热化学循环；CE₁ vs CE₂；RHF/B3LYP/MP2/ZPE；自研 2007-GL | 本层级判定已闭合 | **非一致** | ~98% |
| **P2** | 共轭去稳定且驱动畸变（**元命题汇总**） | 汇总 P3/P5/P6；LFMO-lite | 本层级判定已闭合 | **一致** | ~94% |
| **P3** | NBA：大扭转角拥挤构象最稳定 | 级联预松弛 + 双向 θ 约束 + B3LYP；拓扑门禁 | 已结案 | **一致** | ~95% |
| **P4** | 沿冻结 BLA 路径：δ=0 处 E_tot/EN 最低；ΔEN/ΔEe 符号达阈 | BLA 扫描；Ee/EN | 本层级判定已闭合 | **一致** | ~90% |
| **P5** | 已测体系 ΔEAm>0 且 Δr>0（测试域） | GL/GE-m 独立实现 | 本层级判定已闭合 | **一致** | ~96% |
| **P6** | 无参数 ESE：苯 ≈ −36.3；环丁二烯 ≈ +53–55 | 定域几何 + 可加性 + VR | 本层级判定已闭合 | **一致** | ~98% |
| **P7** | C₁₂H₆：PLG 后 BLA 显著坍塌（单体系） | G vs PLG 的 Δr | 本层级已闭合 | **一致** | ~90% |
| **P8** | 呋喃类用 LDE，介于芳香/非芳香 | ΔEA / ΣΔEAm 符号模式 | 本层级已闭合 | **一致** | ~96% |
| **P9** | N=8–18：4n/4n+2 符号正确，gap 收敛趋势（2007 proxy） | VDE 随 N 变化 | 本层级判定已闭合 | **一致** | ~94% |

---

## 三、分阶段路线图（结题状态 · 2026-08-25）

> 下列勾选反映 **本支路已冻结** 的事实。未完成项不阻断结题。

### Phase 0 — 基础设施
- [x] PDF 文本导出与知识库（**仅本地** `data/`；不入库）
- [x] 命题清单 `docs/propositions.md`（P1–P9）
- [x] 本方案更新；环境只读核查（PySCF 已就绪，**不改配置**）
- [x] 建立 `src/` 公共库与九命题包（仅用现有 PySCF）
- [x] `results/<pid>/` 约定 + Canvas 总看板（Canvas 本机托管）

### Phase 1 — 形式化
- [ ] 独立 `docs/formalism.md`（**不阻断结题**；公式散见各 `report.md` 与 `src/localization/`）
- [x] 各命题专家异议预案（见各 `report.md` / 计算与论证表）

### Phase 2 — 引擎能力验证
- [x] 在 WSL/现有 PySCF 上完成苯、丁二烯、环丁二烯等体系计算
- [x] 验证可改写 S / Fock（2007-GL 等已用于 P1/P5/P6/P8/P9）

### Phase 3 — 逐命题验证（顺序已完成）
1. P1 → 2. P3 → 3. P4 → 4. P6+P5 → 5. P7 → 6. P8、P9 → 7. P2（总判）

每条完成标准（**九命题均已满足**）：
- 可复现脚本（`src/`）+ `results/Pn/` 数据
- **质控闸 G1–G5 通过**；作废数据在 `invalid_*/`
- 单元正式交付：`VERDICT.md` + `report.md` + `evidence/`
- Cursor Canvas 专题 + `canvas_link.md`（无本机绝对路径）
- 回填：`implementation/`、命题总表、`docs/propositions.md`

### Phase 4 — 专家异议回应包装
- [x] 各命题 report / VERDICT 含异议回应
- [ ] 与 Shaik/Hiberty、BLW、NICS 等主流的**专章对照文**（可选深化，不改已定判定）

### Phase 5 — 总交付
- [x] `deliverables/final/`：摘要、判定总表、总报告、证据包、检查清单（**已冻结**）
- [x] 总 Canvas：`propositions-board`（本机）

---

## 四、目录结构与交付落盘

权威说明见 [`directory_structure.md`](directory_structure.md)（**GitHub 公开口径**）。摘要：

```
AI4S-OrgChem/
├─ README · NOTICE · LICENSE · CITATION.cff
├─ docs/                      # 命题、门禁、方案、结构说明
├─ implementation/            # 预注册计划表与跟踪
├─ src/                       # 自研代码
├─ results/Pn/tables/         # 入库；raw/ 不入库
├─ deliverables/{unit,final,papers}/
├─ tools/
├─ source/README.md           # 仅说明；PDF 不入库
└─ data/README.md             # 仅说明；抽文本不入库
```

| 交付层级 | 路径 | 核心文件 |
|----------|------|----------|
| 单元计算 | `results/Pn/` | tables（含 `quality_gate`）、invalid_* |
| 单元正式 | `deliverables/unit/Pn/` | **`VERDICT.md`**、`report.md`、`evidence/` |
| 最终 | `deliverables/final/` | `EXECUTIVE_SUMMARY.md`、`VERDICT_TABLE.md`、`FULL_REPORT.md`、`checklist.md`、`evidence_pack/` |

Canvas：Cursor 本机托管（不入库）；链接写入各 `deliverables/unit/Pn/canvas_link.md`。
