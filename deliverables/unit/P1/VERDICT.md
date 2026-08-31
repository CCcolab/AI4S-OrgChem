# P1 判定（VERDICT）

[English](VERDICT.en.md) | **中文**

> **冻结判定（L1 唯一权威）**：**非一致** · 冻结日 **2026-08-25**。  
> 下文「子估计量」「阅读辅助」**不改变**本判定。见 [`docs/FROZEN_VERDICT_AUTHORITY.md`](../../../docs/FROZEN_VERDICT_AUTHORITY.md)。

- **判定：非一致**
- 日期：2026-08-23（ZPE + L3 自研 GL 补强后维持 L1/L2 判定）
- 状态说明：L1/L2 主判据 + RHF/MP2/ZPE 敏感性已闭合；`quality_gate` 已落盘；L3 自研 2007-GL 已测并对「+1.4」给出**子结论**（不回溯改写 L1/L2）。完成度 **~98%**（2011/2014 交换积分删除未实现，属 P5/P6 深化）。
- 方法：RHF / B3LYP / MP2（MP2 在 RHF 优化几何上单点），基组 6-31G*；SciPy-BFGS + PySCF 解析梯度（未装 geometric/berny，未改环境）
- 定义（与教材氢化热差一致，稳定化为负）：

```
ΔH_hyd(M) = E(n-butane) − E(M) − n_H2 · E(H2)
CE_ref    = ΔH_hyd(butadiene) − 2 · ΔH_hyd(ref)
```

## 子估计量（阅读辅助；不改 L1 非一致）

| ID | 估计量 | 阅读辅助 | 与 L1 关系 |
|----|--------|----------|------------|
| P1a | 实验 ΔH₂₉₈ 氢化热差 | `INCOMPARABLE` | 未作主交付；不得与 P1b 混读 |
| P1b | E_ele / E₀K（B3LYP 等）CE₁/CE₂ | `OPPOSED` | 支撑 L1 **非一致**（无符号翻转） |
| P1c | 2007-GL +4.06 vs 书 +1.4（2014） | `INCOMPARABLE` | 定义未等价；**不回溯**改写 P1b/L1 |

## 判据对照（主交付：B3LYP）

| 观测量 | 一致阈值 | 本项目结果 | 是否满足 |
|--------|----------|------------|----------|
| CE1（vs 1-丁烯） | < 0（稳定化） | **+9.331 kcal/mol** | 否 |
| CE2（vs trans-2-丁烯） | > 0（去稳定） | **+1.861 kcal/mol** | 是 |
| CE2 量级 | 约 +1–3 kcal/mol | +1.86 | 是 |
| **符号翻转** CE1<0 且 CE2>0 | 必须 | **未出现**（二者皆正） | **否** |

## 方法敏感性（同基组，不改判）

| 层次 | CE1 | CE2 | 符号翻转 |
|------|-----|-----|----------|
| RHF/6-31G* | +8.130 | +2.449 | 否 |
| B3LYP/6-31G* | +9.331 | +1.861 | 否 |
| MP2//RHF/6-31G* | +8.503 | +1.858 | 否 |
| B3LYP + ZPE（E0K） | +8.540 | +0.511 | 否 |

三层电子能及 ZPE 校正后 **CE1、CE2 皆为正**；敏感性**未**复现预设翻转，故维持 **非一致**。

## 结论要点

1. 在本项目预设计算判据下，**不能**判定为「一致」：未复现「换参考分子则共轭能符号翻转」。
2. CE2 在 B3LYP/MP2 下均接近书中 +1.9 kcal/mol，但 CE1 远离经典 −3.9 且与 CE2 同号，不能据此支持完整主张。
3. RHF→DFT→MP2→ZPE 同号结构稳健，削弱「仅因 B3LYP / 缺 ZPE 偶然失败」的疑虑。
4. **补充（不改判定）**：文献实验氢化热算术可出现参考态叙事变化；正式判定仍绑定独立 QC 复现判据。
5. **L3 子结论**（独立于 L1/L2）：按 Ch6 公开 2007 定义自研 GL（Fock+重叠 π–π 块删除），丁二烯 ΔE=E(G)−E(GL)≈**+4.06 kcal/mol**（去稳定同号），**未**落入书中 **+1.4**（Ch12 归因于 2014 法）容差 → 对「≈+1.4」**非一致**；交换积分删除未做。

## 敏感性 / 局限

- 已完成：RHF / B3LYP / MP2 电子能；ZPE（谐波，落盘几何）；L3 2007-GL（Fock+S 主协议 + Fock-only 对照）。
- 未含：2011/2014 交换积分删除；完整笛卡尔 GL 解析梯度优化（本层用 r23 一维扫描）。
- 无 geometric/berny，优化器为 SciPy BFGS。

## 质控闸（对照 `docs/quality_gates.md`，不改判）

| 闸 | 结果 | 证据 |
|----|------|------|
| G1 拓扑/异构体 | **通过** | 五分子化学式正确；丁二烯双–单–双；trans-2-丁烯二面角 180° |
| G2 几何健全 | **通过** | 无异常短接触 |
| G3 收敛 | **通过** | B3LYP \|g\|_max ≤ 9.9e−5；`ce_*.json` 含 `quality_gate.passed=true`（含事后回填） |
| G4 能量尺度 | **通过** | CE ~1–9 kcal；L3 Fock+S ΔE ~4 kcal（Fock-only 诊断为病态对照，不作主证据） |
| G5 路径洁净 | **通过** | 无 `invalid_*/`；CE / ZPE / GL 可独立重算 |

与 P3「未收敛+错误分子」误判**不同类**。

审计：`results/P1/tables/p1_quality_audit.json`

## 证据路径

- 汇总：`results/P1/tables/sensitivity_methods.json`、`ce_summary.csv`
- 分项：`ce_RHF_6-31gs.json`、`ce_B3LYP_6-31gs.json`、`ce_MP2_6-31gs.json`
- ZPE：`ce_zpe_B3LYP_6-31gs.json`
- L3 GL：`gl2007_butadiene_B3LYP_6-31gs.json`
- 报告：`deliverables/unit/P1/report.md`
- 几何：`results/P1/raw/{RHF,B3LYP,MP2}_6-31gs/`
