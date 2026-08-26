# P4 技术报告：苯 D₆ₕ 与核排斥 / BLA

## 1. 命题

| | |
|--|--|
| 原著主张 | 沿 dra=−drb 等键长化时核排斥最小导致 D₆ₕ；电子能倾向键长交替；GL→G 叙述 ΔEe&gt;0、ΔEN&lt;0。 |
| 传统对立 | π 离域驱动等键长与芳香稳定。 |
| 本项目判据 | E 与 EN 最低在 δ≈0；δ_max→0 时 ΔEN&lt;0 且 ΔEe&gt;0 → 一致。 |

## 2. 方法

- B3LYP/6-31G*（WSL PySCF）
- D₆ₕ：SciPy BFGS + 解析梯度，平面约束
- BLA：δ = r_a − r_b，r_a = r0+δ/2，r_b = r0−δ/2（0–0.12 Å，13 点），单点
- Ee = E − EN；`check_topology` 每点强制

```bash
PYTHONPATH=. python3 -m src.p4_benzene.run --basis '6-31g*'
```

## 3. 结果

- r0 = 1.3969 Å；E/EN 最低均在 δ=0；Ee 最低在 δ=0.12  
- δ=0.12→0：ΔE=−7.75，ΔEN=−96.0，ΔEe=+88.3 kcal/mol  
- `quality_gate.passed=true`，`agree=true`

## 4. 异议预案

| 异议 | 回应 |
|------|------|
| 未用原著 GL | L1/L2 不依赖其程序；完整 GL→G 另计。 |
| 与 Shaik 相同？ | 只验原著 Ee/EN + dra=−drb 叙述。 |
| ΔEN 近百 kcal | 大形变下 EN 变化本就大；判定看符号与极值，ΔE 跨度仅 ~8 kcal。 |

## 5. 判定

见 [`VERDICT.md`](VERDICT.md)：**一致**。
