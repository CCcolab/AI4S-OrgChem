# 总证据包

跨命题关键数据，每条命题一份主证据文件。原始出处见各文件对应的 `results/Pn/tables/`。

| 文件 | 命题 | 内容 | 原始路径 |
|------|------|------|----------|
| `P1_ce_summary.csv` | P1 | CE₁/CE₂ 四层次汇总（RHF/B3LYP/MP2/+ZPE） | `results/P1/tables/ce_summary.csv` |
| `P1_gl2007.json` | P1 | 自研 2007-GL（Fock+S）丁二烯 ΔE=+4.06 | `results/P1/tables/gl2007_butadiene_B3LYP_6-31gs.json` |
| `P2_lfmo_lite.json` | P2 | LFMO-lite EV / Enσσ / Eπσ 扭转扫描 | `results/P2/tables/p2_v2_lfmo_lite.json` |
| `P3_pes_tight.csv` | P3 | 加严双向 θ 约束 PES（θ_min=44.9°） | `results/P3/tables/pes_tight_*.csv` |
| `P4_bla_scan.json` | P4 | 苯 BLA 扫描与 Ee/EN 分解 | `results/P4/tables/bla_scan_B3LYP_6-31gs.json` |
| `P5_v5c.json` | P5 | ΔEAm 与 Δr（v5c 三异议闭合） | `results/P5/tables/p5_v5c_objections_*.json` |
| `P6_v4.json` | P6 | 苯 ESE 与环丁二烯 ΔEA（v4 三异议闭合） | `results/P6/tables/p6_v4_objections_*.json` |
| `P7_v1c.json` | P7 | Δr(G) → Δr(PLG) 坍塌 | `results/P7/tables/p7_v1c_B3LYP_6-31gs.json` |
| `P8_v2.json` | P8 | 呋喃/吡咯/噁唑 LDE 与苯 ESE 对照 | `results/P8/tables/p8_v2_B3LYP_6-31gs.json` |
| `P9_v2b.json` | P9 | VDE/π 序列与 gap 收敛（含 2011-lite、BLA） | `results/P9/tables/p9_v2b.json` |

## 读法

- 每个 JSON 含 `quality_gate` 字段；`passed` 必须为 `true` 才被判定引用。
- 能量单位：JSON 内部为 Hartree，`*_kcal` 字段为 kcal/mol。稳定化为负、去稳定为正。
- 拒收数据**不在**本包内（清单见 [`../checklist.md`](../checklist.md) §3）。
