# P5 技术报告：局部共轭 ΔEAm（v5c）

## 1. 命题与对立

| | |
|--|--|
| 原著主张 | 局部双键对间共轭能 ΔEAm 恒去稳定（>0），并使中间单键伸长（Δr>0）。 |
| 传统对立 | 局部共轭稳定并使单键缩短。 |
| 本项目判据 | 一致：ΔEAm 成体系为正，且关键多烯 Δr>0；非一致：ΔEAm 系统性为负或 Δr 系统性缩短。 |

## 2. 方法

- 定域：2007 GL/GE-m（片段间 π–π Fock + 重叠块删除）
- 层次：B3LYP/6-31G*
- 体系：丁二烯（接 P1 GL）；己三烯（Ci 绝热 + 非对称桥密扫）；苯（Kekulé 垂直，局部 vs 全局对照）

```bash
PYTHONPATH=. python3 -m src.p5_local.run_v4 --basis '6-31g*'
PYTHONPATH=. python3 -m src.p5_local.run_v5b --basis '6-31g*'
PYTHONPATH=. python3 -m src.p5_local.run_v5c_objections --basis '6-31g*'
```

定判以 **v5c** 为准；更早的 `p5_pilot` / `v4` / `v5b` JSON 保留为版本链，不作终裁。

## 3. 结果摘要（kcal/mol）

| 体系 | 量 | 结果 |
|------|-----|------|
| 丁二烯 | ΔEAm / Δr | **+4.057** / **+0.018 Å** |
| 己三烯 Ci 绝热 | ΔEAm | **+0.616** |
| 己三烯非对称密扫 | ΔEAm / Δr* | **+0.982** / **+0.0045 Å** |
| 苯局部 | ΔEAm×3 | **+9.63×3**（全正） |
| 苯全局 | ΔEA / ESE_proxy | **−6.54** / **−35.44**（→ P6） |

**6/6 ΔEAm>0**；丁二烯与己三烯 Δr>0。

## 4. 三条异议（v5c 已闭合）

| 异议 | 处置 |
|------|------|
| 贴边几何出现负 ΔEAm | 判定为伪 PES；内点/垂直对照为正；贴边数据入 `invalid_*` |
| 己三烯 Δr 未硬证 | 密扫得 Δr*=+0.0045 Å |
| 苯 ΔEA<0 是否推翻 P5 | 局部仍去稳定；全局稳定属 ESE（P6），非 P5 失败 |

## 5. 判定

见 [`VERDICT.md`](VERDICT.md)：**一致**（~96%）。

## 6. 证据

- `results/P5/tables/p5_v5c_objections_B3LYP_6-31gs.json` · `summary_p5_v5c.md`
- `src/p5_local/run_v5c_objections.py`
