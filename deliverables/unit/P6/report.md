# P6 技术报告：无参数 ESE（v4）

## 1. 命题与对立

| | |
|--|--|
| 原著主张 | 无参数 ESE：苯 ≈ −36.3 kcal/mol（稳定）；环丁二烯共轭/反芳香能 ≈ +53–55（去稳定）。 |
| 传统对立 | 芳香稳定化能依赖经验参数或等键反应选取。 |
| 本项目判据 | 一致：苯 ESE 落在约 −30…−40，环丁二烯 ΔEA 落在约 +45…+70，且符号/基组稳定；非一致：符号翻转或严重偏离窗。 |

## 2. 方法

```bash
PYTHONPATH=. python3 -m src.p6_ese.run_v2 --basis '6-31g*' --basis2 '6-31g'
PYTHONPATH=. python3 -m src.p6_ese.run_v3 --basis '6-31g*'
PYTHONPATH=. python3 -m src.p6_ese.run_v4_objections --basis '6-31g*'
```

- 定域：2007 Fock+S；主证据 B3LYP
- 加深：半绝热 BLA；CBD 2D 绝热；RHF/2011-lite 符号检验
- **定判以 v4 为准**；`p6_v2`/`p6_v3` 保留为版本链

## 3. 结果摘要（kcal/mol）

| 观测量 | 本项目 | 原著/窗 |
|--------|--------|---------|
| 苯 ESE（垂直 6-31G*） | **−35.44** | ≈−36.3 |
| 苯 ESE（半绝热 BLA） | **−40.61** | 同号、量级合理 |
| CBD vert@G* ΔEA | **+53.98** | ≈+53.6 |
| CBD 2D 绝热 ΔEA | **+65.51** | 窗 45–70 |
| CBD 于 G* 的 ESE | **≈0** | 两双键 ⇒ G≡GE |

## 4. 三条异议（v4 已闭合）

| # | 异议 | 状态 |
|---|------|------|
| O1 | 仅 1D 半绝热 | CBD 2D 非贴边，ΔEA=+65.51 → **闭合** |
| O2 | RHF-2007 爆炸 | 病理入 `invalid_rhf2007_benzene/`；B3LYP/2011-lite 符号 OK → **闭合** |
| O3 | 「CBD 应为 ESE=+53」 | 正确口径是 ΔEA；vert@G*=+53.98 → **闭合** |

## 5. 判定

见 [`VERDICT.md`](VERDICT.md)：**一致**（~98%）。

## 6. 证据

- `results/P6/tables/p6_v4_objections_B3LYP_6-31gs.json` · `summary_p6_v4.md`
- `src/p6_ese/run_v4_objections.py` · `src/localization/gl_2007.py`
