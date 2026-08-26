# results/ — 单元计算结果

每个命题一个子目录 `P1` … `P9`：

```
Pn/
├── tables/        # ★ 入库：csv/json（含 quality_gate）
├── invalid_*/     # ★ 入库：作废轨迹（禁支撑 VERDICT）
├── logs/          # 可选运行日志
└── raw/           # 本机日志/检查点 — gitignore，不入库
```

**定判文件**（勿用更早版本冒充终裁）：

| Pn | 主证据 |
|----|--------|
| P1 | `ce_*.json` · `gl2007_*.json` · `sensitivity_methods.json` |
| P2 | `p2_v2_aggregate.json` · `p2_v2_lfmo_lite.json` |
| P3 | `pes_tight_*.csv/json`（**不要**用 `invalid_wrong_geometry/`） |
| P4 | `bla_scan_B3LYP_6-31gs.json` |
| P5 | `p5_v5c_objections_*.json` |
| P6 | `p6_v4_objections_*.json` |
| P7 | `p7_v1c_*.json` |
| P8 | `p8_v2_*.json` |
| P9 | `p9_v2b.json` |

正式判定见 `deliverables/unit/Pn/VERDICT.md`。质控见 [`docs/quality_gates.md`](../docs/quality_gates.md)。
