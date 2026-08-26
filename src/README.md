# src/ — 自研可复现代码

本目录存放 **AI4S Agent 独立实现** 的计算与分析脚本，不复制、不翻译原著程序。

## 布局

```
src/
├── common/              # 单位、I/O、Ee/EN、几何
├── localization/        # GL / GE-m / ESE / PLG / LFMO-lite / hetero_gl
├── p1_butadiene/        # 定判入口：run.py · run_zpe.py · run_gl.py
├── p2_aggregate/        # summarize.py · run_lfmo_lite.py
├── p3_nba/              # 定判入口：run_tight.py（勿用早期欠收敛协议定判）
├── p4_benzene/          # run.py（BLA）
├── p5_local/            # 定判入口：run_v5c_objections.py
├── p6_ese/              # 定判入口：run_v4_objections.py
├── p7_strained/         # 定判入口：run_v1c.py
├── p8_furan/            # 定判入口：run_v2.py
├── p9_annulene/         # 定判入口：run_v2b.py
└── README.md
```

更早的 `run_v1` / `run_pilot` 等保留为可复现版本链，**定判入口以上表为准**。

## 约定

- WSL 示例：`cd /mnt/d/AI4S-OrgChem && PYTHONPATH=. python3 -m src.p6_ese.run_v4_objections --basis '6-31g*'`
- **质控（强制）**：生成/优化几何须有拓扑断言；结果 JSON 含 `quality_gate`，未通过不得暗示可写 VERDICT。见 [`docs/quality_gates.md`](../docs/quality_gates.md)。
- 作废数据写入 `results/Pn/invalid_<原因>/`，不得留在主 `tables/` 冒充正式证据。
- 全套复现命令：[deliverables/final/FULL_REPORT.md](../deliverables/final/FULL_REPORT.md) §7。
