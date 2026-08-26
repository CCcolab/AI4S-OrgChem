# deliverables/ — 正式交付区

本目录为采用最新 AI 技术的 **AI for Science（AI4S）Agent** 项目的正式科学产出（判定书、报告、证据），不是聊天记录。

## unit/ — 单元交付（每命题）

`unit/Pn/` 齐套后方可将该命题标为完成：

| 文件 | 含义 |
|------|------|
| `VERDICT.md` | **一致** 或 **非一致**（唯一正式判定；须先过 [`docs/quality_gates.md`](../docs/quality_gates.md)） |
| `report.md` | 方法、数据、敏感性、异议回应 |
| `evidence/` | 精选表图 |
| `canvas_link.md` | Cursor Canvas 路径 |

**禁止**：未收敛、错误分子、或 `results/Pn/invalid_*/` 中的数据支撑终裁；证据不足时保持待定。

索引：见 [unit/INDEX.md](unit/INDEX.md)。

## final/ — 最终交付（**已于 2026-08-25 冻结**）

九命题单元判定齐套（9/9，G1–G5 全过）后冻结：

| 文件 | 含义 | 基调 |
|------|------|------|
| `EXECUTIVE_SUMMARY.md` | 一页总览：8 一致 + 1 非一致 | 成果综述 |
| `VERDICT_TABLE.md` | 逐命题判据 / 实测 / 达阈 / 门禁 | 纯审计 |
| `FULL_REPORT.md` | 方法、结果综述、异议回应、限定、复现命令 | 独立复现研究 |
| `checklist.md` | 判定齐套 / 门禁 / 作废声明 / 合规 / 文档同步 | 纯审计 |
| `evidence_pack/` | 跨命题主证据（每命题一份） | 数据 |

## papers/ — 学术产出骨架

| 文件 | 内容 |
|------|------|
| `chem_preprint_outline.md` | 化学篇：八条确证与一条异议（含撰稿纪律） |
| `ai4s_method_outline.md` | AI4S 方法篇：证据门禁与 P3 判定反转案例 |

详细约定见 [`docs/directory_structure.md`](../docs/directory_structure.md)。
