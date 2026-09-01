# 最终交付区（已冻结）

[English](README.en.md) | **中文**

**冻结日期**：2026-08-25 · 前置条件（P1–P9 单元 `VERDICT.md` 齐套 + G1–G5 全过）已满足。  
**冻结判定零歧义**：[`docs/FROZEN_VERDICT_AUTHORITY.md`](../../docs/FROZEN_VERDICT_AUTHORITY.md) · [English](../../docs/FROZEN_VERDICT_AUTHORITY.en.md)
**原著作者**：**虞忠衡教授**（Professor Zhongheng Yu）  
**项目性质**：**Cursor-assisted AI4S** 证据支路（本仓库）。  
**第三方复审**：2026-09-01 · 本地报告 `AI4S_Post_Remediation_Review_2026-09-01.pdf`（§4 AI4S-OrgChem；不入库）

| 文件 | 内容 | 基调 |
|------|------|------|
| [`EXECUTIVE_SUMMARY.md`](EXECUTIVE_SUMMARY.md) · [`EXECUTIVE_SUMMARY.en.md`](EXECUTIVE_SUMMARY.en.md) | 一页结论：**L0 五态矩阵** + estimand 分层 | 面向读者的成果综述 |
| [`../../docs/expert_quick_review_guide.md`](../../docs/expert_quick_review_guide.md) | **量子化学专家快速审阅指南（中英）** | 专家 15–30 分钟速读 |
| [`VERDICT_TABLE.md`](VERDICT_TABLE.md) · [`VERDICT_TABLE.en.md`](VERDICT_TABLE.en.md) | 逐命题预注册判据 / 实测值 / 达阈 / 质控闸 | **纯审计**，无理论评价 |
| [`FULL_REPORT.md`](FULL_REPORT.md) · [`FULL_REPORT.en.md`](FULL_REPORT.en.md) | 体裁定位、方法、结果综述、专家异议回应、限定与剩余、复现命令 | **独立复现研究** + AI4S Agent |
| [`checklist.md`](checklist.md) · [`checklist.en.md`](checklist.en.md) | 结题核查：判定齐套 / G1–G5 / 作废声明 / 合规 / 文档同步 | 纯审计 |
| [`evidence_pack/`](evidence_pack/) | 跨命题主证据（每命题一份） | 数据 |
| [`README.en.md`](README.en.md) | English index for this folder | — |

## 引用口径

**对外科学状态（L0）为主键**：各 `deliverables/unit/Pn/VERDICT.md` 的 **对外状态（L0）** 与本目录 [`VERDICT_TABLE.md`](VERDICT_TABLE.md) §1。  
**L1 预注册阈值审计（2026-08-25 快照）**：历史记录，**不是**科学总分。改判/复审响应见 [`FROZEN_VERDICT_AUTHORITY.md`](../../docs/FROZEN_VERDICT_AUTHORITY.md)。

「一致 / 非一致」是**可复现性**判定：一致意味着在原著自身定义与参考态下结果可被第三方独立复现，**不**意味着该能量分解是唯一正确的物理图像，**也不**意味着传统共轭/芳香性理论被推翻。详见 [`FULL_REPORT.md`](FULL_REPORT.md) §1.1。

## 合规

未使用、复制、反编译或翻译原著任何程序代码；主证据引擎仅为 WSL 中既有 PySCF 2.14.0，未改动环境配置。详见 [`FULL_REPORT.md`](FULL_REPORT.md) §6。
