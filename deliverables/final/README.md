# 最终交付区（已冻结）

**冻结日期**：2026-08-25 · 前置条件（P1–P9 单元 `VERDICT.md` 齐套 + G1–G5 全过）已满足。  
**项目性质**：采用最新 AI 技术的 **AI for Science（AI4S）Agent** 验证支路（Cursor · Grok）。

| 文件 | 内容 | 基调 |
|------|------|------|
| [`EXECUTIVE_SUMMARY.md`](EXECUTIVE_SUMMARY.md) | 一页结论：AI4S Agent 定位 + 8 一致 + 1 非一致 | 面向读者的成果综述 |
| [`VERDICT_TABLE.md`](VERDICT_TABLE.md) | 逐命题预注册判据 / 实测值 / 达阈 / 质控闸 | **纯审计**，无理论评价 |
| [`FULL_REPORT.md`](FULL_REPORT.md) | 体裁定位、方法、结果综述、专家异议回应、限定与剩余、复现命令 | **独立复现研究** + AI4S Agent |
| [`checklist.md`](checklist.md) | 结题核查：判定齐套 / G1–G5 / 作废声明 / 合规 / 文档同步 | 纯审计 |
| [`evidence_pack/`](evidence_pack/) | 跨命题主证据（每命题一份） | 数据 |

## 引用口径

正式判定只以本目录与各 `deliverables/unit/Pn/VERDICT.md` 为准。

「一致 / 非一致」是**可复现性**判定：一致意味着在原著自身定义与参考态下结果可被第三方独立复现，**不**意味着该能量分解是唯一正确的物理图像，**也不**意味着传统共轭/芳香性理论被推翻。详见 [`FULL_REPORT.md`](FULL_REPORT.md) §1.1。

## 合规

未使用、复制、反编译或翻译原著任何程序代码；主证据引擎仅为 WSL 中既有 PySCF 2.14.0，未改动环境配置。详见 [`FULL_REPORT.md`](FULL_REPORT.md) §6。
