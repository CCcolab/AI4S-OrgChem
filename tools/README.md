# tools/ — 运维与环境脚本

**入库脚本**：

| 文件 | 用途 |
|------|------|
| `check_wsl_env.sh` | 只读核查 python3 / PySCF / GPU |
| `inventory_env.sh` | 扩展环境盘点 |
| `wsl_snapshot.sh` | 刷新 `WSL2/` 归档（环境指纹 + 小体积 home 拷贝） |
| `run_p1.sh` | WSL 便捷启动 P1（路径可按本机修改） |
| `extract_pdf.py` | **仅本地**从 PDF 抽文本到 `data/`；**输出不得提交** |

一次性调试探针已删除，不进入 GitHub。

| `validate_repo.py` | 冻结 L1 一致性、JSON 解析、关键 evidence 的 `quality_gate` 检查 |
