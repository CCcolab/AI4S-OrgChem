# src/ — 自研可复现代码

本目录存放**独立实现**的计算与分析脚本，不复制、不翻译原著程序。

## 建议布局

```
src/
├── common/           # 单位、I/O、Ee/EN、几何
├── p1_butadiene/     # P1 入口与反应能
├── p3_nba/           # …
├── localization/     # GL/GE-m/ESE/PLG/LFMO（P5+ 共用）
└── README.md
```

## 约定

- 每个命题包提供 `run.py` 或 `README.md` 写明：WSL 启动命令、写入哪个 `results/Pn/`。
- 默认：`wsl -e bash -lc 'cd /mnt/d/AI4S-OrgChem && python3 -m src.p1_butadiene.run'`
