# GitHub 上传实施方案

**文档日期**：2026-08-26  
**状态**：准备就绪，**尚未** `git init` / 建仓 / 推送（须你确认 Owner 与可见性后执行）  
**项目性质**：采用最新 AI 技术的 **AI for Science（AI4S）Agent**  
**本支路**：Cursor · Grok 独立验证（三大 LLM Agent 计划之一）  
**结题冻结**：2026-08-25 · **8 一致 + 1 非一致（P1）**

权威目录树见 [`directory_structure.md`](directory_structure.md)。计算平面见 [`../WSL2/README.md`](../WSL2/README.md)。

---

## 0. 目标与成功标准

| 目标 | 成功标准 |
|------|----------|
| 公开可审计的复现制品 | 他人克隆后能读懂背景/价值/成就，并按 `FULL_REPORT` §7 在自有 WSL+PySCF 上对照 |
| 零版权泄露 | 仓内无原著 PDF、无全书抽文本 |
| 叙事清晰 | README 标明 AI4S Agent、三模型计划、Cursor Grok 支路、WSL2 计算 |
| 可引用 | `LICENSE` + `CITATION.cff`（真实 URL）+ Release tag |

---

## 1. 已拍板决策（勿临时改口径）

| 项 | 决定 |
|----|------|
| 仓库名 | `AI4S-OrgChem`（可改，需同步 `CITATION.cff`） |
| 可见性 | **待你确认**：Public，或先 Private 审阅再转 Public |
| License | MIT（`LICENSE`）+ `NOTICE` |
| `source/` PDF | **不上传**（仅 `source/README.md`） |
| `data/` 抽文本 | **不上传**（仅 `data/README.md`） |
| `.cursor/` / Canvas | **不上传** |
| `results/**/raw/` | **不上传** |
| `WSL2/` | **上传**（inventory + 小体积 home 拷贝；不含 8.9G micromamba） |
| 多模型 | 本仓 = Cursor Grok 一支；B/C 另仓 |

---

## 2. 公开目录树（应推送内容）

```
AI4S-OrgChem/
├── README.md · README.zh-CN.md · NOTICE · LICENSE · CITATION.cff · .gitignore
├── docs/                 # 含本文件、propositions、quality_gates、directory_structure…
├── implementation/       # 预注册计划表
├── src/                  # 自研九命题包 + localization
├── results/Pn/           # tables + invalid_*（无 raw）
├── deliverables/         # unit / final / papers
├── WSL2/                 # 计算平面归档 + 详细计算介绍
├── tools/                # check_wsl_env · wsl_snapshot · …
├── source/README.md
└── data/README.md
```

**禁止出现在首次 commit 中**：`*.PDF`、`data/book_*`、`.cursor/`、`.specstory/`、`*.canvas.tsx`、`__pycache__/`、`results/**/raw/**`、多 GB env。

---

## 3. 分期实施

### 阶段 A — 上传前终检（约 15–30 分钟，人工）

| # | 动作 | 判据 |
|---|------|------|
| A1 | 填 `CITATION.cff` 的 `REPLACE_OWNER` | 变为真实 `https://github.com/<owner>/AI4S-OrgChem` |
| A2 | 选定可见性 | Public / Private 二选一 |
| A3 | 确认 `gh auth status` 已登录目标账号 | 能 `gh repo create` |
| A4 | 版权抽检 | 见 §4 清单全部勾选 |
| A5 | 大文件抽检 | 拟提交树无 >5 MB 单文件（PDF 必须被 ignore） |
| A6 | 口径抽检 | README 含 AI4S Agent、8+1、WSL2 链接；P3 evidence 仅为 `pes_tight_*` |

可选：在 WSL 刷新一次环境指纹（不改软件）：

```bash
wsl -d Ubuntu-24.04 -e bash /mnt/d/AI4S-OrgChem/tools/wsl_snapshot.sh
```

### 阶段 B — 本地 Git 初始化（约 5 分钟）

在 **PowerShell**（或 WSL）于仓库根执行：

```powershell
cd D:\AI4S-OrgChem

git init -b main
git add -A
git status
```

**硬停**：目视 `git status`：

- 必须看到：`README.md`、`src/`、`deliverables/`、`WSL2/`、`implementation/`…
- 不得看到：任何 `.PDF`、`book_full.txt`、`book_text/`、`.cursor/`

辅助核对（PowerShell）：

```powershell
# 应均返回路径（表示被 ignore）
git check-ignore -v source\*.PDF 2>$null
git check-ignore -v data\book_full.txt

# 暂存区不得含 PDF
git diff --cached --name-only | Select-String -Pattern '\.PDF$|book_full|book_text|\.cursor'
# 无输出 = 通过
```

通过后再提交：

```powershell
git commit -m @"
Initial public release: AI4S Agent (Cursor Grok) P1-P9 freeze.

Eight agree, one dissent (P1). No copyrighted book PDF or extracted text.
Includes WSL2 compute-plane archive and deliverables/final pack.
"@
```

### 阶段 C — 创建远程并推送（约 5 分钟）

**方案 C1 — `gh` 一键（推荐）**

```powershell
# Public：
gh repo create AI4S-OrgChem --public --source=. --remote=origin --push

# 或先 Private：
# gh repo create AI4S-OrgChem --private --source=. --remote=origin --push
```

**方案 C2 — 网页建空仓再推**

1. GitHub 新建空仓库（不要勾选自动加 README）  
2. 本地：

```powershell
git remote add origin https://github.com/<OWNER>/AI4S-OrgChem.git
git push -u origin main
```

### 阶段 D — 打标签与 Release（约 10 分钟）

```powershell
git tag -a v1.0.0-freeze-2026-08-25 -m "P1-P9 frozen: 8 agree + 1 dissent (P1)"
git push origin v1.0.0-freeze-2026-08-25

gh release create v1.0.0-freeze-2026-08-25 `
  --title "v1.0.0 — P1–P9 freeze (Cursor Grok arm)" `
  --notes @"
## Summary
- AI4S Agent independent replication of Yu's nine propositions
- Result: **8 agree**, **1 dissent (P1)**
- Compute plane: WSL2 Ubuntu 24.04 + PySCF 2.14.0 (GPU-visible)
- No copyrighted book PDF or full-text extracts in this repository

## Start here
- README.md
- deliverables/final/EXECUTIVE_SUMMARY.md
- WSL2/README.md
"@
```

### 阶段 E — 仓库抛光（建仓后当天）

| # | 动作 |
|---|------|
| E1 | About：`Independent PySCF / AI4S Agent replication (8 agree / 1 dissent)` |
| E2 | Topics：`ai-for-science` `ai4s` `reproducibility` `pyscf` `quantum-chemistry` `organic-chemistry` |
| E3 | 默认分支保护（可选）：禁止 force push `main` |
| E4 | 确认 GitHub 文件浏览器中无 PDF、有 `WSL2/inventory/` |
| E5 | 把真实 URL 写回本地 `CITATION.cff`（若阶段 A 已填则跳过）并再推一小 commit |

### 阶段 F — 后续（非阻塞首发）

| 项 | 说明 |
|----|------|
| 支路 B/C | 另仓独立实现；齐套后可加 `deliverables/final/CROSS_LLM_TABLE.md` |
| 预印本 | 按 `deliverables/papers/` 骨架扩写；出 arXiv 后回填 README |
| Zenodo | 可选：Release 同步 DOI |
| CI | 可选轻量：禁止提交 `*.PDF`；**不要**在 CI 全量跑 SCF |

---

## 4. 上传前检查清单（打印勾选）

```text
□ 仓库尚未含 .git，或已确认将在干净树上 init
□ CITATION.cff 中 REPLACE_OWNER 已替换
□ 可见性已选定（Public / Private）
□ gh / git 凭据可用
□ git check-ignore 挡住 source PDF 与 data/book_full.txt
□ git status / staged 无 .PDF、无 book_text、无 .cursor
□ README：AI4S Agent + 背景/价值/成就 + WSL2 链接 + 8+1
□ NOTICE：不提供原著材料
□ deliverables/unit/P3/evidence 仅为 pes_tight_*（非作废几何）
□ canvas_link.md 无 C:\Users\… 绝对路径
□ WSL2/README.md 含详细计算介绍；inventory 在
□ 拟提交无 >5 MB 可疑大文件
□ 承诺：未获确认前不 push；push 后不 force 改写已公开历史以藏 PDF
```

---

## 5. 风险与应急

| 风险 | 预案 |
|------|------|
| 误把 PDF 推进去 | **立刻**使仓库 Private → 从历史删除并轮换（公开后仅改最新 commit 不够）；预防靠 A4/A5 |
| `gh` 未登录 | `gh auth login` 后重试 C |
| 路径含中文文件名 | Git for Windows 默认可处理；若乱码，统一 UTF-8，勿改文件名除非必要 |
| WSL 符号链接 | `WSL2/home_cascade/ai4orgchem_env` 已改为说明文件，避免 GitHub 断链 |
| 审稿/合作者未齐 | 用 **Private** 首发，再改 Public |

---

## 6. 角色分工建议

| 角色 | 职责 |
|------|------|
| 你（Owner） | 定 Public/Private、填 GitHub 用户名、执行或授权 `gh repo create` |
| Agent（本对话） | 已完成内容清洗与文档；**仅在你明确说「初始化并推送」后**执行 B–D |
| 审阅者（可选） | Private 阶段读 `EXECUTIVE_SUMMARY` + `checklist` |

---

## 7. 当前缺口（阻塞项仅一项）

| 缺口 | 是否阻塞首发 |
|------|----------------|
| `CITATION.cff` 的 `REPLACE_OWNER` | **是**（至少推送前改掉） |
| 可见性未口头确认 | **是** |
| `docs/formalism.md` 未建 | 否 |
| 支路 B/C 未启动 | 否 |
| GitHub Actions | 否 |

---

## 8. 一句话执行令（给你复制）

当你准备好时，回复例如：

> Owner=`<你的GitHub用户名>`，可见性=`Public`（或 Private），请初始化并推送。

即可按本方案阶段 B→E 执行。在此之前**不**运行 `git init` / `gh repo create` / `git push`。
