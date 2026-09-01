# AI4S-OrgChem

[![release](https://img.shields.io/badge/release-v1.0.1-orange)](https://github.com/CCcolab/AI4S-OrgChem/releases/tag/v1.0.1-post-freeze-2026-08-28)
[![verdict](https://img.shields.io/badge/L0-estimand_layered-blue)](docs/FROZEN_VERDICT_AUTHORITY.en.md)
[![quality gates](https://img.shields.io/badge/quality_gates-G1--G5_passed-success)](docs/quality_gates.en.md)
[![license](https://img.shields.io/badge/license-MIT-brightgreen)](LICENSE)
[![evidence](https://img.shields.io/badge/evidence-P1--P9_evidence_pack-blue)](deliverables/unit/INDEX.en.md)
[![WSL2](https://img.shields.io/badge/WSL2-Ubuntu_24.04-blue)](WSL2/README.md)
[![PySCF](https://img.shields.io/badge/PySCF-2.14.0-blue)](WSL2/README.md)
[![AI4S](https://img.shields.io/badge/AI4S-Cursor--assisted_workflow-purple)](docs/expert_quick_review_guide.md)

**English** | [中文](README.zh-CN.md)

> **This repository is only** [`CCcolab/AI4S-OrgChem`](https://github.com/CCcolab/AI4S-OrgChem) (P1–P9 PySCF replication arm).  
> Do **not** confuse with sibling arms under the same owner:  
> [`AI4OrgChem`](https://github.com/CCcolab/AI4OrgChem) · [`AI4S-AI4OrgChem`](https://github.com/CCcolab/AI4S-AI4OrgChem) — separate repos, separate implementations, **not** linked from this codebase.

**This is a Cursor-assisted AI-for-Science (AI4S) workflow**: pre-registered criteria and hard evidence gates drive proposition decomposition, independent coding, QC orchestration, and formal verdicts—not chat commentary on the book. (Until agent runtime traces are public, we use **Cursor-assisted workflow**; see [`docs/FROZEN_VERDICT_AUTHORITY.en.md`](docs/FROZEN_VERDICT_AUTHORITY.en.md).)

**Scientific task**: Independent third-party replication of nine core propositions (P1–P9) from **Professor Zhongheng Yu**’s **Questioning Fundamental Principles of Organic Chemistry** that challenge traditional organic structural theory. Each proposition receives an **Agree / Disagree** verdict. The primary computation engine is **PySCF** in a fixed existing environment.

> **This repository does not redistribute the book PDF or full-text extract.** Please obtain the monograph yourself through lawful channels. Local workspace folders `source/` and `data/` are **never** uploaded to GitHub.

---

## Background

Conjugation stabilization and aromaticity energy criteria are among the bedrock tenets of organic structural theory. **Professor Zhongheng Yu**'s monograph systematically questions several of these "fundamental principles"—e.g., local/pairwise conjugation being destabilizing relative to localized references, geometric distortions attributable to π rather than strain, large annulenes approaching polyene character, and furan-type rings requiring distinct criteria from benzene-like ESE. These propositions rely heavily on author-developed localized energy decompositions, **previously lacking independent third-party testing implemented from public definitions alone**.

This project operationalized those challenges into nine testable proposition units (P1–P9): locking criteria and thresholds before execution (pre-registration), having an **Agent independently write** PySCF scripts using only public mathematical definitions, orbital diagrams, and energy formulas from the PDF, and enforcing hard evidence quality gates (topology / geometry / convergence / scale / clean path) before recording formal verdicts.

**Genre**: An independent replication/verification study embedded in an AI4S Agent workflow—neither a polemic nor promotional material. "Agree" indicates replicability under the monograph's own definitions; **it does not imply** the energy decomposition is the sole physically correct picture, **nor does it imply** textbook organic chemistry is overturned.

---

## What this is: an AI4S Agent, not “AI-written reports”

| Dimension | This project |
|-----------|--------------|
| **Paradigm** | **AI for Science (AI4S)**: frontier LLM agents drive an auditable research workflow |
| **Execution** | **Agent loop**: read specs → write/edit code → run PySCF → pass G1–G5 gates → write `VERDICT`; on failure, isolate voided data and rerun |
| **Stack** | Latest-generation LLM Agent + electronic-structure compute (PySCF) + project-level rules/gates (not one-off prompts) |
| **Output** | Reproducible scripts, tables with `quality_gate`, unit verdicts, freeze pack—all runnable and auditable by third parties |
| **Deliberately omitted** | No ML potentials/symbolic regression in place of primary evidence; no chat statements treated as formal verdicts; no copying author code |

Scientific authority derives from **gated evidence chains**; the value of the AI4S Agent is executing, logging, and exposing this chain.

---

## Research value

1. **Auditable research paradigm for AI4S Agents**  
   Demonstrates how modern AI can execute a **genuine research loop** (code—compute—gate—verdict) rather than unverified narratives. If gates fail, the status remains "undetermined"—providing a valid exit path for insufficient evidence and curbing confident errors.

2. **Replicability before correctness judgments**  
   First answering whether claims can be reproduced under public definitions by a third party before debating theoretical merits.

3. **Pre-registered criteria + hard evidence gates**  
   Criteria are locked before runs; G1–G5 failures block verdicts. P3 nearly received a wrong **Disagree** from a misbuilt molecule and under-converged PES; gates stopped it, voided data were isolated, and a rerun flipped to **Agree**—verdicts follow the evidence process, not model verbal preference.

4. **Independent implementation from public definitions only**  
   No contact with or copying of the author’s programs; testing definition-level replicability is stronger than “same code, same number.”

5. **Restrained academic framing**  
   Reference-state dependence of energy decomposition, explicit caveats, and one formal **Disagree** (P1) are stated upfront. A replication with dissent is more credible than unanimous agreement.

---

## Final verdict — replication of Professor Yu’s P1–P9

> **Public scientific status (L0) is primary**: see [`docs/FROZEN_VERDICT_AUTHORITY.en.md`](docs/FROZEN_VERDICT_AUTHORITY.en.md) §2. **Do not** use “8 Agree · 1 Disagree” as a scientific score.  
> **L1 pre-registered threshold audit (2026-08-25 snapshot)** is historical only—see [`deliverables/final/VERDICT_TABLE.en.md`](deliverables/final/VERDICT_TABLE.en.md).

**Repository role (post 2026-09-01 review)**: P1–P9 **Cursor-assisted evidence arm**; results layered by **estimand and scope**. **8 independent propositions + P2 derived index** (P2 excluded from independent verification count).

| L0 status | Propositions |
|-----------|--------------|
| `SUPPORTED_WITHIN_SCOPE` | P3, P5, P6, P8 |
| `PARTIAL` | P1 (composite estimands), P4, P7, P9 |
| `DERIVED` | P2 (meta aggregate) |

**Read the conclusions (start here):**

| Doc | Link |
|-----|------|
| **Verdict table (audit)** | [deliverables/final/VERDICT_TABLE.en.md](deliverables/final/VERDICT_TABLE.en.md) · [中文](deliverables/final/VERDICT_TABLE.md) |
| **Frozen verdict authority** | [**docs/FROZEN_VERDICT_AUTHORITY.en.md**](docs/FROZEN_VERDICT_AUTHORITY.en.md) · [中文](docs/FROZEN_VERDICT_AUTHORITY.md) |
| **Executive summary** | [deliverables/final/EXECUTIVE_SUMMARY.en.md](deliverables/final/EXECUTIVE_SUMMARY.en.md) · [中文](deliverables/final/EXECUTIVE_SUMMARY.md) |
| **Unit verdicts P1–P9** | [deliverables/unit/INDEX.en.md](deliverables/unit/INDEX.en.md) · [中文](deliverables/unit/INDEX.md) |

> **Agree / Disagree = replicability under Yu’s definitions**, not “the book is right/wrong” or “textbook theory is overturned.” See the executive summary for caveats.

### Three strongest positive lines of evidence

| # | Prop. | Key result |
|---|-------|------------|
| 1 | **P6** | Parameter-free ESE: benzene **−35.44** kcal/mol (Yu ≈ −36.3); cyclobutadiene vertical ΔEA **+53.98** (Yu ≈ +53.6)—**both opposite-sign benchmarks hit** |
| 2 | **P7** | Strained aromatic: after cutting center–periphery π coupling, central-ring BLA **+0.207 Å → +0.020 Å**, contrary to Mills–Nixon angle-strain dominance |
| 3 | **P3** | NBA relaxed PES minimum **θ ≈ 44.9°**; unconstrained optimization independently **≈34.5°**—challenges the common inference that conjugation stabilization makes coplanarity most stable |

### Methodological assets (AI4S Agent pipeline)

- **Cursor-assisted workflow**: proposition specs → homemade localization (`src/localization/`) → batch compute → gates → L0/L1 verdicts, fully traceable
- Nine unit packs: `VERDICT.md` + `report.md` + `evidence/` + quality-gate fields
- Voided-data isolation: `results/P3/invalid_wrong_geometry/` (auditable; must not support final verdicts)—a template for when an agent is wrong but fluent
- Paper skeletons: [deliverables/papers/](deliverables/papers/) (chemistry preprint + **AI4S methods**: evidence gates and verdict reversal)

### Caveats to keep in view

- Verdicts are **replicability** judgments, not claims of unique physical truth.
- Energy decompositions are **reference-state dependent**; protocols are logged with sensitivities.
- Yu himself uses split definitions: local conjugation destabilization (ΔEAm>0) and benzene extra delocalization stabilization (ESE<0) must not be collapsed into “all conjugation destabilizes.”
- This repo did **not** use, copy, or translate any of the author’s program code; primary evidence is **PySCF** in the existing environment only.

---

## What is / is not in this repository

| Included (public) | Not included (local or obtain yourself) |
|-------------------|-------------------------------------------|
| `src/` homemade scripts | **`source/` book PDF** |
| `docs/`, `implementation/` specs & pre-registration | **`data/` full-text extract** |
| `results/Pn/tables/` and voided-zone notes | Author appendix / program code |
| `deliverables/` verdicts, reports, evidence pack | Local Cursor Canvas binaries |
| `tools/` env checks, etc. | Assumptions that require changing system package versions |

Searchable title: **Questioning Fundamental Principles of Organic Chemistry** (**Professor Zhongheng Yu**). Once obtained, public definitions and numerical tables can be checked against this repo’s criteria and results; **no** electronic copy is provided here.

---

## Key entry points

| Document | Description |
|----------|-------------|
| [docs/propositions.en.md](docs/propositions.en.md) | Authoritative proposition specs ([中文](docs/propositions.md)) |
| [docs/quality_gates.en.md](docs/quality_gates.en.md) | Hard evidence gates G1–G5 ([中文](docs/quality_gates.md)) |
| [**docs/FROZEN_VERDICT_AUTHORITY.en.md**](docs/FROZEN_VERDICT_AUTHORITY.en.md) | **Frozen verdict authority (zero ambiguity)** ([中文](docs/FROZEN_VERDICT_AUTHORITY.md)) |
| [**docs/expert_quick_review_guide.md**](docs/expert_quick_review_guide.md) | **Quick review guide for quantum-chemistry experts (EN/中文)** |
| [deliverables/final/README.en.md](deliverables/final/README.en.md) | Freeze pack index ([中文](deliverables/final/README.md)) |
| [deliverables/unit/INDEX.en.md](deliverables/unit/INDEX.en.md) | Unit verdicts P1–P9 ([中文](deliverables/unit/INDEX.md)) |
| [docs/research_plan.md](docs/research_plan.md) | Implementation plan (中文) |
| [docs/directory_structure.md](docs/directory_structure.md) | Directory & deposit conventions (中文) |
| [WSL2/README.md](WSL2/README.md) | **WSL2 compute plane** (Ubuntu + GPU archive) |
| [docs/github_upload_plan.md](docs/github_upload_plan.md) | GitHub release plan (中文) |
| [docs/RELEASE_NOTES.md](docs/RELEASE_NOTES.md) | **v1.0.0 freeze Release notes** |
| [implementation/](implementation/) | Plan tables, env/tools, AI4S mapping |
| [NOTICE](NOTICE) | Copyright & independent-implementation notice |
| [README.zh-CN.md](README.zh-CN.md) | Chinese homepage |

---

## Directory layout (GitHub)

Authoritative note: [docs/directory_structure.md](docs/directory_structure.md)

```
AI4S-OrgChem/
├── README.md · README.zh-CN.md · NOTICE · LICENSE · CITATION.cff · .gitignore
├── docs/                 # specs, quality gates, plans, structure
├── implementation/       # pre-registration plans & tracking
├── src/                  # homemade code: common / localization / p1…p9
├── results/Pn/           # tables (tracked) · invalid_* (audit) · raw (not tracked)
├── deliverables/
│   ├── unit/Pn/          # VERDICT · report · evidence · canvas_link
│   ├── final/            # five-piece freeze pack
│   └── papers/           # chemistry / AI4S method skeletons
├── WSL2/                 # ★ compute plane: Ubuntu-24.04 + GPU archive
├── tools/                # env checks · wsl_snapshot.sh
├── source/README.md      # pointer only; PDF not in repo
└── data/README.md        # pointer only; extract not in repo
```

---

## Reproduction (environment)

**Primary compute is on WSL2** (Ubuntu 24.04, PySCF 2.14.0, GPU visible). The Windows tree and WSL mount are the same directory: `D:\AI4S-OrgChem` ≡ `/mnt/d/AI4S-OrgChem`. See [`WSL2/README.md`](WSL2/README.md).

```bash
# Inside WSL2
cd /mnt/d/AI4S-OrgChem
PYTHONPATH=. python3 -m src.p6_ese.run_v4_objections --basis '6-31g*'
```

Full entry commands: [FULL_REPORT.md §7](deliverables/final/FULL_REPORT.md).  
If `quality_gate.passed != true`, scripts exit non-zero and do not suggest a verdict.

---

## Citation and license

- **Book author (subject of verification)**: **Professor Zhongheng Yu**
- **Copyright (c) 2026 Xiao Chen** · e-mail: [chenxiao0101@gmail.com](mailto:chenxiao0101@gmail.com)
- Cite this replication arm: [`CITATION.cff`](CITATION.cff)
- Code license: [`LICENSE`](LICENSE) (MIT)
- Copyright in the monograph remains with **Professor Zhongheng Yu** and the publisher; this repository does not redistribute the book materials.
