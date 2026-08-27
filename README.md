# AI4S-OrgChem

**English** | [中文](README.zh-CN.md)

**This is an AI-for-Science (AI4S) Agent project** using state-of-the-art LLM agent technology: an agent end-to-end decomposes propositions, independently writes code, orchestrates quantum-chemical calculations, enforces hard evidence gates, and writes formal verdicts—rather than using a chat model to “comment” on a book.

**Scientific task**: Independent third-party replication of nine core propositions (P1–P9) from **Professor Zhongheng Yu**’s *Questioning Fundamental Principles of Organic Chemistry* that challenge traditional organic structural theory. Each proposition receives an **Agree / Disagree** verdict. The primary computation engine is **PySCF** in a fixed existing environment.

> **This repository does not redistribute the book PDF or full-text extract.** Please obtain the monograph yourself through lawful channels. Local workspace folders `source/` and `data/` are **never** uploaded to GitHub.

---

## What this is: an AI4S Agent, not “AI-written reports”

| Dimension | This project |
|-----------|--------------|
| **Paradigm** | **AI for Science (AI4S)**: frontier LLM agents drive an auditable research workflow |
| **Execution** | **Agent loop**: read specs → write/edit code → run PySCF → pass G1–G5 gates → write `VERDICT`; on failure, isolate voided data and rerun |
| **Stack** | Latest-generation LLM agent (this arm: **Cursor · Grok**) + electronic structure (PySCF) + project Rules/gates (not one-shot prompts) |
| **Outputs** | Reproducible scripts, result tables with `quality_gate`, unit verdicts, five-piece freeze pack—all re-runnable by others |
| **Deliberately not done** | No ML potentials / symbolic regression as primary evidence; no chat conclusions as formal verdicts; no copying of the author’s programs |

Authority of scientific claims comes from a **gate-authorized evidence chain**. The AI4S Agent’s value is running that chain end-to-end, leaving a trail, and making it reviewable.

---

## Background

Conjugation stabilization and aromaticity energy criteria are foundational statements in organic structural theory. **Professor Zhongheng Yu**’s monograph systematically challenges several of these “fundamentals”—e.g. local/pairwise conjugation as destabilizing under a localized reference, geometric distortion attributable to π rather than strain, large annulenes trending polyene-like, and furan-like systems not being scored with benzene-like ESE. These claims rest heavily on the author’s localization-based energy decompositions, and **had lacked independent third-party tests implemented from published definitions**.

This project turns those challenges into nine AI4S-Agent-executable propositions (P1–P9): criteria and thresholds are locked before computing (pre-registration); PySCF scripts are **written independently by the agent** from public math, orbital diagrams, and energy expressions in the PDF; formal verdicts are allowed only after hard evidence gates (topology / geometry / convergence / energy scale / clean paths).

**Genre**: a replication / verification study embedded in an AI4S Agent workflow—not advocacy and not promotional material for the book. **Agree** means reproducible under Yu’s own definitions; it does **not** mean that decomposition is the unique correct physical picture, nor that textbook theory is overturned.

### Three-LLM independent verification plan

The overall design is **three frontier LLM Agents**, each completing a full or equivalent verification **independently**—no shared implementation code or intermediate verdicts—then cross-comparing afterward, to reduce single-model bias and “pleasing the source” risk. That design is itself an AI4S credibility experiment.

| Arm | Environment | Status |
|-----|-------------|--------|
| **This repo** | **Cursor · Grok** AI4S Agent, independent | **Frozen** (2026-08-25): 8 Agree + 1 Disagree |
| Arm B | Another LLM Agent (separate repo / implementation) | Planned |
| Arm C | Third LLM Agent (separate repo / implementation) | Planned |

This repository is only the complete audit trail for the **Cursor Grok** arm. A cross-model consensus table will be published after all three arms are complete.

---

## Why it matters

1. **Auditable AI4S research workflow**  
   Shows how state-of-the-art AI can close a real research loop (code → compute → QC → verdict), not produce unverifiable prose. Failed gates force **Pending**—a legitimate “insufficient evidence” exit that shrinks confident wrong conclusions.

2. **Replicability before adjudicating theory**  
   First ask whether claims are third-party reproducible under published definitions; theory choice comes later. For heterodox frameworks that long relied on a single program without independent replication, that step was missing.

3. **Pre-registered criteria + hard evidence gates**  
   Criteria are locked before runs; G1–G5 failures block verdicts. P3 nearly received a wrong **Disagree** from a misbuilt molecule and under-converged PES; gates stopped it, voided data were isolated, and a rerun flipped to **Agree**—verdicts follow the evidence process, not model verbal preference.

4. **Independent implementation from public definitions only**  
   No contact with or copying of the author’s programs; testing definition-level replicability is stronger than “same code, same number.”

5. **Multi-model independent verification design**  
   Three LLM Agents in separate repos on the same proposition set provide a citable AI4S cross-check: the contest is signs and magnitudes under gates, not writing style.

6. **Restrained academic framing**  
   Reference-state dependence of energy decomposition, explicit caveats, and one formal **Disagree** (P1) are stated upfront. A replication with dissent is more credible than unanimous agreement.

---

## Results (this arm · Cursor Grok AI4S Agent)

### Overall verdict (frozen 2026-08-25)

**8 Agree + 1 Disagree.**

| Verdict | Propositions |
|---------|--------------|
| **Agree** | P2, P3, P4, P5, P6, P7, P8, P9 |
| **Disagree** | **P1** (butadiene hydrogenation-heat reference sign flip: CE₁ positive at four levels; no flip) |

One-page summary · table · full report:

- [deliverables/final/EXECUTIVE_SUMMARY.en.md](deliverables/final/EXECUTIVE_SUMMARY.en.md) ([中文](deliverables/final/EXECUTIVE_SUMMARY.md))
- [deliverables/final/VERDICT_TABLE.en.md](deliverables/final/VERDICT_TABLE.en.md) ([中文](deliverables/final/VERDICT_TABLE.md))
- [deliverables/final/FULL_REPORT.en.md](deliverables/final/FULL_REPORT.en.md) ([中文](deliverables/final/FULL_REPORT.md))
- [deliverables/final/checklist.en.md](deliverables/final/checklist.en.md) ([中文](deliverables/final/checklist.md))

### Three strongest positive lines of evidence

| # | Prop. | Key result |
|---|-------|------------|
| 1 | **P6** | Parameter-free ESE: benzene **−35.44** kcal/mol (Yu ≈ −36.3); cyclobutadiene vertical ΔEA **+53.98** (Yu ≈ +53.6)—**both opposite-sign benchmarks hit** |
| 2 | **P7** | Strained aromatic: after cutting center–periphery π coupling, central-ring BLA **+0.207 Å → +0.020 Å**, contrary to Mills–Nixon angle-strain dominance |
| 3 | **P3** | NBA relaxed PES minimum **θ ≈ 44.9°**; unconstrained optimization independently **≈34.5°**—challenges the common inference that conjugation stabilization makes coplanarity most stable |

### Methodological assets (AI4S Agent pipeline)

- **Agent end-to-end**: proposition specs → homemade localization (`src/localization/`) → batch compute → gates → verdicts, fully traceable
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

Searchable title: *Questioning Fundamental Principles of Organic Chemistry* (**Professor Zhongheng Yu**). Once obtained, public definitions and numerical tables can be checked against this repo’s criteria and results; **no** electronic copy is provided here.

---

## Key entry points

| Document | Description |
|----------|-------------|
| [docs/propositions.en.md](docs/propositions.en.md) | Authoritative proposition specs ([中文](docs/propositions.md)) |
| [docs/quality_gates.en.md](docs/quality_gates.en.md) | Hard evidence gates G1–G5 ([中文](docs/quality_gates.md)) |
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

- **Book author (subject of verification)**: **Professor Zhongheng Yu** (虞忠衡教授)
- **Copyright (c) 2026 Xiao Chen** · e-mail: [chenxiao0101@gmail.com](mailto:chenxiao0101@gmail.com)
- Cite this replication arm: [`CITATION.cff`](CITATION.cff)
- Code license: [`LICENSE`](LICENSE) (MIT)
- Copyright in the monograph remains with **Professor Zhongheng Yu** and the publisher; this repository does not redistribute the book materials.
