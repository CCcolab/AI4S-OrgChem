# Release notes — v1.0.0-freeze-2026-08-25

In-repo copy of the GitHub Release. Update the live page with:

```bash
gh release edit v1.0.0-freeze-2026-08-25 --title "v1.0.0 Freeze 2026-08-25 — 8 Agree + 1 Disagree (P1)" --notes-file docs/RELEASE_NOTES.md
```

(Omit this header block when pasting if you prefer; GitHub body below starts at the `##` heading.)

---

## AI4S-OrgChem v1.0.0 — Freeze 2026-08-25

**AI-for-Science (AI4S) Agent** arm (this repository)  
Independent third-party **PySCF** replication of nine propositions (P1–P9) from **Professor Zhongheng Yu**, **Questioning Fundamental Principles of Organic Chemistry**.

| | |
|--|--|
| **Verdict (this arm)** | **8 Agree + 1 Disagree (P1)** |
| **Engine** | PySCF 2.14.0 (WSL2), mainly B3LYP/6-31G* |
| **Genre** | Replication / verification — not advocacy; not overturning textbook theory |
| **Tag** | `v1.0.0-freeze-2026-08-25` |
| **Release** | https://github.com/CCcolab/AI4S-OrgChem/releases/tag/v1.0.0-freeze-2026-08-25 |

### One-sentence result

Under Professor Yu’s published definitions and reference-state protocols, **eight claims are independently reproducible; one (P1 thermochemical sign-flip narrative) is not.** Agree/Disagree means **replicability**, not unique physical truth.

### Start here

| Doc | Link |
|-----|------|
| Homepage (EN) | [README.md](../README.md) · [中文](../README.zh-CN.md) |
| Executive summary | [EXECUTIVE_SUMMARY.en.md](../deliverables/final/EXECUTIVE_SUMMARY.en.md) |
| Verdict table (audit) | [VERDICT_TABLE.en.md](../deliverables/final/VERDICT_TABLE.en.md) |
| Full report | [FULL_REPORT.en.md](../deliverables/final/FULL_REPORT.en.md) |
| Unit verdicts P1–P9 | [INDEX.en.md](../deliverables/unit/INDEX.en.md) |
| QC expert guide (15–30 min) | [expert_quick_review_guide.md](expert_quick_review_guide.md) |
| Proposition specs | [propositions.en.md](propositions.en.md) |
| Quality gates G1–G5 | [quality_gates.en.md](quality_gates.en.md) |

### Scoreboard

| Prop. | Verdict | Highlight |
|-------|---------|-----------|
| **P1** | **Disagree** | No CE₁/CE₂ sign flip at four levels |
| P2 | Agree | Umbrella: NBA + localization classes |
| P3 | Agree | NBA θ_min ≈ 44.9°; free opt ≈ 34.5° |
| P4 | Agree | D₆ₕ: EN favorable / Ee favors BLA |
| P5 | Agree | ΔEAm 6/6 > 0; single-bond lengthening |
| P6 | Agree | Benzene ESE −35.44; CBD ΔEA +53.98 |
| P7 | Agree | BLA +0.207 → +0.020 Å after PLG |
| P8 | Agree | Furan-like LDE vs benzene ESE sign split |
| P9 | Agree | \|VDE/π\| gap converges with N |

### What’s included

- Homemade `src/` scripts (independent of the author’s programs)
- Pre-registered criteria, hard evidence gates G1–G5, WSL2 compute notes
- Unit packs: `VERDICT.md` + `VERDICT.en.md` + `report.md` + `evidence/`
- Freeze pack in `deliverables/final/` (with English mirrors)
- Result tables under `results/Pn/tables/` (+ voided-data isolation where applicable)

### What’s not included

- Book PDF or full-text extract (`source/`, `data/` never uploaded)
- Author’s appendix / program code
- LLM Arms B & C (planned; cross-model consensus later)

### Caveats (please cite with results)

1. Replicability ≠ unique physical correctness  
2. Energy decompositions are reference-state dependent  
3. Local ΔEAm > 0 and benzene ESE < 0 are different protocols — not “all conjugation destabilizes”  
4. This repository is **one** of three planned independent LLM verification arms  

### Copyright & contact

- **Book author (subject of verification):** Professor Zhongheng Yu  
- **Software:** Copyright (c) 2026 Xiao Chen · chenxiao0101@gmail.com  
- MIT License · See `LICENSE`, `NOTICE`, `CITATION.cff`

### Reproduce (WSL2)

```bash
cd /mnt/d/AI4S-OrgChem
PYTHONPATH=. python3 -m src.p6_ese.run_v4_objections --basis '6-31g*'
```

Full entry commands: [FULL_REPORT.en.md §7](../deliverables/final/FULL_REPORT.en.md). Scripts exit non-zero if `quality_gate.passed != true`.

---

中文读者可从 [README.zh-CN.md](../README.zh-CN.md) 与结题摘要 [EXECUTIVE_SUMMARY.md](../deliverables/final/EXECUTIVE_SUMMARY.md) 开始。
