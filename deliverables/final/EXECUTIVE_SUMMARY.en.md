# Executive Summary: AI4S Agent Independent Replication of Professor Zhongheng Yu’s Destabilizing-Delocalization Claims

**English** | [中文](EXECUTIVE_SUMMARY.md)

**Date**: 2026-08-25 · **Verdicts**: [`VERDICT_TABLE.en.md`](VERDICT_TABLE.en.md) · **Frozen authority**: [`docs/FROZEN_VERDICT_AUTHORITY.en.md`](../../docs/FROZEN_VERDICT_AUTHORITY.en.md) · **Methods**: [`FULL_REPORT.en.md`](FULL_REPORT.en.md)  
**Book author (subject of verification)**: **Professor Zhongheng Yu**  
**Project type**: **Cursor-assisted AI4S workflow** under pre-registered criteria and hard evidence gates—not chat commentary.  
**Arm role (post 2026-09-01 review)**: P1–P9 **companion independent implementation / QC case-study arm**

---

## One-sentence conclusion

Within fixed systems, state definitions, energy **estimands**, and computational paths, several anti-traditional claims receive auditable within-scope or partial support; **do not** state “8/9 book propositions verified.” **L0 five-state labels are primary** (see [`FROZEN_VERDICT_AUTHORITY.en.md`](../../docs/FROZEN_VERDICT_AUTHORITY.en.md) §2); **8 independent propositions + P2 derived** (P2 excluded from independent verification count).

## Scoreboard (L0 public · primary)

| Prop. | Claim (short) | **L0 public** | L1 audit | Completeness |
|-------|---------------|---------------|----------|--------------|
| **P1** | Composite butadiene estimands | **`PARTIAL`** | P1b Disagree | ~98% |
| **P2** | Meta aggregate (**derived**) | **`DERIVED`** | Agree | ~94% |
| **P3** | Crowded large-twist NBA minimum | `SUPPORTED_WITHIN_SCOPE` | Agree | ~95% |
| **P4** | Frozen BLA signs; causal not tested | `PARTIAL` | Agree | ~90% |
| **P5** | ΔEAm>0 in tested set | `SUPPORTED_WITHIN_SCOPE` | Agree | ~96% |
| **P6** | Parameter-free ESE dual benchmarks | `SUPPORTED_WITHIN_SCOPE` | Agree | ~98% |
| **P7** | C₁₂H₆ single-system BLA | `PARTIAL` | Agree | ~90% |
| **P8** | Furan LDE vs benzene ESE | `SUPPORTED_WITHIN_SCOPE` | Agree | ~96% |
| **P9** | N=8–18 trend; **O1 open** | **`PARTIAL`** | N=8–18 sub-criterion Agree | ~94% |

## Three strongest positive lines among the eight Agrees

1. **P6 parameter-free ESE benchmarks hit.** Benzene ESE = **−35.44** kcal/mol (Yu −36.3; experimental scale ≈ −36); cyclobutadiene vertical ΔEA = **+53.98** (Yu +53.6). Two opposite-sign benchmarks hit with the same protocol—no empirical parameters or isodesmic reaction choice. Strongest self-consistency check of the method family.

2. **P7 direct counter to Mills–Nixon.** Central-ring BLA of benzotricyclobutadiene is **+0.207 Å** in the full geometry; after cutting center–periphery π coupling per the published definition, BLA collapses to **+0.020 Å**. If angle strain drove BLA, cutting π coupling should not erase it.

3. **P3 crowded-conformation minimum is reproducible.** NBA relaxed-PES energy minimum at **θ ≈ 44.9°**; unconstrained B3LYP independently **≈34.5°**, consistent with crystallographic scale; near-planar is not the most stable. Directly challenges the textbook inference that conjugation stabilization makes coplanarity most stable.

## The one Disagree: P1

Yu’s entry narrative requires **CE₁ < 0 and CE₂ > 0** (sign flip when the reference molecule changes). At RHF, B3LYP, MP2, and B3LYP+ZPE, this project finds CE₁ always positive (+8.1 to +9.3 kcal/mol)—**no flip**. Homemade 2007-GL per Ch6 published definitions gives ΔE = +4.06 kcal/mol (same destabilizing sign) but does not fall in Yu’s +1.4 tolerance.

Hence P1 is **Disagree**. This does not overturn localization-level evidence (GL remains destabilizing), but means the **classical hydrogenation-heat entry argument that “conjugation energy can be positive or negative” does not reproduce in our independent QC replication**.

Keeping this dissent is part of credibility: a unanimous-pass pipeline cannot prove it is not endorsement.

## Three caveats that must travel with any citation

1. **Replicability, not physical uniqueness.** Agree means signs and magnitudes reproduce under Yu’s definitions—not that the decomposition is unique or optimal.
2. **Reference-state dependence of energy decomposition** is a known domain weakness; all protocols (GL / GE-m / ESE / VR / PLG) are logged with level/basis sensitivities.
3. **Yu himself splits definitions.** Local/pairwise conjugation destabilization (ΔEAm>0) and benzene extra delocalization stabilization (ESE<0) are different things. Reading this project as “all conjugation destabilizes” is a misreading.

## Methodological by-product (AI4S Agent)

Auditable pipeline: **criteria first → hard gates G1–G5 → voided-data isolation → verdict files only** (chat is not a formal verdict).

P3 is the key case: a ring-building bug (seven-membered N-heterocycle) and under-converged PES nearly produced a wrong “Disagree”; gates stopped it, bad data went to `invalid_wrong_geometry/`, and after correction the verdict became Agree. **The flip came from fixing methodology errors, not reinterpreting the same evidence.**

## Compliance

No use, copying, decompilation, or translation of the author’s programs; scripts implement published math, orbital diagrams, and energy expressions only. Primary engine: existing WSL **PySCF 2.14.0** with no environment changes. No ML potentials or symbolic regression as verdict evidence.

**Copyright (c) 2026 Xiao Chen** · chenxiao0101@gmail.com
