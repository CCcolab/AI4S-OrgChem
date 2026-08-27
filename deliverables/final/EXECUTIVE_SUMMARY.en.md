# Executive Summary: AI4S Agent Independent Replication of Yu’s Destabilizing-Delocalization Claims

**English** | [中文](EXECUTIVE_SUMMARY.md)

**Date**: 2026-08-25 · **Verdicts**: [`VERDICT_TABLE.en.md`](VERDICT_TABLE.en.md) · **Methods**: [`FULL_REPORT.en.md`](FULL_REPORT.en.md)  
**Project type**: An **AI-for-Science (AI4S) Agent** project using state-of-the-art AI—an agent executes verification end-to-end under pre-registered criteria and hard evidence gates, not verbal commentary.  
**Verification arm**: The **Cursor · Grok** arm of a planned three-LLM-Agent independent verification (this summary covers **this arm only**)

---

## One-sentence conclusion

Under Yu’s own published definitions and reference-state protocols, **eight of nine core propositions are independently third-party reproducible (Agree); one is not (Disagree)**. This neither overturns traditional conjugation/aromaticity theory nor establishes Yu’s energy decomposition as the unique correct physical picture.

## Scoreboard

| Prop. | Claim (short) | Verdict | Completeness |
|-------|---------------|---------|--------------|
| **P1** | Butadiene conjugation energy flips sign when the reference molecule changes | **Disagree** | ~98% |
| **P2** | Conjugation destabilizes and drives distortion (umbrella claim) | **Agree** | ~94% |
| **P3** | Crowded large-twist NBA-type conformations can be most stable | **Agree** | ~95% |
| **P4** | Benzene D₆ₕ is driven by nuclear-repulsion minimization | **Agree** | ~90% |
| **P5** | Local conjugation energy ΔEAm is always destabilizing and lengthens the single bond | **Agree** | ~96% |
| **P6** | Parameter-free ESE: benzene ≈ −36; cyclobutadiene ≈ +54 | **Agree** | ~98% |
| **P7** | Strained-aromatic BLA arises from π delocalization, not angle strain | **Agree** | ~90% |
| **P8** | Furan-like systems should use LDE, not benzene-like ESE for aromaticity | **Agree** | ~96% |
| **P9** | Large [N]annulenes trend polyene-like | **Agree** | ~94% |

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
