# Full Report: Independent Replication Study

**English** | [中文](FULL_REPORT.md)

**Subject**: Nine core propositions (P1–P9) in **Professor Zhongheng Yu**’s **Questioning Fundamental Principles of Organic Chemistry** that challenge traditional organic structural theory  
**Book author (subject of verification)**: **Professor Zhongheng Yu**  
**Execution**: AI Agent + WSL PySCF 2.14.0; pre-registered criteria + hard evidence gates  
**Freeze date**: 2026-08-25  
**Verdicts**: [`VERDICT_TABLE.en.md`](VERDICT_TABLE.en.md) · **Summary**: [`EXECUTIVE_SUMMARY.en.md`](EXECUTIVE_SUMMARY.en.md) · **Checklist**: [`checklist.en.md`](checklist.en.md)

---

## 1. Genre and positioning

This is an **AI-for-Science (AI4S) Agent** project using state-of-the-art AI **and** an independent **replication / verification study**—not advocacy, not book promotion, and not “LLM verbal agreement or disagreement.”

**AI4S Agent** means a frontier LLM agent, under project Rules and hard evidence gates, decomposes propositions, writes code independently, orchestrates PySCF, enforces QC, and writes `VERDICT` files. Scientific authority comes from the gate-authorized evidence chain, not prose style.

It is one arm of a **three-LLM-Agent independent verification** plan: all verdicts and implementations in this repo were produced independently by **this arm’s AI4S Agent**; the other two arms should redo the work in separate repos with different agents, then cross-compare. Numbers and conclusions here **bind this arm only**.

Three deliberate design choices:

| Design | Purpose |
|--------|---------|
| **Criteria first** | Variables, observables, thresholds, and Agree/Disagree conditions written in `implementation/Pn_*` plans **before** computing; no post-hoc loosening |
| **Independent implementation from published definitions only** | No contact with the author’s programs; tests **definition-level** replicability, not re-running existing code |
| **Hard evidence gates** | G1–G5 (topology/geometry/convergence/scale/paths): failure → **Pending** only; no verdicts on unconverged or wrong molecules |

### 1.1 Verdict semantics (critical)

Agree / Disagree are **replicability** judgments:

- **Agree** = under Yu’s own definitions and reference-state protocols, an independent implementation reproduces signs and magnitudes and withstands challenges on reference/method choice.
- **Disagree** = independent results conflict in sign, magnitude, or key conclusion, or the claim is not reproducible under its own definitions.

Agree does **not** mean the decomposition is uniquely correct physically, nor that textbook theory is overturned. That distinction runs through the whole report.

## 2. Methods

### 2.1 Computational settings

| Item | Setting |
|------|---------|
| Engine | PySCF 2.14.0 in **WSL2** (Ubuntu 24.04); GPU (RTX 4060 Laptop) visible; workdir `/mnt/d/AI4S-OrgChem` ≡ Windows project root. Fingerprints in `WSL2/inventory/`. No software/config changes |
| Main level | B3LYP/6-31G* (P2 LFMO-lite: RHF/STO-3G, matching Yu’s early π–σ decomposition level) |
| Cross-checks | RHF, MP2//RHF, ZPE (P1); dual basis (P8); 2011-lite matrix zeroing (P5/P6/P9) |
| Units | Internal Hartree; report kcal/mol; stabilization negative, destabilization positive |
| Optimization | BFGS + analytic gradients; localized SCF: DIIS → Newton fallback |

### 2.2 Homemade localization family

All implemented from published PDF definitions in `src/localization/`:

| Module | Published definition | Used in |
|--------|----------------------|---------|
| `gl_2007.py` | Ch6 2007: π-AO ID + inter-fragment π–π Fock/overlap block deletion; GE-m merges one allowed conjugated double-bond pair | P1, P5, P6, P9 |
| `hetero_gl.py` | Heteroatom (O/N) π-fragment localization | P8 |
| `plg.py` | Partial localization: exclude specified ring–ring π Fock/S/K couplings | P7 |
| `lfmo_ao_proxy.py` | **AO proxy** of Ch4–5 Fig. 5-15 four-state definitions (not a Kost LFMO reproduction) | P2 |

`zero_exchange` implements “2011-lite”: additionally zero inter-fragment exchange blocks for sign tests. Strict 2011/2014 two-electron exchange deletion is **not** implemented (shared residual).

### 2.3 Quality gates G1–G5

| Gate | Standard | On failure |
|------|----------|------------|
| G1 Topology | Bonding, coordination, ring size asserted; empty defect list | Abort; void old results |
| G2 Geometry | No atom clash (dmin ≥ 0.85 Å); key bonds in chemical window | Same |
| G3 Convergence | Key points converged **or** independent proxy (e.g. bidirectional PES hysteresis ≤1 kcal/mol) | No verdicts from unconverged curves |
| G4 Energy scale | ΔE in pre-declared chemical window | Reject this round |
| G5 Clean paths | Formal tables in `results/Pn/tables/`; voided data only in `invalid_*/` | Migrate; VERDICT must not cite invalid |

Any failure → remain **Pending** (**Pending ≠ verification failure**).

## 3. Results overview

Measured values, thresholds, and pass/fail: [`VERDICT_TABLE.en.md`](VERDICT_TABLE.en.md). Here only cross-proposition structure.

### 3.1 Successful evidence clusters

**Cluster A — Parameter-free energy benchmarks (P6, P8)**  
Strongest cluster. P6 hits two **opposite-sign** benchmarks with one localization + virtual-reference pipeline: benzene ESE = −35.44 (stabilizing); CBD ΔEA = +53.98 (destabilizing). P8 extends to heterocycles: furan-like ΔEA>0 vs benzene ΔEA<0. Persuasion comes from “one protocol, multiple sign directions, no tunable parameters.”

**Cluster B — Geometric attribution (P4, P7)**  
Energy components and localization answer what drives geometry. P4: toward equal bonds ΔEN = −96, ΔEe = +88 kcal/mol (nuclear repulsion favors D₆ₕ; electronic term favors alternation). P7: after cutting π coupling, central BLA 0.207 → 0.020 Å—contrary to Mills–Nixon angle-strain dominance.

**Cluster C — Conformational counterexamples (P3, P2)**  
P3: crowded large twist can be most stable (θ_min≈44.9°; free opt ≈34.5°). P2: on NBA, destabilizing terms fall with twist (dEV/dθ<0, dEnσσ/dθ<0)—they not only destabilize but **drive** distortion.

**Cluster D — Trend convergence (P5, P9)**  
P5: ΔEAm 6/6 positive with single-bond lengthening; P9: \|VDE/π\| 4n/4n+2 gap converges 2.836 → 0.424 with N. Pattern evidence, not single-point luck.

### 3.2 The non-reproduced claim: P1

Yu’s entry argument: classical CE ≈ −3.9 vs two 1-butenes (stabilizing); vs trans-2-butene CE ≈ +1.9 (destabilizing)—so “conjugation-energy sign depends on reference choice.”

This project (B3LYP/6-31G*, kcal/mol):

| Quantity | This work | Yu |
|----------|-----------|-----|
| CE₁ (vs 1-butene) | **+9.331** | −3.9 |
| CE₂ (vs trans-2-butene) | **+1.861** | +1.9 |

CE₂ is close; CE₁ is **same sign (positive)**—no flip. Four levels agree → not a method accident.

Structure: flip requires `2|ΔH(trans-2)| < |ΔH(butadiene)| < 2|ΔH(1-butene)|`. Electronic + ZPE energies here keep \|ΔH(butadiene)\| below both doubled hydrogenation heats, so CE₁ and CE₂ share sign. Literature **experimental** heats can reproduce the narrative arithmetic (see `unit/P1/report.md` §5), but project boundaries bind formal verdicts to independent QC → **Disagree**.

L3: homemade 2007-GL ΔE = +4.06 kcal/mol (destabilizing) but not ≈+1.4 (Yu Ch12 notes +1.4 from 2014 with exchange deletion) → magnitude Disagree as well.

**P1 vs P2**: P1 rejects thermochemical **entry narrative** replicability, not localization-level destabilization. P2 treats P1 as a **disclosed gap**, neither using it to overturn localization evidence nor rewriting P1 from that evidence.

## 4. Expert objections and responses

### 4.1 “Energy decomposition is non-unique; your protocol favors Yu”

**Accept the premise; reject the inference.** Reference-state dependence is a known weakness. Our handling:

1. All protocols (GL / GE-m / ESE / VR / PLG / LFMO-lite) are explicit in scripts and JSON;
2. Each proposition carries level/basis sensitivity; signs must not flip under sensitivity to count as Agree;
3. Verdict semantics are limited to “reproducible under that protocol,” not “protocol is true.”

Changing protocol may change numbers; overturning sign-level conclusions here requires an equally coherent, parameter-free alternative that **also** explains both opposite P6 benchmarks—a testable challenge, not rhetoric.

### 4.2 “Eight Agrees—is the AI pleasing the book?”

Three counters:

1. **P1 is Disagree**, and it is Yu’s entry argument;
2. **P3 was rejected by our own gates**—wrong-geometry “Disagree” readings were voided, archived, and rerun with a full trail;
3. Each proposition keeps a **rejection list** (e.g. P2 θ=45° EV; P9 N=16 BLA)—not every number is evidence.

### 4.3 “Does this overturn conjugation stabilization theory?”

**No; this report does not claim that.** Traditional evidence includes spectroscopy, barriers, bond-length trends, reactivity, etc. This project only tests the part Yu challenges—mainly **thermochemical and energy-decomposition stabilization arguments**. What can be said:

- The common inference “conjugation stabilization makes coplanar most stable” fails for NBA, polyenes, and strained aromatics as treated here;
- Local/pairwise conjugation under a localized reference is independently reproducible as destabilizing;
- Benzene-level extra delocalization remains stabilizing (ESE<0)—Yu treats this separately.

### 4.4 “STO-3G / 6-31G* too small—basis artifact?”

P8 dual basis (\|ΔLDE\|=0.68 kcal); P9 2011-lite and BLA sensitivities (no sign flip); P1 four-level cross-check. P2 LFMO-lite uses STO-3G to match Yu’s early level; its main evidence (P3–P9) stands independently at B3LYP/6-31G*. Larger basis / higher correlation remain residuals; no sign of basis-dependent signs.

### 4.5 “Without the author’s code, how can you claim replication?”

Precisely because of that. Author’s code only tests “same code, same number.” Independent implementation from published definitions tests “claims hold at the definition level.” Cost: digit-by-digit match is not expected (P2 LFMO-lite explicitly not vs Kost tables; P9 VDE is 2007 ESE proxy, not EV(2011))—caveats are written in each `VERDICT.md`.

## 5. Limitations and residuals

### 5.1 Shared residuals

- **Strict 2011/2014 exchange deletion not implemented.** Several props use 2011-lite for sign tests. P1 Disagree on +1.4 magnitude is directly related.
- Higher correlation (CCSD etc.) and larger bases not systematically rolled out.
- Experimental comparisons only as P1 appendix thermochemical arithmetic.

### 5.2 Per-proposition residuals (do not change verdicts)

| Prop. | Residual |
|-------|----------|
| P1 | 2011/2014 exchange deletion |
| P2 | Eπσ vs θ resistance slope not isolated; true LFMO not implemented |
| P3 | — |
| P4 | Optional full GL→G and 2011 localization contrasts |
| P5 | Exchange deletion / formal ESE moved to P6 |
| P6 | Minor deepenings |
| P7 | More molecules; fully free PLG opt; denser Δr grid |
| P8 | Broader series (e.g. imidazole); strict 2011 exchange deletion |
| P9 | N=20–26 full B3LYP; nonplanar large rings |

## 6. Compliance boundaries

- **No** use, copy, decompile, or translate of the author’s programs; published math, orbital diagrams, energy expressions, and numerical tables only.
- **No** changes to WSL/Windows software versions or GPU/driver config; **no** installs or upgrades.
- All primary verdict energies from electronic structure; no ML potentials (MACE etc.) or symbolic regression (PySR).
- Scope locked to P1–P9; PBH-type extrapolations not promoted to separate verdicts.
- Formal verdicts only in `deliverables/unit/Pn/VERDICT.md` and this `final/` directory; chat does not count.

## 7. How to reproduce

```bash
# Environment: WSL, PySCF 2.14.0; no new packages required
cd /mnt/d/AI4S-OrgChem

PYTHONPATH=. python3 -m src.p1_butadiene.run --sensitivity --basis '6-31g*'
PYTHONPATH=. python3 -m src.p1_butadiene.run_zpe
PYTHONPATH=. python3 -m src.p1_butadiene.run_gl
PYTHONPATH=. python3 -m src.p3_nba.run_tight
PYTHONPATH=. python3 -m src.p4_benzene.run
PYTHONPATH=. python3 -m src.p5_local.run_v5c_objections
PYTHONPATH=. python3 -m src.p6_ese.run_v4_objections --basis '6-31g*'
PYTHONPATH=. python3 -m src.p7_strained.run_v1c --basis '6-31g*'
PYTHONPATH=. python3 -m src.p8_furan.run_v2 --basis '6-31g*' --basis2 '6-31g'
PYTHONPATH=. python3 -m src.p9_annulene.run_v2b
PYTHONPATH=. python3 -m src.p2_aggregate.run_lfmo_lite
PYTHONPATH=. python3 -m src.p2_aggregate.summarize
```

Scripts write `results/Pn/tables/*.json` (with `quality_gate`) and `summary*.md`. If `quality_gate.passed != true`, exit non-zero; no verdict suggestion.

## 8. Document map

| Layer | Location |
|-------|----------|
| Proposition specs & verdict rules | `docs/propositions.en.md` / `propositions.md` |
| Quality-gate specs | `docs/quality_gates.en.md` / `quality_gates.md` |
| Plans & tracking | `docs/research_plan.md`, `implementation/` |
| Book claim map | `docs/book_overview.md` |
| Unit verdicts & technical reports | `deliverables/unit/Pn/{VERDICT,report}.md` + `evidence/` |
| Interactive viz | Cursor Canvas (local; not on GitHub) |

**Copyright (c) 2026 Xiao Chen** · chenxiao0101@gmail.com
