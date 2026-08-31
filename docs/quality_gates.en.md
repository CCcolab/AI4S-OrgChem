# Evidence quality gates

**English** | [中文](quality_gates.md)

> **Purpose**: Prevent Agree/Disagree verdicts built on “unconverged + wrong molecule.”  
> **Trigger**: P3 (2026-08-23)—under-converged false PES nearly wrote Disagree; `_phenyl` ring bug meant the system was not NBA.  
> **Force**: Mandatory checklist before writing `deliverables/unit/Pn/VERDICT.md`; Cursor Agent rule `.cursor/rules/ai4s-quality-gates.mdc` (`alwaysApply: true`).

---

## 1. One sentence

> **If the molecule is wrong or the calculation is not defensibly converged, no final verdict—Pending only; reject/archive bad data.**

---

## 2. Five gates before a verdict (all required)

| # | Gate | Pass standard (examples; props may tighten) | On failure |
|---|------|-----------------------------------------------|------------|
| G1 | **Chemical graph / topology** | Bonding, coordination, ring size match target; `check_topology` (or equivalent) returns empty defects | Abort; fix geometry; mark old results void |
| G2 | **Geometry integrity** | No atom clash (suggest dmin ≥ 0.85 Å); key bonds in a reasonable window | Same |
| G3 | **Defensible convergence** | Key points converged; **or** independent proxy (e.g. bidirectional PES hysteresis ≤ 1 kcal/mol) | Raise maxiter / change protocol; **do not** verdict from unconverged curves |
| G4 | **Energy scale** | Conformational/control ΔE in chemical window (declared in plan; typical ≲ tens of kcal; 10²-scale → reject) | Reject this PES/table; check unrelaxed bonds / wrong molecule |
| G5 | **Clean evidence paths** | Formal tables in `results/Pn/tables/`; voided data only in `results/Pn/invalid_*/` | Migrate bad data; VERDICT must not cite invalid |

**Only when G1–G5 all pass** may `VERDICT.md` state **Agree** or **Disagree**.  
Any failure → stay **in progress / Pending** (Pending ≠ verification failure).

---

## 3. Negative example (P3)

| Stage | Error | Consequence |
|-------|-------|-------------|
| First relaxed run | maxiter too short; ΔE ~10² kcal | Read E_min≈75° → nearly wrong **Disagree** |
| Tightened round 1 | Ring bug → seven-membered N-heterocycle | E_min≈90°, free opt ≈ perpendicular — **wrong molecule** |
| After fix | Six-membered + topology gate + bidirectional hysteresis + ΔE~2 kcal | E_min≈44.9°, free B3LYP≈34.5° → **Agree** |

Lesson: **the flip came from correcting methodology errors, not reinterpreting the same evidence.**

---

## 4. Code and deposit conventions

### 4.1 Every geometry-related proposition package must

1. Implement topology/bonding assertions (P3 reference: `src/p3_nba/geometry.py` → `check_topology`).  
2. Call them **after build, after pre-relax, after every scan point**; failure → `SystemExit` (or hard fail).  
3. Result JSON includes `quality_gate` (or equivalent) with at least:  
   - `passed: bool`  
   - topology / dmin / θ or coordinate error / ΔE span / convergence or hysteresis  
4. When `passed == false`: `agree` (or auto-verdict flag) must be `null`; scripts must **not** print “suggest Agree/Disagree.”

### 4.2 Voided data

```
results/Pn/invalid_<reason>/   # e.g. invalid_wrong_geometry/
```

- After migration, main `tables/` / `summary*.md` must not treat that data as formal evidence.  
- `VERDICT.md` / `report.md` may **name** the void and explain why, but must not use its numbers for the final verdict.

### 4.3 Unit definition-of-done add-on

On top of `docs/directory_structure.md` DoD:

- [ ] G1–G5 checked in compute/argumentation tables  
- [ ] `quality_gate.passed == true` (or written equivalent)  
- [ ] If invalid archives exist, VERDICT states they are not cited  

---

## 5. Minimum geometry self-checks (starter list)

| Prop. | Minimum checks |
|-------|----------------|
| P1 | Atom count/formula; key C–C windows; no overlap after opt |
| P3 | Both phenyls six-membered; ipso coordination=3; N–C / C=N present; θ-defining atoms not collinear |
| P4 | Benzene six-membered; BLA deformation does not break topology |
| P5–P9 | Localization/macrocycles: ring size and valence; virtual references must not silently become another molecule |

When starting a new proposition, copy this row into that plan’s “quality gates” section.

---

## 6. Agent / human operating discipline

1. Completeness slogans (e.g. “80%”) **cannot** replace G1–G5.  
2. Chat “looks supportive/opposed” is **not** a verdict; only `VERDICT.md` counts.  
3. Topology or scale anomalies → **stop, archive, fix gates**, then rerun; never “verdict first.”  
4. Missing packages (e.g. `geometric`) → record the gap and ask the user; do not skip G3 because an optimizer “won’t install.”

---

## 7. Maintenance

| Document | Role |
|----------|------|
| **This file** | Authoritative quality-gate spec |
| `.cursor/rules/ai4s-quality-gates.mdc` | Agent hard reminder |
| `docs/propositions.md` §1 | Verdict rules cite these gates |
| `implementation/` AI4S mapping | QC step in the seven-step flow |
| Each `implementation/Pn_*` | Proposition-level thresholds & checkboxes |

Sync the table above on revision; major incidents should add a “negative example” subsection.

---

## 8. Frozen verdicts vs gates

- G1–G5 are **preconditions** for writing `VERDICT.md` (**Agree / Disagree**), not verdicts themselves.  
- The **2026-08-25 frozen tally** and lawful change procedure: [`FROZEN_VERDICT_AUTHORITY.en.md`](FROZEN_VERDICT_AUTHORITY.en.md).  
- Future G6–G11 (provenance, formal consistency) strengthen engineering; they **do not** auto-overturn frozen L1 without §5 procedure.

---

## 9. GL / PLG energy convention (formal note)

Homemade `gl_2007.py`, `plg.py`, etc. solve modified SCF with altered Fock/overlap blocks; reported energies use converged `kernel()` values.

**In this project**:

- Treat as **operational estimators** per the book’s public matrix-deletion rules;
- **Not** automatically variational constrained energies of the standard Hamiltonian unless a full Lagrangian and consistent gradients are proven separately;
- L1 **Agree** when preregistered thresholds are met; interpretive text must separate “reproducing the book’s estimator” from “unique physical causation.”

This section **does not change** frozen L1 verdicts (2026-08-25).
