# Freeze checklist

**English** | [中文](checklist.md)

> Audit view. Item-by-item verifiable; no evaluative prose.  
> Freeze date: 2026-08-25 · Specs: [`docs/quality_gates.en.md`](../../docs/quality_gates.en.md), [`docs/directory_structure.md`](../../docs/directory_structure.md)

---

## 1. Verdict completeness

| # | Check | Result |
|---|-------|--------|
| 1.1 | P1–P9 each have `deliverables/unit/Pn/VERDICT.md` | ✅ 9/9 |
| 1.2 | Each unit has **L0 public status** + **L1 audit snapshot**; no Pending left | ✅ see [`FROZEN_VERDICT_AUTHORITY.en.md`](../../docs/FROZEN_VERDICT_AUTHORITY.en.md) §2 |
| 1.3 | No early final verdicts with insufficient evidence | ✅ |
| 1.4 | Formal conclusions only in `VERDICT.md` / `final/` | ✅ |

## 2. Quality gates G1–G5 (per proposition)

| Prop. | G1 Topology | G2 Geometry | G3 Convergence | G4 Scale | G5 Paths | `quality_gate.passed` |
|-------|-------------|----------------|----------------|----------|----------|------------------------|
| P1 | ✅ | ✅ | ✅ | ✅ | ✅ | true (incl. backfill) |
| P2 | ✅ | ✅ | ✅ | ✅ | ✅ | true |
| P3 | ✅ | ✅ | ✅ (hysteresis 0.28 ≤ 1) | ✅ | ✅ | true |
| P4–P9 | ✅ | ✅ | ✅ | ✅ | ✅ | true |

**Conclusion**: No verdict rests on unconverged calculations or wrong molecules.

## 3. Voided / rejected data declaration

| Prop. | Item | Disposition | Cited in verdict? |
|-------|------|-------------|-------------------|
| P3 | All old scans/free opts from ring-building bug (once read E_min≈75°/90°) | Archived `results/P3/invalid_wrong_geometry/` | **No** |
| P2 | θ=45° EV=+55.7 kcal (π/σ AO assignment unstable) | Logged but excluded from EV window (0–30°) | **No** |
| P2 | Early constrained-SCF routes (FUD/DSI nonconvergent at θ>0) | Never written to formal tables | **No** |
| P9 | N=16 BLA variant (GL Newton stall abort) | Excluded from formal tables | **No** |
| P1 | Fock-only GL diagnostic (−71.8 kcal, pathological control) | Diagnostic annotation only | **No** |

## 4. Unit delivery completeness

| Prop. | `VERDICT.md` | `report.md` | `evidence/` | `canvas_link.md` |
|-------|--------------|-------------|-------------|-------------------|
| P1–P9 | ✅ | ✅ | ✅ | ✅ (no local absolute paths) |
| P3 note | — | — | ✅ **only** `pes_tight_*`; voided geometry removed | ✅ |

## 5. Final delivery completeness

| # | File | Status |
|---|------|--------|
| 5.1 | `EXECUTIVE_SUMMARY.md` (+ `.en.md`) | ✅ |
| 5.2 | `VERDICT_TABLE.md` (+ `.en.md`) | ✅ |
| 5.3 | `FULL_REPORT.md` (+ `.en.md`) | ✅ |
| 5.4 | `checklist.md` (+ `.en.md`) | ✅ (this file) |
| 5.5 | `evidence_pack/` | ✅ |

## 6. Compliance checks

| # | Check | Result |
|---|-------|--------|
| 6.1 | No use/copy/decompile/translate of author’s programs | ✅ Homemade from published definitions |
| 6.2 | No WSL/Windows package or GPU config changes | ✅ No install/upgrade |
| 6.3 | Primary verdict energies all from electronic structure | ✅ PySCF only |
| 6.4 | No MACE/PySR as verdict basis | ✅ Not used |
| 6.5 | Scope locked to P1–P9 | ✅ PBH etc. remain extrapolations |
| 6.6 | Caveats (“not Kost LFMO”, “2007 ESE proxy”) written | ✅ See P2/P9 VERDICT & report |

## 7. Document sync

| # | Document | Status |
|---|----------|--------|
| 7.1–7.4 | `implementation/` · `docs/propositions.md` · `research_plan.md` · `unit/INDEX.md` | ✅ Aligned with VERDICTs |
| 7.5 | Cursor Canvas (P1–P9 + boards) | ✅ Local (11) |
| 7.6 | `docs/formalism.md` | ⬜ Not built (non-blocking) |

## 8. Open items (explicit; do not change closed verdicts)

| # | Item | Impact |
|---|------|--------|
| 8.1 | Strict 2011/2014 exchange deletion not implemented | Related to P1 magnitude Disagree vs +1.4; others used 2011-lite for signs |
| 8.2 | P2 Eπσ vs θ slope not isolated | Only planar Eπσ(0)=0.00; EV/Enσσ channels unaffected |
| 8.3 | P9 no N=20–26; planar Kekulé only | Threshold N≳16–18 already covered |
| 8.4 | P7 single molecule C₁₂H₆ | Conclusion limited to that molecule |
| 8.5 | CCSD etc. / larger bases not systematic | Sensitivities show no sign dependence |
