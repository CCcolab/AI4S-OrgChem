# P1 Verdict

**English** | [中文](VERDICT.md)

- **Verdict: Disagree**
- Date: 2026-08-23 (maintained after ZPE + L3 homemade GL)
- Status: L1/L2 main criteria + RHF/MP2/ZPE sensitivity closed; `quality_gate` logged; L3 homemade 2007-GL measured with **sub-conclusion** on “+1.4” (does not rewrite L1/L2). Completeness **~98%** (2011/2014 exchange deletion not implemented → P5/P6 residual).
- Methods: RHF / B3LYP / MP2 (MP2 SP on RHF geometries), 6-31G*; SciPy-BFGS + PySCF analytic gradients (no geometric/berny; no env changes)
- Definition (same sign as textbook hydrogenation differences; stabilization negative):

```
ΔH_hyd(M) = E(n-butane) − E(M) − n_H2 · E(H2)
CE_ref    = ΔH_hyd(butadiene) − 2 · ΔH_hyd(ref)
```

## Criteria vs results (primary: B3LYP)

| Observable | Agree threshold | This work | Met? |
|------------|-----------------|-----------|------|
| CE1 (vs 1-butene) | < 0 (stabilizing) | **+9.331 kcal/mol** | No |
| CE2 (vs trans-2-butene) | > 0 (destabilizing) | **+1.861 kcal/mol** | Yes |
| CE2 magnitude | ~ +1–3 kcal/mol | +1.86 | Yes |
| **Sign flip** CE1<0 and CE2>0 | Required | **Absent** (both positive) | **No** |

## Method sensitivity (same basis; does not change verdict)

| Level | CE1 | CE2 | Sign flip |
|-------|-----|-----|-----------|
| RHF/6-31G* | +8.130 | +2.449 | No |
| B3LYP/6-31G* | +9.331 | +1.861 | No |
| MP2//RHF/6-31G* | +8.503 | +1.858 | No |
| B3LYP + ZPE (E0K) | +8.540 | +0.511 | No |

Electronic energies at three levels and after ZPE: **CE1 and CE2 both positive**; sensitivity does **not** reproduce the pre-registered flip → **Disagree**.

## Key points

1. Under pre-registered QC criteria, **cannot** Agree: no “CE sign flips when reference molecule changes.”
2. CE2 near Yu’s +1.9 at B3LYP/MP2, but CE1 far from classical −3.9 and same sign as CE2 → full claim unsupported.
3. RHF→DFT→MP2→ZPE same-sign structure is robust (not a B3LYP/ZPE accident).
4. **Aside (no verdict change)**: literature experimental heats can change the narrative arithmetically; formal verdict stays bound to independent QC.
5. **L3 sub-conclusion** (independent of L1/L2): homemade GL per Ch6 2007 (Fock+overlap π–π block deletion), butadiene ΔE=E(G)−E(GL)≈**+4.06 kcal/mol** (destabilizing), **not** within Yu’s **+1.4** (Ch12: 2014 method) → **Disagree** on “≈+1.4”; exchange deletion not done.

## Quality gates

| Gate | Result |
|------|--------|
| G1–G5 | **Pass** (see Chinese VERDICT for evidence rows) |

Audit: `results/P1/tables/p1_quality_audit.json`

## Evidence paths

- `results/P1/tables/sensitivity_methods.json`, `ce_summary.csv`
- `ce_RHF_6-31gs.json`, `ce_B3LYP_6-31gs.json`, `ce_MP2_6-31gs.json`
- ZPE: `ce_zpe_B3LYP_6-31gs.json`
- L3 GL: `gl2007_butadiene_B3LYP_6-31gs.json`
- Report: `deliverables/unit/P1/report.md`
