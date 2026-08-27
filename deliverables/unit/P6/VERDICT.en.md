# P6 Verdict

**English** | [中文](VERDICT.md)

- **Verdict: Agree**
- Date: 2026-08-24 (v4 objections closed)
- Completeness: **~98%**
- Methods: 2007 Fock+S (B3LYP primary); CBD 2D adiabatic; RHF/2011-lite sensitivity

## Criteria vs results

| Observable | Threshold | This work | Met? |
|------------|-----------|-----------|------|
| Benzene ESE (vertical B3LYP) | ≈−36.3, window 30–40 | **−35.44** (6-31G*) / **−36.12** (6-31G) | **Yes** |
| Benzene ESE (semi-adiabatic BLA) | <0, \|ESE\|≈25–50 | **−40.61** | **Yes** |
| CBD conjugation/antiaromaticity | ≈+53–55, window 45–70 | vertical **ΔEA=+57.67**; **vert@G\*=+53.98**; 2D adiabatic **+65.51** | **Yes** |
| CBD ESE at G* | ≡0 (two double bonds) | **≈0** | **Yes** |
| Sign / basis | Stable | B3LYP both bases + adiabatic same sign | **Yes** |
| Butadiene additivity | ESE≈0 | **≈0** | **Yes** |

## Three objections (v4 closed)

| # | Objection | Evidence | Status |
|---|-----------|----------|--------|
| O1 | Only 1D semi-adiabatic; incomplete | CBD **2D** (r_d×r_s): G\*=(1.34,1.56), GL\*=(1.34,1.44), **not boundary-hit**; ΔEA=+65.51 (in window) | **Closed** |
| O2 | RHF-2007 blow-up / 2011 insufficient | B3LYP 07/11 ESE≈−35.44/−35.48; RHF-2011 benzene ESE=−45.46 (sign OK); RHF-2007 **pathological void** | **Closed** |
| O3 | “CBD should be ESE=+53” | Two double bonds ⇒ G≡GE ⇒ **ESE≡0**; Yu +53–55 **≡ ΔEA**; vert@G\* **+53.98** (vs Yu 53.6, error <0.4) | **Closed** |

## Quality gates

| Gate | Result |
|------|--------|
| G1–G5 | **Pass** |

## Residuals (no verdict change)

- Full Cartesian adiabatic opts; strict 2011/2014 AO exchange integrals (not K-block approx).

## Evidence

- `results/P6/tables/p6_v4_objections_B3LYP_6-31gs.json` · `summary_p6_v4.md`
- `src/p6_ese/run_v4_objections.py` · `src/localization/gl_2007.py`
