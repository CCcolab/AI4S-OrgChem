# P8 Verdict

**English** | [中文](VERDICT.md)

- **Verdict: Agree**
- Date: 2026-08-24 (v2 deepened)
- Completeness: **~96%**
- Methods: 2007 Fock+S; hetero O/N π fragments + 3×GE-m; B3LYP; vertical + semi-adiabatic + basis sensitivity

## Criteria vs results

| System | Observable | Threshold | This work | Met? |
|--------|------------|-----------|-----------|------|
| Furan (vertical) | ΔEA>0; LDE≈−39.3 | In window | ΔEA=+28.61; **LDE=−37.72** | **Yes** |
| Furan (semi-adiabatic G*) | Pattern unchanged; closer to Yu | — | **LDE=−39.04** | **Yes** |
| Pyrrole (vertical) | LDE≈−49 | — | **−48.20** | **Yes** |
| Pyrrole (semi-adiabatic) | Pattern unchanged | — | **−51.35** | **Yes** |
| Oxazole (extension) | ΔEA>0; LDE<0 (Yu≈−36.3) | — | ΔEA=+24.59; **LDE=−41.17** | **Yes** |
| Benzene (control) | ΔEA<0; ESE≈−36 | — | ΔEA=−6.54; **ESE=−35.44** | **Yes** |
| Basis sensitivity | Sign/magnitude stable | \|ΔLDE\|≪10 | 6-31G*/6-31G: **0.68** | **Yes** |
| Cross-class | Furan-like ≠ benzene-like | Sign split | **Clear split** | **Yes** |

## Three deepenings (v2 closed)

| # | Gap | Evidence | Status |
|---|-----|----------|--------|
| O1 | Only furan/pyrrole | Oxazole same pattern: ΔEA>0, LDE=−41.2 | **Closed** |
| O2 | Vertical model geometries only | After planar G* opt, vertical LDE: furan −39.04 (vs Yu −39.3) | **Closed** |
| O3 | Basis unchecked | 6-31G* vs 6-31G: \|ΔLDE\|=0.68 kcal | **Closed** |

## Physical reading

Furan-like (furan/pyrrole/oxazole) consistently satisfy **ΔEA>0 and ΔEA<ΣΔEAm**; difference LDE is negative but **should not be called ESE**; benzene **ΔEA<0** uses ESE protocol. Semi-adiabatic brings furan LDE nearly onto Yu −39.3.

## Quality gates

| Gate | Result |
|------|--------|
| G1–G5 | **Pass** |

## Residuals (no verdict change)

- Broader series (imidazole); strict 2011 exchange deletion; NICS only as side evidence.

## Evidence

- `results/P8/tables/p8_v2_B3LYP_6-31gs.json` · `summary_p8_v2.md`
- `src/p8_furan/run_v2.py` · `src/localization/hetero_gl.py`
