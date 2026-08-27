# P5 Verdict

**English** | [中文](VERDICT.md)

- **Verdict: Agree**
- Date: 2026-08-24 (v5c objections closed)
- Completeness: **~96%** (three objections closed; 2014 exchange deletion & formal ESE → **P6**)
- Methods: 2007 Fock+S; B3LYP/6-31G*

## Criteria vs results

| Observable | Threshold | This work | Met? |
|------------|-----------|-----------|------|
| ΔEAm sign | Predominantly positive | **6/6 pairs positive** | Yes |
| Butadiene Δr | >0 | **+0.018 Å** | Yes |
| Hexatriene adiabatic ΔEAm | >0 | **+0.62 / +0.97 / dense +0.98** | Yes |
| Hexatriene Δr | >0 | **dense Δr\*=+0.0045 Å** (interior) | Yes |
| Benzene local ΔEAm | >0 | **three pairs each +9.63** | Yes |

## Three objections closed (v5c)

| Objection | Evidence | Status |
|-----------|----------|--------|
| v5 negative ΔEAm | Boundary-hit geom ΔEAm=**−0.043**; interior **+0.602**; vertical **+0.446** | **Closed** (false PES; isolated) |
| Hexatriene Δr not hard-proven | Asymmetric-bridge dense scan step=0.005: r_GL\*=1.443 / r_GE\*=1.448, **Δr\*=+0.0045 Å**, ΔEAm\*=+0.982 | **Closed** |
| Benzene ΔEA<0 should Disagree | ΔEA=**−6.54**; ΣΔEAm=**+28.90**; ESE_proxy=ΔEA−Σ=**−35.44** ≈ Yu −36.3 | **Closed** (belongs to P6, not P5 failure) |

## Main numbers (kcal/mol)

| System | Quantity | Result |
|--------|----------|--------|
| Butadiene (P1) | ΔEAm / Δr | **+4.057** / **+0.018 Å** |
| Hexatriene Ci adiabatic (v5b) | ΔEAm | **+0.616** |
| Hexatriene asymmetric dense (v5c) | ΔEAm / Δr | **+0.982** / **+0.0045 Å** |
| Benzene Kekulé vertical | ΔEAm×3 / ΔEA / ESE_proxy | **+9.63×3** / **−6.54** / **−35.44** |

## Quality gates

| Gate | Result |
|------|--------|
| G1–G5 | **Pass** (v5 → `invalid_multibond_bound_hit/`) |

## Evidence

- `results/P5/tables/p5_v5c_objections_B3LYP_6-31gs.json` · `summary_p5_v5c.md`
- `src/p5_local/run_v5c_objections.py` · `deliverables/unit/P5/report.md`
