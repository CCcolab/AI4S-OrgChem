# P3 Verdict

**English** | [中文](VERDICT.md)

> **Public scientific status (L0)**: `SUPPORTED_WITHIN_SCOPE` · NBA tested path θ_min≈44.9°.  
> **Pre-registered threshold audit (L1, 2026-08-25 snapshot)**: **Agree** (historical audit; **not** the public scientific score).

- **Public status (L0)**: `SUPPORTED_WITHIN_SCOPE`
- **Threshold audit (L1)**: **Agree**
- Date: 2026-08-23
- Status: Tight rescan quality gates PASSED; independent free B3LYP/6-31G* optimum angle

## Summary

Under the agreed protocol (RHF/3-21G bidirectional sequential θ-constrained relax + B3LYP/6-31G* single points):

1. **E(θ) minimum near 44.9°** (inside 30°–60° window)  
2. **Vs near-planar**: at that point EN↑, Ee↓ (qualitative signs match Yu)  
3. **Near-planar is not the global minimum** (ΔE span only ~1.9 kcal/mol—chemically reasonable)  
4. **Unconstrained B3LYP/6-31G* free optimum**: signed θ=145.5° → fold from planar **~34.5°**, near crystallographic 36°–55°  

Relative to the traditional foil “sterics raise large-angle forms,” this tier supports Yu’s falsifiable core: crowded large-twist conformations can be most stable.

## Key numbers

| Item | Result |
|------|--------|
| Scan E_min | **44.9°** (ΔE = −0.96 kcal/mol vs θ≈0°) |
| Free B3LYP optimum (fold) | **~34.5°** |
| Book θ_exp | 36°–55° |
| ΔE span (0–90°) | 1.94 kcal/mol |
| Bidirectional hysteresis max | 0.28 kcal/mol |
| max \|θ−θ_target\| | 0.19° |
| Topology (six-membered) | Intact at all points |

Data: `results/P3/tables/pes_tight_B3LYP_6-31gs_on_RHF_3-21g_both.json`  
Free opt: `results/P3/free_b3lyp/seed_prerelaxed.xyz`

## Sensitivities / caveats

- **Ring-building was once badly wrong** (`_phenyl` → star carbon cluster → seven-membered N-heterocycle after relax). Voided data in `results/P3/invalid_wrong_geometry/` — **not used**. Current protocol has hard `check_topology`.  
- Ee/EN: signs EN↑/Ee↓ hold, but \|ΔEe\|, \|ΔEN\| still ~10³ kcal while total ΔE ~1 kcal — RHF geom + DFT SP bond-length micro-difference cancellation artifact. **Verdict rests mainly on E(θ) location and free optimum angle**, not absolute decomposition magnitudes.  
- Two angles (θ=0, 30) missed gradient thresholds at maxiter=40; bidirectional hysteresis ≤0.28 kcal/mol substitutes as G3 proxy.  

## Criteria checklist

| Criterion | This round |
|-----------|------------|
| E_min ∈ ~30°–60° | ✓ 44.9° |
| That region: EN↑, Ee↓ vs near-planar | ✓ qualitative |
| Near-planar not global min | ✓ |

## Evidence

- `results/P3/` · `deliverables/unit/P3/report.md` · Canvas `p3-nba-twist`
