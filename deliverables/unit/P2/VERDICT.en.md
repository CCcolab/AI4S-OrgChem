# P2 Verdict

**English** | [中文](VERDICT.md)

> **Public scientific status (L0)**: **`DERIVED`** · meta aggregate of P1/P3–P9 + LFMO-lite; **excluded** from independent verification count (8 independent + this derived index).  
> **Pre-registered threshold audit (L1, 2026-08-25 snapshot)**: aggregate threshold met → **Agree** (historical audit; **not** the public scientific score).

- **Public status (L0)**: `DERIVED`
- **Threshold audit (L1)**: **Agree** (aggregate layer)
- Date: 2026-08-24 (v2: aggregate + LFMO-lite)
- Completeness: **~94%**
- Methods: Inherits gated P1/P3–P9 verdicts; independent AO-proxy on NBA reproduces Fig. 5-15 EV / Enσσ signs and slopes

## Protocol (read first)

This proposition is **not** “all π delocalization always destabilizes.” Yu himself separates:

- **Local/pairwise** conjugation (ΔEAm, polyene GL ΔE, furan ΔEA) → destabilizing (positive)
- **Benzene extra** delocalization (ESE) → stabilizing (negative)
- **Geometry effects** (BLA, single-bond lengthening, NBA twist) → π/electronic terms favor distortion; benzene geometry on **P4** frozen BLA path (description only—not causal verdict here)

The falsifiable core of the umbrella claim: traditional “conjugation stabilization makes coplanar most stable” fails on the molecule classes already tested.

LFMO-lite is **not** a bit-for-bit Kost localization reproduction; it is an AO-block / Fock proxy of the published state definitions (Fig. 5-15).

## Criteria vs results

| Observable | Threshold | This work | Met? |
|------------|-----------|-----------|------|
| Polyene/aromatic localization class | Local destabilization **and** drives distortion | P5: ΔEAm 6/6>0, Δr>0; P4: Ee favors BLA | **Yes** |
| NBA-type | Crowded large twist can be most stable | P3: E_min≈44.9°; free B3LYP≈34.5° | **Yes** |
| Both classes | Plan Agree criterion | Rows above | **Yes** |
| EV (π–π) | >0 and dEV/dθ<0 (0–30°) | 6.07→…→0.43 kcal; slope −0.17 | **Yes** |
| Enσσ (nonbond σ–σ) | >0 and dEnσσ/dθ<0 (0–45°) | 38.3→…→5.5 kcal; slope −0.71 | **Yes** |
| Eπσ (π–σ) | Planar ≈0 (Table 5-16) | Eπσ(0)=0.00 | **Yes** (slope not isolated) |
| Unit QC | Source G1–G5 | All inherited pass; LFMO-lite own gate pass | **Yes** |

## How unit verdicts feed P2

| ID | Unit verdict | Contribution to P2 |
|----|--------------|--------------------|
| P1 | **Disagree** | Thermochemical CE1/CE2 **no** flip; GL ΔE≈+4.06 still destabilizing. Treated as a **gap**, does not overturn localization evidence |
| P3 | Agree | Conformational counterexample |
| P4 | Agree | π(Ee) drives BLA; EN dominates equal bonds |
| P5 | Agree | Local ΔEAm destabilizing + single-bond lengthening |
| P6 | Agree | Benzene ESE≈−35.4 vs CBD ΔEA≈+54; split protocols required |
| P7 | Agree | BLA collapses after cutting center–periphery π |
| P8 | Agree | Furan-like ΔEA>0; not benzene ESE |
| P9 | Agree | Large-ring 4n/4n+2 gap converges |

## Quality gates

| Gate | Result |
|------|--------|
| G1–G5 | **Pass** (`p2_v2_lfmo_lite.json` · `p2_v2_aggregate.json`) |

θ=45° EV=55.7 excluded from EV criterion (π assignment unstable).

## Residuals (no verdict change)

- π–σ resistance slope vs θ not independently isolated with AO proxy
- P1 classical hydrogenation sign flip remains unit **Disagree**

## Evidence

- `results/P2/tables/p2_v2_lfmo_lite.json` · `p2_v2_aggregate.json`
- `src/p2_aggregate/run_lfmo_lite.py` · `src/localization/lfmo_ao_proxy.py`
- Unit `deliverables/unit/Pn/VERDICT.md`
