# P7 Verdict

**English** | [中文](VERDICT.md)

- **Verdict: Agree**
- Date: 2026-08-24 (v1c)
- Completeness: **~90%**
- Methods: B3LYP/6-31G*; benzotricyclobutadiene C₁₂H₆; G = D₃ₕ four-parameter opt; PLG = after excluding central-ring(A)–peripheral C=C(B) π Fock/S/K coupling, fix mean bond length and scan Δr

## Criteria vs results

| Observable | Threshold | This work | Met? |
|------------|-----------|-----------|------|
| Δr(G)=r_endo−r_exo | Large positive (Yu 0.177) | **+0.2065 Å** | **Yes** |
| Δr(PLG) | Near 0 (Yu −0.002) | **scan min +0.020 Å** | **Yes** |
| Significant \|Δr\| drop | drop ≫ 0.05 Å | **drop = +0.1865 Å** | **Yes** |
| PLG prefers equalization | E(PLG@Δr≈0) < E(PLG@Δr(G)) | **Yes** (ΔE ≈ −11.2 kcal/mol) | **Yes** |

## Physical reading

After cutting center–periphery π coupling, central-ring BLA falls from **0.21 Å to ~0.02 Å**, contrary to angle-strain (Mills–Nixon) dominance; supports Yu: BLA mainly from π delocalization.

## Quality gates

| Gate | Result |
|------|--------|
| G1 Topology | **Pass** (C₁₂H₆; three 4-rings + central 6-ring) |
| G2–G5 | **Pass** (G and all PLG scan points SCF-converged; PLG DIIS→Newton) |

## Residuals (no verdict change)

- More molecules (benzocyclobutadiene); fully free endo/exo PLG opts; denser Δr grid (could push PLG* from +0.02 toward Yu −0.002).

## Evidence

- `results/P7/tables/p7_v1c_B3LYP_6-31gs.json` · `summary_p7_v1c.md`
- `src/p7_strained/run_v1c.py` · `src/localization/plg.py`
