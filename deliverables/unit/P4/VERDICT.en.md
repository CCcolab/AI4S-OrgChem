# P4 Verdict

**English** | [中文](VERDICT.md)

> **Frozen verdict (L1 authority)**: **Agree** · freeze **2026-08-25**.  
> Path-description scope only; causal “driving” language is **not** part of L1. See [`docs/FROZEN_VERDICT_AUTHORITY.en.md`](../../../docs/FROZEN_VERDICT_AUTHORITY.en.md).

- **Verdict: Agree**
- Date: 2026-08-23
- Completeness: **~90%** (optional full GL→G / 2011 localization contrasts not done; non-blocking)
- Methods: B3LYP/6-31G*; D₆ₕ SciPy-BFGS then scan δ=r_a−r_b (r_a=r0+δ/2, r_b=r0−δ/2); Ee = E − EN

## Criteria vs results

| Observable | Agree threshold | This work | Met? |
|------------|-----------------|-----------|------|
| E_tot minimum | δ≈0 (equal bonds) | **δ = 0.000 Å** | Yes |
| EN minimum | δ≈0 | **δ = 0.000 Å** | Yes |
| δ_max→0: ΔEN | < 0 (nuclear repulsion favorable) | **−96.0 kcal/mol** | Yes |
| δ_max→0: ΔEe | > 0 (electronic unfavorable) | **+88.3 kcal/mol** | Yes |
| Ee preference | Favors alternation (min not at δ=0) | Ee min at **δ=0.120** | Yes (qualitative) |

Conjunction → **Agree**.

## Main numbers (relative to δ=0)

| δ / Å | ΔE | ΔEN | ΔEe |
|-------|-----|------|------|
| 0.00 | 0 | 0 | 0 |
| 0.06 | +1.91 | +23.97 | −22.06 |
| 0.12 | +7.75 | +96.03 | −88.28 |

- r0 (D₆ₕ) = **1.3969 Å**; `opt_success=true`, BLA≈0  
- ΔE span 7.75 kcal (<40 → G4 pass)

## Key points

1. Under published dra=−drb scan, **equal-bond geometry minimizes both total energy and EN** (path description on frozen BLA coordinates).  
2. From strong alternation toward equal bonds: ΔEN < 0, ΔEe > 0 (signs pass threshold); Ee = E − EN is an identity decomposition—**do not** infer causal “nuclear repulsion drives D₆ₕ” from this alone (causal claim not tested).  
3. This tier does **not** reproduce Yu’s full GL→G localization tables or 2011 program decimals.  
4. Large absolute \|ΔEN\|, \|ΔEe\|: verdict trusts **signs and extremum locations** (same discipline as P3 on decomposition artifacts).

## Quality gates

| Gate | Result |
|------|--------|
| G1–G5 | **Pass** → `results/P4/tables/bla_scan_B3LYP_6-31gs.json` |

## Evidence

- `results/P4/tables/summary_bla.md` · `bla_scan_B3LYP_6-31gs.json`
- `src/p4_benzene/` · `deliverables/unit/P4/report.md`
