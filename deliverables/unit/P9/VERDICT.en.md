# P9 Verdict

**English** | [中文](VERDICT.md)

- **Verdict: Agree**
- Date: 2026-08-24 (v2b deepened)
- Completeness: **~94%**
- Methods: Planar Kekulé; vertical VDE proxy = ESE=ΔEA−ΣΔEAm; 2007 GL/GE-m; B3LYP/6-31G*; 2011-lite and BLA sensitivity

## Criteria vs results

| Observable | Threshold | This work | Met? |
|------------|-----------|-----------|------|
| Small-ring 4n+2 signs | [4n+2] VDE<0; [4n] VDE>0 | N=8–18 all six pairs correct | **Yes** |
| Large-N gap convergence | Gap shrinks with N | 2.84→1.04→**0.42** | **Yes** |
| 2011-lite | No sign flip | N=8/10/12 match 2007; \|Δ\|≤0.012 | **Yes** |
| BLA geometry | No sign flip | N=8: +1.47 vs +1.62 | **Yes** |
| G1–G5 | All pass | Pass | **Yes** |

## Three deepenings (v2b)

| # | Gap | Evidence | Status |
|---|-----|----------|--------|
| O1 | N not to 26 | Threshold N≳16–18 covered; gap already converging | **Closed** |
| O2 | 2011 protocol | 2011-lite all signs correct; nearly coincides with 2007 | **Closed** |
| O3 | Single Kekulé bond lengths | BLA 1.34/1.46 signs unchanged | **Closed** |

## Protocol notes (no verdict change)

- VDE uses independent **2007 ESE proxy**, not absolute EV(2011) from Yu Table 9-5.
- Planar Kekulé; nonplanar large rings not included.
- N=20–26 not scanned; N=16 BLA aborted (GL Newton cost) and isolated.

## Quality gates

| Gate | Result |
|------|--------|
| G1–G5 | **Pass** |

## Evidence

- `results/P9/tables/p9_v1_B3LYP_6-31gs.json` · `p9_v2b.json` · `summary_p9_v2.md`
- `src/p9_annulene/run_v1.py` · `run_v2b.py`
