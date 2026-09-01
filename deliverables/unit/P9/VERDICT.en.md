# P9 Verdict

**English** | [中文](VERDICT.md)

> **Public scientific status (L0)**: **`PARTIAL`** · N=8–18 planar Kekulé trend; **O1 open** (N≥20, nonplanar TBD) · see [`docs/FROZEN_VERDICT_AUTHORITY.en.md`](../../../docs/FROZEN_VERDICT_AUTHORITY.en.md) §4.  
> **Pre-registered threshold audit (L1, 2026-08-25 snapshot)**: N=8–18 sub-criterion met → **Agree** (historical audit; **not** the public scientific score).

- **Public status (L0)**: `PARTIAL`
- **Threshold audit (L1)**: **Agree** (N=8–18 sub-criterion)
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
| O1 | N not to 26 / exact onset | N=8–18 gap falls; **N=20–26 not scanned**; N=16 BLA aborted; nonplanar excluded | **Open** |
| O2 | 2011 protocol | 2011-lite all signs correct; nearly coincides with 2007 | **Closed** |
| O3 | Single Kekulé bond lengths | BLA 1.34/1.46 signs unchanged | **Closed** |

## Protocol notes (scope limits; L0 remains PARTIAL)

- VDE uses independent **2007 ESE proxy**, not absolute EV(2011) from Yu Table 9-5.
- Planar Kekulé; nonplanar large rings not included.
- N=20–26 not scanned; N=16 BLA aborted (GL Newton cost) and isolated.
- **O1 open**: N=8–18 trend must **not** substitute for N≥20 exact onset or nonplanar verification.
- **These limits narrow interpretive scope; L1 records N=8–18 sub-criterion only; public L0 is PARTIAL.**

## Quality gates

| Gate | Result |
|------|--------|
| G1–G5 | **Pass** |

## Evidence

- `results/P9/tables/p9_v1_B3LYP_6-31gs.json` · `p9_v2b.json` · `summary_p9_v2.md`
- `src/p9_annulene/run_v1.py` · `run_v2b.py`
