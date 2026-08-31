# P2 v2 — 共轭去稳定总命题（汇总 + LFMO-lite）

- gates=True lfmo_gate=True agree=True completion~94%
- two_class=True LFMO_two_channel=True
- unit verdicts **read from** `deliverables/unit/Pn/VERDICT.md` (L1)

| ID | 判定 | 去稳定 | 驱动畸变 | 要点 |
|----|------|--------|----------|------|
| P1 | 非一致 | 是 | — | CE2=+1.86；GL ΔE=+4.06 |
| P3 | 一致 | — | 是 | θ_min=44.9°；ΔE(0–90°)≈1.9 kcal |
| P4 | 一致 | 是 | 是 | E_min 与 EN_min 均在 δ=0 |
| P5 | 一致 | 是 | 是 | ΔEAm 全正；Δr>0 |
| P6 | 一致 | 是 | — | 苯 ESE=−35.44；CBD vert@G*=+53.98 |
| P7 | 一致 | — | 是 | drop Δr=+0.1865 Å |
| P8 | 一致 | 是 | — | 呋喃 ΔEA=+28.6；LDE≈−39 |
| P9 | 一致 | — | — | gap(16–18)=0.424 |

## LFMO-lite（RHF/STO-3G）

| θ | EV | Enσσ | Eπσ |
|---|----|------|-----|
| -0.0 | 6.07 | 38.31 | 0.00 |
| 10.0 | 5.60 | 34.96 | -7.52 |
| 17.0 | 4.65 | 29.43 | -21.85 |
| 30.0 | 0.43 | 16.82 | -65.72 |
| 45.0 | 55.67 | 5.48 | -188.23 |

EV 窗 0–30°；Enσσ 窗 0–45°；θ=45° 的 EV 因 π 指派失稳不纳入 EV 判据。
