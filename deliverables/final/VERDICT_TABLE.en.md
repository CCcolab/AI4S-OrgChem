# P1–P9 Verdict Table (audit view)

**English** | [中文](VERDICT_TABLE.md)

> **Pure audit** record: pre-registered criteria, measured values, threshold met?, quality gates.  
> No theory-layer commentary (see [`FULL_REPORT.en.md`](FULL_REPORT.en.md)).  
> Authoritative sources: each `deliverables/unit/Pn/VERDICT.md`. Freeze date: 2026-08-25.

## 1. Totals

| Verdict | Count | Propositions |
|---------|-------|--------------|
| **Agree** | 8 | P2, P3, P4, P5, P6, P7, P8, P9 |
| **Disagree** | 1 | P1 |
| Pending | 0 | — |

All nine unit `VERDICT.md` files present; G1–G5 passed → freeze prerequisites for final delivery met.

## 2. By proposition

### P1 — Butadiene CE reference-molecule sensitivity

| Item | Content |
|------|---------|
| Pre-registered criterion | Agree requires **CE₁<0 and CE₂>0** (CE₂ ≈ +1–3 kcal/mol) |
| Measured (B3LYP/6-31G*) | CE₁ = **+9.331**; CE₂ = **+1.861** |
| Sensitivity | RHF +8.130/+2.449; MP2//RHF +8.503/+1.858; B3LYP+ZPE +8.540/+0.511 — **no flip at four levels** |
| Sub-criterion (L3) | Homemade 2007-GL (Fock+S) ΔE = **+4.057** (destabilizing, same sign); not within Yu’s +1.4 |
| Threshold met | **No** (CE₁ sign) |
| G1–G5 | Pass |
| **Verdict** | **Disagree** (~98%) |

### P2 — Conjugation destabilizes and drives distortion (umbrella)

| Item | Content |
|------|---------|
| Pre-registered criterion | Agree requires **both** NBA-type and polyene/aromatic-localization classes to reproduce “destabilize and/or drive distortion” |
| Measured (two classes) | P3 θ_min≈44.9° (NBA); P5 ΔEAm 6/6>0 and Δr>0; P4 Ee favors BLA (localization) |
| Measured (LFMO-lite, RHF/STO-3G) | EV: +6.07→+0.43 kcal (0–30°), dEV/dθ = −0.17; Enσσ: +38.31→+5.48 kcal (0–45°), dEnσσ/dθ = −0.71; Eπσ(0°) = 0.00 (vs Table 5-16) |
| Rejected | θ=45° EV=+55.7 (π assignment unstable); Eπσ vs θ slope not isolated |
| Threshold met | **Yes** |
| G1–G5 | Pass |
| **Verdict** | **Agree** (~94%) |

### P3 — NBA large-twist crowded conformation most stable

| Item | Content |
|------|---------|
| Pre-registered criterion | Agree requires E_min ∈ [30°, 60°] and EN↑/Ee↓ in that region |
| Measured | θ_min = **44.9°**, ΔE = −0.96 kcal/mol; free B3LYP fold **≈34.5°** |
| Convergence proxy | Bidirectional scan hysteresis max **0.28** kcal/mol (threshold ≤1) |
| Threshold met | **Yes** |
| G1–G5 | Pass (**previously rejected**: wrong-ring data in `results/P3/invalid_wrong_geometry/`) |
| **Verdict** | **Agree** (~95%) |

### P4 — Benzene D₆ₕ driven by nuclear repulsion

| Item | Content |
|------|---------|
| Pre-registered criterion | Agree requires E_tot minimum near equal bond lengths, with EN favorable / Ee unfavorable for equalization |
| Measured (B3LYP/6-31G*) | Toward equal bonds: ΔEN = **−96**, ΔEe = **+88** kcal/mol; E_tot and EN minima at δ=0; Ee minimum at δ=**0.12 Å** |
| Threshold met | **Yes** (signs and extremum locations) |
| G1–G5 | Pass |
| **Verdict** | **Agree** (~90%) |

### P5 — Local ΔEAm always destabilizing; single-bond lengthening

| Item | Content |
|------|---------|
| Pre-registered criterion | Agree requires systematic ΔEAm>0 and intervening single-bond Δr>0 |
| Measured | ΔEAm **6/6 positive**; butadiene Δr = **+0.018 Å**; hexatriene Δr* = **+0.0045 Å** |
| Threshold met | **Yes** |
| G1–G5 | Pass |
| **Verdict** | **Agree** (~96%, v5c three objections closed) |

### P6 — Parameter-free ESE benchmarks

| Item | Content |
|------|---------|
| Pre-registered criterion | Agree requires benzene ESE<0 with \|ESE\| ≈ 30–40; CBD >0 and ~45–60; not overly level/basis sensitive |
| Measured | Benzene ESE = **−35.44** (book −36.3); CBD vertical@G* = **+53.98** (book +53.6); CBD 2D ΔEA = +65.51 |
| Threshold met | **Yes** |
| G1–G5 | Pass |
| **Verdict** | **Agree** (~98%, v4 three objections closed) |

### P7 — Strained-aromatic BLA attributed to π delocalization

| Item | Content |
|------|---------|
| Pre-registered criterion | Agree requires large positive Δr(G) and significant collapse after PLG (drop ≫ 0.05 Å) |
| Measured (C₁₂H₆) | Δr(G) = **+0.2065 Å** (book +0.177) → Δr(PLG) = **+0.020 Å** (book −0.002); drop = **0.1865 Å** |
| Side evidence | PLG state ~11.2 kcal/mol lower near Δr≈0 than at Δr(G) |
| Threshold met | **Yes** |
| G1–G5 | Pass |
| **Verdict** | **Agree** (~90%) |

### P8 — Furan-like systems use LDE protocol

| Item | Content |
|------|---------|
| Pre-registered criterion | Agree requires furan-like ΔEA>0 and LDE<0, sign split vs benzene (ΔEA<0), basis-stable |
| Measured | Furan ΔEA=+28.61, LDE=−37.72 (semi-adiabatic **−39.04**, book −39.3); pyrrole −48.20; oxazole ΔEA=+24.59, LDE=−41.17; benzene ΔEA=**−6.54**, ESE=−35.44 |
| Sensitivity | 6-31G* vs 6-31G: \|ΔLDE\| = **0.68** kcal/mol |
| Threshold met | **Yes** |
| G1–G5 | Pass |
| **Verdict** | **Agree** (~96%, v2 three objections closed) |

### P9 — Large [N]annulenes trend polyene-like

| Item | Content |
|------|---------|
| Pre-registered criterion | Agree requires clear 4n/4n+2 sign split for small N and converging gap as N grows |
| Measured (N=8–18) | Signs **6/6** match 4n+2; \|VDE/π\| gap **2.836 → 1.040 → 0.424** |
| Sensitivity | 2011-lite all signs correct (\|Δ\|≤0.012); BLA geometry N=8 +1.472 vs +1.619, no sign flip |
| Caveat | Proxy is 2007 ESE, not Yu EV(2011) absolute values; planar Kekulé; N=20–26 not scanned |
| Threshold met | **Yes** |
| G1–G5 | Pass (N=16 BLA abort isolated) |
| **Verdict** | **Agree** (~94%, v2b three objections closed) |

## 3. Sign and unit conventions

- Stabilization **negative**; destabilization **positive** (same as Yu)
- Internal Hartree; report kcal/mol; 1 Ha = 627.5095 kcal/mol
- Primary engine: WSL PySCF 2.14.0; no environment changes
