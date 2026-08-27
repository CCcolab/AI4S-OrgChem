# Proposition specifications (P1–P9 · frozen)

**English** | [中文](propositions.md)

> **Project**: **AI-for-Science (AI4S) Agent** using state-of-the-art AI; third-party verification with WSL/PySCF (no changes to the existing environment) of claims in **Professor Zhongheng Yu**’s **Questioning Fundamental Principles of Organic Chemistry** that challenge traditional organic structural theory.  
> **Freeze (2026-08-25)**: **8 Agree + 1 Disagree (P1)**. Formal verdicts: `deliverables/unit/Pn/VERDICT.md` and `deliverables/final/`.  
> **Delivery**: each proposition gets **Agree** or **Disagree**; conclusions must withstand expert review.  
> **Compliance**: published math, orbital diagrams, energy expressions, and numerical tables only—no use/copy/decompile/translate of the author’s programs. Book PDF / extract **not** in the repo.  
> **Tracking**: [`implementation/`](../implementation/). Expert guide: [`expert_quick_review_guide.md`](expert_quick_review_guide.md).

---

## 0. Completeness of the proposition set

After cross-checking public sources, **P1–P9 are the complete set of propositions that must receive verdicts**:

| Source | Role |
|--------|------|
| [OSF project abstract](https://doi.org/10.17605/OSF.IO/CAKFZ) | Lists core heterodox claims and three major benchmarks by chapter |
| ScienceNet / monograph preface (author’s blog) | Author’s starting points and method evolution |
| Amazon / Academia Proof abstracts | Consistent with OSF: butadiene, NBA, ESE, benzene D₆ₕ, furan, [N]annulene, strained aromatics |

**Book → proposition map**:

| Book point (public abstracts) | This project |
|-------------------------------|--------------|
| Ch1 amino-polyene aldehydes: conjugation causes distortion | Folded into **P2/P3** evidence (not a separate delivery) |
| Ch2 NBA large twist; crowded conformation most stable | **P3** |
| Ch3–5 LFMO: π–π / σ–σ / π–σ all destabilizing; role split | **P2** mechanism subtask + **P3** |
| Ch6/8/12 butadiene CE destabilizing; artificial reference choice | **P1** |
| Local double-bond-pair conjugation additive (polyenes) / not (aromatics) | **P5** (additivity) → **P6** (ESE) |
| Three benchmarks: benzene −36.3; CBD +53.6 | **P6** |
| Benzene D₆ₕ by nuclear-repulsion minimization; π favors D₃ₕ | **P4** |
| Strained-aromatic BLA from conjugation, not angle strain | **P7** |
| Furan-like between aromatic/nonaromatic (LDE) | **P8** |
| Large [N]annulene (N>18) trends polyene-like | **P9** |
| PBH / substituted butadienes etc. | **Extended checks** (not separate fundamental-principle deliveries) |

**Why not separate**: PBH and substituted butadienes are method extrapolations; Ch1 and LFMO support P2/P3 and are executed as subtasks in those plans to avoid double-counting “fundamental principle” votes.

---

## 1. Uniform verdict rules

For each proposition `P`, after independent computation, one verdict:

| Verdict | Meaning |
|---------|---------|
| **Agree** | Under the agreed level, basis, and reference definitions, independent results support Yu’s claim in sign, magnitude, and key numbers, and withstand expert challenges on reference/method choice |
| **Disagree** | Independent results conflict in sign, magnitude, or key conclusion, or the claim is not reproducible under its own definitions / rests on indefensible conventions |

**Hard requirements (defensibility)**:
1. Each proposition must be a **falsifiable statement**: variables, observables, thresholds, controls.  
2. Report both **same-protocol results as Yu** and **reference/method sensitivity** (if sensitivity flips the conclusion, do not lightly Agree unless Yu claims that reference choice itself).  
3. Numerical near-coincidence with Yu is not sufficient; chemical meaning must hold.  
4. Final verdicts only **Agree / Disagree**; if evidence is insufficient, remain “in progress,” do not early-judge.  
5. **Hard quality gates (mandatory)** before `VERDICT.md`: topology/correct molecule, geometry integrity, convergence (or hysteresis proxy), energy scale, voided-data isolation. **Never** verdict on “unconverged + wrong molecule.” Spec: [`quality_gates.en.md`](quality_gates.en.md).

**Sign convention** (same as Yu): stabilization negative; destabilization positive. Report unit: kcal/mol.

---

## 2. Master table

| ID | One-line claim (Yu vs tradition) | Traditional foil | Tools | Difficulty | Verdict | Status |
|----|----------------------------------|------------------|-------|------------|---------|--------|
| **P1** | Classical butadiene “conjugation stabilization” evidence is unreliable; vs trans-2-butene CE is destabilizing | CE ≈ −3.9 kcal/mol vs two 1-butenes proves conjugation | Thermochemical cycles; CE₁ vs CE₂; RHF/B3LYP/MP2/ZPE/6-31G*; homemade 2007-GL | Low | **Disagree** | Closed (~98%) |
| **P2** | Conjugation (π delocalization) is often **destabilizing** and a **distortion driver**, not the main stabilizer of coplanarity | Conjugation stabilization makes coplanar most stable | Aggregate P3/P5/P6; LFMO/π–σ (from published defs); PES | High | **Agree** | Closed (~94%, v2 + LFMO-lite) |
| **P3** | Most sterically crowded conformation can be the energy minimum (NBA-type large twist) | Steric destabilization → crowded conformations unstable | Cascaded pre-relax + bidirectional θ constraints + B3LYP SP; free B3LYP check | Mid | **Agree** | Closed (~95%) |
| **P4** | Benzene D₆ₕ mainly from nuclear-repulsion minimization; π energy favors bond alternation (D₃ₕ) | π delocalization drives equal bonds and aromatic stability | BLA scan; Ee/EN/E_tot; optional localized GL | Mid | **Agree** | Closed (~90%) |
| **P5** | Local pair conjugation ΔEAm always destabilizing and lengthens the intervening single bond | Local conjugation stabilizes and shortens the single bond | GL/GE-m localized geometries; ΔEAm & Δr stats; additivity | High | **Agree** | Closed (~96%, v5c) |
| **P6** | Parameter-free ESE: benzene ≈ −36.3 (stable); CBD conjugation/antiaromaticity ≈ +53–55 (destabilizing) | ASE depends on empirical/semi-empirical params or isodesmic choice | Localized geom + ΣΔEAm + VR; ESE benchmarks; level/basis sensitivity | Very high | **Agree** | Closed (~98%, v4) |
| **P7** | Strained-aromatic central-ring BLA from π delocalization, not angle strain (Mills–Nixon) | Mills–Nixon: angle strain causes central BLA | G vs PLG; Δr = r_endo−r_exo; cut center–periphery π coupling | High | **Agree** | Closed (~90%, v1c) |
| **P8** | Furan-like systems should not be scored aromatic via benzene-like ESE; use LDE; intermediate state | Furan often listed as aromatic heterocycle | ΔEA / ΣΔEAm sign patterns; furan·pyrrole vs benzene; LDE definition | Mid–high | **Agree** | Closed (~96%, v2) |
| **P9** | Large [N]annulenes (N≳16–18) trend polyene-like in delocalization energy; aromaticity/antiaromaticity weaken | 4n+2 still substantively aromatic/antiaromatic for large rings | Series VDE (or equivalent) vs N; 4n vs 4n+2 gap convergence | Mid–high | **Agree** | Closed (~94%, v2b) |

> Draft C1–C7 were reorganized into **P1–P9**. Implementation sync: [`implementation/命题总表.md`](../implementation/命题总表.md).

---

## 3. Detailed specifications

### P1 — Butadiene CE sign depends on reference molecule

**Yu**: Literature CE ≈ −3.9 vs two 1-butenes “proves” conjugation; vs **trans-2-butene**, CE ≈ **+1.9 (destabilizing)**. GL ≈ +1.4. If classical evidence stands, the book framework is empty.

**Tradition**: Conjugation stabilization is fundamental; short butadiene single bond and negative hydrogenation difference are classic evidence.

**Falsifiable test**:
1. Same level (start B3LYP/6-31G*; cross RHF/MP2):  
   - CE₁ = ΔH_hyd(butadiene) − 2·ΔH_hyd(1-butene)  
   - CE₂ = ΔH_hyd(butadiene) − 2·ΔH_hyd(trans-2-butene)  
2. **Agree**: CE₁ < 0 and CE₂ > 0 (flip); CE₂ in +1–3 kcal/mol.  
3. **Disagree**: CE₂ still negative, or flip not reproducible under a defensible thermochemical definition.

**Tools**: PySCF opt + reaction energies; ZPE; homemade 2007-GL (`src/localization/`).  
**Verdict**: L1/L2 **Disagree** (no flip); L3 vs book +1.4 **Disagree (magnitude)** (homemade ΔE≈+4.06).

---

### P2 — π delocalization destabilizes and drives distortion (umbrella)

**Yu**: π–π and nonbonded σ–σ destabilize and drive distortion; π–σ also destabilizes but resists twist; large twist is a compromise. Conjugation-as-fundamental-principle should be questioned.

**Tradition**: Conjugation makes coplanar most stable; conformation balances conjugation vs sterics.

**Falsifiable test** (aggregate; depends on P3/P5/P6):
1. Under NBA-type or localized-geometry frameworks, key interaction signs match Yu (π–π, σ–σ destabilizing; distortion driver correctly attributed).  
2. **Agree**: at least two key classes (NBA-type + polyene/aromatic localization) reproduce “π delocalization destabilizes / drives distortion.”  
3. **Disagree**: mainstream defensible decompositions systematically reverse signs or mis-attribute the driver.

---

### P3 — NBA-type: large-twist crowded conformation most stable

**Yu**: Experimental twist ~36°–55°; on relaxed PES, energy lowest near θ≈θ_exp with high EN / low Ee—**crowded conformation most stable**.

**Tradition**: Steric destabilization → large-twist/crowded forms rise in energy.

**Falsifiable test**:
1. Relaxed PES in θ for Ph–N=CH–Ph (and 1–2 hetero analogs).  
2. **Agree**: E(θ) minimum in large-angle region (~30°–60°, crystallographic scale); EN up / Ee down vs planar nearby.  
3. **Disagree**: global min near planar (θ≲15°), or large-angle only a much higher secondary min.

**Tools**: PySCF relaxed scan; Ee/EN; **topology + quality gates**.

---

### P4 — Benzene D₆ₕ: nuclear repulsion vs π distortion preference

**Yu**: From localized GL → delocalized ground state, ΔEe > 0, ΔEN < 0 with |ΔEN| > ΔEe; equal bond lengths minimize nuclear repulsion → D₆ₕ.

**Tradition**: π delocalization drives equal bonds and aromatic stability.

**Falsifiable test**:
1. Scan benzene along BLA; record E, Ee, EN.  
2. **Agree**: total-energy minimum near equal bonds; along equalization, nuclear-repulsion contribution dominates stabilization (or at least matches Yu qualitatively: Ee unfavorable, EN favorable).  
3. **Disagree**: Ee change dominates stabilization, or equal-bond geometry is not an EN minimum.

**Verdict**: this tier **Agree** (BLA + Ee/EN; full GL→G separate).

---

### P5 — Local ΔEAm always destabilizing; single-bond lengthening

**Yu**: For many double-bond pairs, ΔEAm = E(GE-m)−E(GL) > 0; intervening single bond longer in GE-m than GL; r vs ΔEAm fit possible.

**Tradition**: Local conjugation stabilizes (lowers energy) and shortens the linking single bond.

**Falsifiable test**:
1. Homemade GL/GE-m (published matrix-deletion defs) on butadiene, hexatriene, benzene, etc.: ΔEAm and Δr.  
2. **Agree**: sample ΔEAm predominantly positive; corresponding Δr > 0.  
3. **Disagree**: under a reasonable localization definition, ΔEAm systematically negative or single bonds systematically shorten.

---

### P6 — Parameter-free ESE benchmarks: benzene and cyclobutadiene

**Yu**:  
ESE = ΔEA − ΣΔEAm = E(G) − E(VR).  
Benzene ESE ≈ **−36.3** kcal/mol; CBD destabilizing conjugation/antiaromaticity ≈ **+53–55**; matching experimental sign/scale is necessary for method reasonableness.

**Tradition**: ASE usually depends on isodesmic/homodesmic reactions or empirical parameters.

**Falsifiable test**:
1. Homemade localization + ESE per published 2011/2014-facing definitions (as implemented).  
2. **Agree**: benzene ESE < 0 with \|ESE\| ~30–40; CBD corresponding quantity > 0 and ~45–60; not overly level/basis sensitive.  
3. **Disagree**: wrong sign, or magnitude far from experiment/Yu without a defensible implementation difference.

---

### P7 — Strained aromatics: BLA from π delocalization, not angle strain

**Yu**: In PLG (exclude π interaction between central ring and fused small rings), Δr falls from large positive in the ground state to near zero → BLA attributed to π delocalization.

**Tradition**: Mills–Nixon—angle strain causes central BLA.

**Falsifiable test**:
1. Compare G vs PLG Δr = r_endo − r_exo for benzotricyclobutadiene (or standard models).  
2. **Agree**: after cutting center–periphery π coupling, \|Δr\| drops sharply (→ 0).  
3. **Disagree**: Δr essentially unchanged after cutting π → supports angle-strain dominance.

---

### P8 — Furan-like: LDE not ESE; between aromatic and nonaromatic

**Yu**: Furan-like ΔEA > 0 and ΔEA < ΣΔEAm, so (ΔEA−ΣΔEAm)<0 but should not be called ESE; call it LDE; not directly comparable to benzene-like cases.

**Tradition**: Furan often treated as an aromatic heterocycle.

**Falsifiable test**:
1. Compute ΔEA, ΣΔEAm, and difference sign patterns for furan/pyrrole etc.  
2. **Agree**: reproduce “ΔEA>0 yet difference still negative,” qualitatively distinct from benzene (ΔEA<0).  
3. **Disagree**: furan and benzene fall in the same sign pattern; no LDE/intermediate distinction.

---

### P9 — Large [N]annulenes trend polyene-like

**Yu**: From N=8 to 26, VDE still respects 4n+2; but for N≳16, [4n] and [4n+2] VDE/π approach each other—large rings more like polyenes.

**Tradition**: Large rings should still show strong aromatic/antiaromatic character.

**Falsifiable test**:
1. Series of [N]annulene VDE (or equivalent delocalization metrics) vs N.  
2. **Agree**: small rings show 4n+2 sign separation; at large N, \|VDE/π\| gap converges.  
3. **Disagree**: large N keeps a strong, non-converging aromatic/antiaromatic split.

---

## 4. Recommended verification order (defensibility first)

1. **P1** — most independent; faces “is the book built on cherry-picked references?”  
2. **P3** — experimental conformational fact; no localization method required  
3. **P4** — Ee/EN; aromatic geometry controversy  
4. **P6 + P5** — method core and benchmarks (most expert pushback)  
5. **P7** — strained-aromatic attribution  
6. **P8, P9** — extension and classification  
7. **P2** — umbrella verdict (depends on the above)

---

## 5. Old → new ID map

| Old | New |
|-----|-----|
| C1 | P1 |
| C2 | P3 |
| C3 | P4 |
| C4 | P5 |
| C5 | P6 |
| C6 | P7 |
| C7 | P8 + P9 |
| (umbrella) | P2 |
