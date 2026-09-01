# Verdict authority tiers (2026-09-01 post-review remediation)

**English** | [中文](FROZEN_VERDICT_AUTHORITY.md)

> **Basis**: Post-remediation QC review dated 2026-09-01 (§4 AI4S-OrgChem standalone review, §8 P0)  
> **Goal**: Public-facing **five-state + estimand layering**; **do not** use “8 Agree · 1 Disagree” as a scientific score.  
> **L1 pre-registered threshold audit (2026-08-25 snapshot)** remains historical audit only—not the public scientific conclusion.

---

## 1. Authority sources (priority)

| Tier | File / field | Force |
|------|--------------|-------|
| **L0 (public primary)** | **Public status (L0)** five-state label in each `deliverables/unit/Pn/VERDICT.md` | **Public scientific status**; README, summaries, expert guide **must** follow this |
| **L1 (audit snapshot)** | **Pre-registered threshold audit (L1)** in the same file | **2026-08-25** binary threshold met or not—**not** a scientific score |
| **L2** | `deliverables/final/VERDICT_TABLE.md` | L0 + L1 crosswalk for all nine units |
| **L3** | `EXECUTIVE_SUMMARY` / `FULL_REPORT` | Narrative; must not override L0 |
| **L4** | README, chat, Canvas | Non-verdict; cite L0–L2 only |

**Rule**: Disputes about *how strongly a claim is supported* → **L0 five-state + estimand table only**; L1 answers only whether the frozen pre-registered threshold was met.

---

## 2. Public scientific status matrix (L0 · primary)

**Independent verification count = 8** (P2 is `DERIVED`, **excluded**).  
**Forbidden** publicly: “8/9 book propositions verified” or “8 Agree / 1 Disagree = eight scientific truths.”

| Prop. | **L0 public status** | Estimand / scope | **L1 audit snapshot** (2026-08-25) |
|-------|----------------------|------------------|-------------------------------------|
| **P1** | **`PARTIAL`** (composite) | P1a experimental ΔH298 `INCOMPARABLE`; P1b QC CE `OPPOSED`; P1c GL2007 vs 2014 `INCOMPARABLE` | P1b CE flip criterion → **Disagree** |
| **P2** | **`DERIVED`** | Meta aggregate P1/P3–P9 + LFMO-lite; **no independent SCF** | Aggregate threshold met → **Agree** |
| **P3** | `SUPPORTED_WITHIN_SCOPE` | NBA tested path θ_min≈44.9° | **Agree** |
| **P4** | `PARTIAL` | Frozen BLA path signs pass; **causal claim not tested** | **Agree** |
| **P5** | `SUPPORTED_WITHIN_SCOPE` | Six tested instances ΔEAm>0; not “always universal” | **Agree** |
| **P6** | `SUPPORTED_WITHIN_SCOPE` | Benzene ESE≈−35.44; CBD ≈+53.98 | **Agree** |
| **P7** | `PARTIAL` | C₁₂H₆ single system; BLA collapse under PLG | **Agree** |
| **P8** | `SUPPORTED_WITHIN_SCOPE` | Furan LDE vs benzene ESE sign split | **Agree** |
| **P9** | **`PARTIAL`** | N=8–18 planar Kekulé trend; **O1 open** (N≥20, nonplanar TBD) | N=8–18 sub-criterion met → **Agree** |

### L0 five-state definitions

| State | Meaning | Public “support” wording |
|-------|---------|--------------------------|
| `SUPPORTED_WITHIN_SCOPE` | Same estimand, pre-registered protocol passes | Yes, **with explicit scope** |
| `PARTIAL` | Sub-claim / single system / trend only, or key objection open | “Partial support” only |
| `OPPOSED` | Pre-registered threshold not met under same estimand | Unsupported for that estimand |
| `INDETERMINATE` | Insufficient evidence, convergence failure, missing methods | No strong verdict |
| `INCOMPARABLE` | Different estimand, reference state, or sign convention | No voting, no averaging |
| `DERIVED` | Aggregate only; no new QM data | **Do not** count as independent verification |

---

## 3. P1 must use estimand as primary key

| ID | Estimand | L0 | Note |
|----|----------|-----|------|
| P1a | Experimental ΔH298 hydrogenation differences | `INCOMPARABLE` | Not primary machine delivery; must not be covered by P1b |
| P1b | B3LYP/RHF/MP2/ZPE CE₁/CE₂ | `OPPOSED` | No CE sign flip |
| P1c | GL(2007) +4.06 vs book +1.4 (2014) | `INCOMPARABLE` | Not comparable to GL(2014) |

**P1 overall L0 = `PARTIAL`**: L1 Disagree on P1b must **not** summarize the entire book composite claim.

---

## 4. P9 O1 reopened (2026-09-01 correction)

| Objection | Status | Note |
|-----------|--------|------|
| O1 extend N / exact onset | **Open** | N=20–26 not scanned; N=16 BLA aborted; nonplanar excluded |
| O2 2011-lite | Closed | Signs match 2007 |
| O3 Kekulé bond lengths | Closed | BLA 1.34/1.46 signs unchanged |

**P9 L0 = `PARTIAL`**: Trend layer informative; **cannot** claim full large-ring verification or O1 closed.

---

## 5. What is **not** a verdict reversal

1. **L0 five-state updates** (§2): governance response to third-party review—not new QM reruns.  
2. **L1 snapshot retained**: 2026-08-25 binary audit kept as historical record only.  
3. **Wording / metadata / validator upgrades**: engineering hardening.  
4. **P2 marked `DERIVED`**: identity clarification only.

**True L0 change** requires: new evidence + pre-registration + G1–G5 + update `VERDICT.md` + sync L2/L3 + git audit note.

---

## 6. Quality gates

- G1–G5: hard gates before L1 audit (see [`quality_gates.en.md`](quality_gates.en.md)).  
- `tools/validate_repo.py`: checks L0 presence, JSON schema, no duplicate keys, key `quality_gate`—**does not** treat L1 binary snapshot as public scientific score.

---

## 7. Maintenance

| Item | Convention |
|------|------------|
| Unit VERDICT | Header must include **L0** and **L1** lines |
| README | Show L0 matrix; **forbid** “8/1 scientific score” badge |
| Index | [`deliverables/unit/INDEX.en.md`](../deliverables/unit/INDEX.en.md) |

*2026-08-31 initial · 2026-09-01 post-review remediation*
