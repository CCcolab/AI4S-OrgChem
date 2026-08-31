# Frozen Verdict Authority (Zero Ambiguity)

[中文](FROZEN_VERDICT_AUTHORITY.md) | **English**

> **Freeze date**: 2026-08-25  
> **Frozen tally**: **8 Agree · 1 Disagree** (P1 Disagree; P2–P9 Agree)  
> **Purpose**: Define what counts as a formal verdict, what does not, and when a verdict may change—so wording edits, third-party review responses, or reading-aid labels are not mistaken for reversals.

---

## 1. Sole authoritative sources (priority order)

| Tier | File | Effect |
|------|------|--------|
| **L1 (highest)** | Each `deliverables/unit/Pn/VERDICT.md` | **Sole formal verdict** per proposition; values are **Agree** or **Disagree** only |
| **L2** | `deliverables/final/VERDICT_TABLE.en.md` | Audit table; must match L1 one-to-one |
| **L3** | `EXECUTIVE_SUMMARY.en.md`, `FULL_REPORT.en.md` | Narrative; **must not override L1** |
| **L4** | README, Canvas, chat, review responses | **Not verdicts**; may cite L1–L2 only |

Disputes about “what was decided” are resolved **only from L1**.

---

## 2. Frozen scoreboard (unchanged by wording edits)

| Prop. | **Frozen verdict** | Preregistered gist |
|-------|-------------------|----------------------|
| **P1** | **Disagree** | CE₁<0 and CE₂>0 not both met |
| **P2** | **Agree** | Two-class aggregate (meta-proposition) |
| **P3** | **Agree** | E_min in [30°,60°]; coplanar not global minimum |
| **P4** | **Agree** | On frozen BLA path; E_tot/EN min at δ=0; signs pass |
| **P5** | **Agree** | ΔEAm>0 in tested set; Δr>0 |
| **P6** | **Agree** | Benzene ESE≈−35.44; CBD vertical ΔEA≈+53.98 |
| **P7** | **Agree** | C₁₂H₆: BLA collapse under PLG |
| **P8** | **Agree** | Furan-type vs benzene sign split |
| **P9** | **Agree** | N=8–18: 4n/4n+2 signs; gap trend passes threshold |

Phrases like “trend support” or “partial support” elsewhere are **scope qualifiers (reading aids)** mapping to this table—**not new verdicts**.

---

## 3. What does **not** change the frozen verdict

1. Softer or stricter **wording** (P4 causal language, P9 trend layer, etc.)  
2. **Five-state reading-aid labels** (§4)—parallel to L1, never replace it  
3. **Metadata** / schema / provenance fixes  
4. **Third-party review responses** without §5 procedure  
5. **README / badge / Agent narrative** edits  
6. **P2 meta-proposition clarification**—P2 remains Agree  

---

## 4. Five-state reading aids (not formal verdicts)

| Aid | Meaning | Maps to frozen L1 |
|-----|---------|-------------------|
| `SUPPORTED_WITHIN_SCOPE` | Passes preregistered threshold in scope | **Agree** |
| `PARTIAL` | Passes threshold; limited extrapolation | Still **Agree** |
| `OPPOSED` | Threshold not met | **Disagree** |
| `INDETERMINATE` | Insufficient evidence | Not used at freeze |
| `INCOMPARABLE` | Different estimand | Sub-note only (e.g. P1c) |

---

## 5. Sole lawful verdict-change procedure

All required: new evidence → preregistered criteria → G1–G5 pass → update L1 `VERDICT.md` → sync `VERDICT_TABLE` → audit trail.  
**Forbidden**: changing L1 via chat, review text, reading aids, or README alone.

---

## 6. Three-way caliber (anti-misread; compatible with L1)

Local/pairwise destabilization (ΔEAm) · benzene extra stabilization (ESE) · geometric distortion (P3, P7) are **distinct**. This project is **not** “all conjugation is destabilizing.”

---

## 7. Relation to quality gates

G1–G5 were satisfied at freeze. Future G6–G11 strengthen engineering; they do **not** auto-overturn 2026-08-25 L1 without §5.

---

*Established 2026-08-31 · Zero-ambiguity frozen verdict policy*
