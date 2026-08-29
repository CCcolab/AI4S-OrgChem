# Release notes — v1.0.1-post-freeze-2026-08-28

Patch release aligning the **GitHub Release tarball/zipball** with `main` after the v1.0.0 freeze tag.  
**No change to scientific verdicts** (still **8 Agree + 1 Disagree on P1**).

| | |
|--|--|
| **Tag** | `v1.0.1-post-freeze-2026-08-28` |
| **Base verdict** | Same as `v1.0.0-freeze-2026-08-25` |
| **Prior tag** | `v1.0.0-freeze-2026-08-25` (`aa35800`) |

---

## AI4S-OrgChem v1.0.1 — Post-freeze maintenance (2026-08-28)

### Why this release

`v1.0.0-freeze-2026-08-25` pointed at an early public commit. **`main` moved ahead by 11 commits** (English docs, branding, privacy redaction). Downloaders who used the v1.0.0 Release asset saw an outdated tree (e.g. unredacted WSL hostnames). **v1.0.1 tags current `main`** so Release downloads match the default branch.

### What changed (documentation & hygiene only)

| Area | Change |
|------|--------|
| **Internationalization** | English default `README.md`; mirrors for freeze pack, unit `VERDICT.en.md`, `propositions.en.md`, `quality_gates.en.md`; expert quick-review guide |
| **Attribution / branding** | Book author **Professor Zhongheng Yu**; monograph title formatting; removed product-specific labels from public docs |
| **Privacy** | WSL `inventory/`, `MANIFEST.txt`, `README.md`, `mounts/`, `pyvenv.cfg` redacted (`wsluser` / `wsl-host`); `tools/wsl_snapshot.sh` uses `$HOME` + `redact_inventory()` |
| **`.gitignore`** | Ignore local-only `开发过程全记录.md` and `docs/security_architecture.md` |
| **Citation** | `CITATION.cff` contact email and repository URL finalized |

### What did **not** change

- P1–P9 **Agree / Disagree** outcomes  
- `results/Pn/tables/` primary evidence JSON/CSV used in verdicts  
- `deliverables/unit/Pn/VERDICT.md` scientific conclusions  
- PySCF scripts in `src/` (except `wsl_snapshot.sh` hygiene)

### Start here (unchanged science, improved docs)

| Doc | Link |
|-----|------|
| Homepage (EN) | [README.md](https://github.com/CCcolab/AI4S-OrgChem/blob/main/README.md) · [中文](https://github.com/CCcolab/AI4S-OrgChem/blob/main/README.zh-CN.md) |
| v1.0.0 freeze notes | [RELEASE_NOTES.md](https://github.com/CCcolab/AI4S-OrgChem/blob/main/docs/RELEASE_NOTES.md) |
| Executive summary | [EXECUTIVE_SUMMARY.en.md](https://github.com/CCcolab/AI4S-OrgChem/blob/main/deliverables/final/EXECUTIVE_SUMMARY.en.md) |
| Unit verdicts P1–P9 | [INDEX.en.md](https://github.com/CCcolab/AI4S-OrgChem/blob/main/deliverables/unit/INDEX.en.md) |
| QC expert guide | [expert_quick_review_guide.md](https://github.com/CCcolab/AI4S-OrgChem/blob/main/docs/expert_quick_review_guide.md) |

### Scoreboard (unchanged)

| Prop. | Verdict |
|-------|---------|
| **P1** | **Disagree** |
| P2–P9 | Agree |

### Copyright & contact

- **Book author (subject of verification):** Professor Zhongheng Yu  
- **Software:** Copyright (c) 2026 Xiao Chen · chenxiao0101@gmail.com  
- MIT License · See `LICENSE`, `NOTICE`, `CITATION.cff`

---

中文读者：科学结论与 v1.0.0 冻结版相同；本补丁仅同步文档、隐私脱敏与 `.gitignore`。请从 [README.zh-CN.md](https://github.com/CCcolab/AI4S-OrgChem/blob/main/README.zh-CN.md) 进入。
