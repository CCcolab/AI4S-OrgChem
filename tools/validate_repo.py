#!/usr/bin/env python3
"""Repository consistency checks (post-2026-09-01 remediation). Exit 0 if pass, 1 if fail."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

try:
    import jsonschema
except ImportError:  # pragma: no cover
    jsonschema = None

ROOT = Path(__file__).resolve().parents[1]
AUTHORITY = ROOT / "docs" / "FROZEN_VERDICT_AUTHORITY.md"
SCHEMA_PATH = ROOT / "schemas" / "evidence.schema.json"

L0_EXPECTED = {
    "P1": "PARTIAL",
    "P2": "DERIVED",
    "P3": "SUPPORTED_WITHIN_SCOPE",
    "P4": "PARTIAL",
    "P5": "SUPPORTED_WITHIN_SCOPE",
    "P6": "SUPPORTED_WITHIN_SCOPE",
    "P7": "PARTIAL",
    "P8": "SUPPORTED_WITHIN_SCOPE",
    "P9": "PARTIAL",
}

L1_AUDIT_SNIPPETS = {
    "P1": "非一致",
    "P2": "一致",
    "P3": "一致",
    "P4": "一致",
    "P5": "一致",
    "P6": "一致",
    "P7": "一致",
    "P8": "一致",
    "P9": "一致",
}

L0_LINE = re.compile(
    r"- \*\*(?:对外状态|Public status)[（(]L0[）)]\*\*[：:]\s*`([A-Z_]+)`"
)
L1_LINE = re.compile(
    r"- \*\*(?:预注册阈值审计|Threshold audit)[（(]L1[）)]\*\*[：:].*\*\*(非一致|一致|Disagree|Agree)\*\*"
)

KEY_EVIDENCE = [
    "results/P1/tables/ce_B3LYP_6-31gs.json",
    "results/P1/tables/ce_zpe_B3LYP_6-31gs.json",
    "results/P1/tables/gl2007_butadiene_B3LYP_6-31gs.json",
    "results/P5/tables/p5_v5c_objections_B3LYP_6-31gs.json",
    "results/P9/tables/p9_v2b.json",
    "results/P2/tables/p2_v2_aggregate.json",
]


def fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)


def ok(msg: str) -> None:
    print(f"OK: {msg}")


def load_json_strict(path: Path) -> tuple[dict | None, list[str]]:
    """Parse JSON and report duplicate keys at any object depth."""
    text = path.read_text(encoding="utf-8")
    dups: list[str] = []

    def hook(pairs: list[tuple[str, object]]) -> dict:
        seen: set[str] = set()
        out: dict = {}
        for key, val in pairs:
            if key in seen:
                dups.append(key)
            seen.add(key)
            out[key] = val
        return out

    try:
        data = json.loads(text, object_pairs_hook=hook)
    except json.JSONDecodeError as e:
        return None, [f"JSONDecodeError: {e}"]
    return data, dups


def check_authority_doc() -> bool:
    if not AUTHORITY.is_file():
        fail("missing docs/FROZEN_VERDICT_AUTHORITY.md")
        return False
    text = AUTHORITY.read_text(encoding="utf-8")
    if "L0（对外主键）" not in text or "DERIVED" not in text:
        fail("FROZEN_VERDICT_AUTHORITY.md missing L0/DERIVED remediation")
        return False
    if "8 一致 · 1 非一致" in text and "禁止" not in text:
        fail("authority doc still promotes 8/1 without forbid note")
        return False
    ok("authority doc (L0-primary) present")
    return True


def check_l0_l1_verdicts() -> bool:
    good = True
    for pid, expected_l0 in L0_EXPECTED.items():
        vpath = ROOT / "deliverables" / "unit" / pid / "VERDICT.md"
        if not vpath.is_file():
            fail(f"{pid} missing VERDICT.md")
            good = False
            continue
        text = vpath.read_text(encoding="utf-8")
        m0 = L0_LINE.search(text)
        if not m0 or m0.group(1) != expected_l0:
            fail(f"{pid} L0 expected `{expected_l0}`, got {m0.group(1) if m0 else None}")
            good = False
        else:
            ok(f"{pid} L0 = {expected_l0}")
        m1 = L1_LINE.search(text)
        snippet = L1_AUDIT_SNIPPETS[pid]
        if not m1 or snippet not in m1.group(0):
            fail(f"{pid} L1 audit missing `{snippet}`")
            good = False
        else:
            ok(f"{pid} L1 audit contains {snippet}")
    return good


def check_all_results_json() -> bool:
    good = True
    schema = None
    if SCHEMA_PATH.is_file() and jsonschema is not None:
        schema = jsonschema.Draft202012Validator(
            json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
        )

    tables = sorted((ROOT / "results").glob("P*/tables/*.json"))
    invalid_paths = sorted((ROOT / "results").glob("P*/invalid_*/*.json"))
    all_paths = tables + invalid_paths
    schema_targets = {ROOT / rel for rel in KEY_EVIDENCE}

    for p in all_paths:
        rel = p.relative_to(ROOT).as_posix()
        data, dups = load_json_strict(p)
        if data is None:
            fail(f"invalid JSON {rel}: {dups[0]}")
            good = False
            continue
        if dups:
            fail(f"{rel}: duplicate keys {sorted(set(dups))}")
            good = False
        if "invalid_" not in rel and rel.startswith("results/P") and rel.endswith(".json"):
            audit_only = rel.endswith(
                ("p1_quality_audit.json", "sensitivity_methods.json")
            )
            if not audit_only:
                qg = data.get("quality_gate")
                if not isinstance(qg, dict) or "passed" not in qg:
                    fail(f"{rel}: missing quality_gate.passed")
                    good = False
        if schema is not None and p in schema_targets:
            errors = sorted(schema.iter_errors(data), key=lambda e: e.path)
            if errors:
                fail(f"{rel}: schema {errors[0].message}")
                good = False

    ok(f"checked {len(all_paths)} results JSON files (strict parse)")
    return good


def check_key_evidence() -> bool:
    good = True
    for rel in KEY_EVIDENCE:
        p = ROOT / rel
        if not p.is_file():
            fail(f"missing key evidence {rel}")
            good = False
            continue
        data, dups = load_json_strict(p)
        if dups:
            fail(f"{rel}: duplicate keys {sorted(set(dups))}")
            good = False
        if rel.endswith("p9_v2b.json"):
            o1 = (data or {}).get("objections", {}).get("O1_extend_N", {})
            if o1.get("closed") is not False:
                fail("p9_v2b.json O1_extend_N must be closed=false")
                good = False
            else:
                ok("p9_v2b O1 open")
        ok(f"{rel} key evidence OK")
    return good


def check_forbidden_public_score() -> bool:
    """README must not use 8/1 as scientific score badge."""
    good = True
    for name in ("README.md", "README.zh-CN.md"):
        p = ROOT / name
        if not p.is_file():
            continue
        text = p.read_text(encoding="utf-8")
        if re.search(r"verdict-8[_\s].*1_(非一致|Disagree)", text, re.I):
            fail(f"{name} still uses 8/1 verdict badge")
            good = False
        if "8 一致 · 1 非一致" in text and "科学总分" not in text and "不是" not in text:
            # allow if explicitly negated nearby - check first 120 lines
            head = text[:4000]
            if "8 一致 · 1 非一致" in head and "禁止" not in head and "不是" not in head:
                fail(f"{name} promotes 8/1 without forbid context in header")
                good = False
    if good:
        ok("README avoids 8/1 scientific-score badge")
    return good


def main() -> int:
    checks = [
        check_authority_doc(),
        check_l0_l1_verdicts(),
        check_all_results_json(),
        check_key_evidence(),
        check_forbidden_public_score(),
    ]
    if all(checks):
        print("\nvalidate_repo: ALL PASSED")
        return 0
    print("\nvalidate_repo: FAILED", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
