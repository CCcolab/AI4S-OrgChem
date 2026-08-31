#!/usr/bin/env python3
"""Repository consistency checks (freeze-era). Exit 0 if pass, 1 if failures."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FROZEN = {
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
VERDICT_LINE = re.compile(r"- \*\*判定[：:]\s*(一致|非一致)\*\*")
AUTHORITY = ROOT / "docs" / "FROZEN_VERDICT_AUTHORITY.md"


def fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)


def ok(msg: str) -> None:
    print(f"OK: {msg}")


def parse_verdict(path: Path) -> str | None:
    if not path.is_file():
        return None
    for line in path.read_text(encoding="utf-8").splitlines():
        m = VERDICT_LINE.search(line)
        if m:
            return m.group(1)
    return None


def check_frozen_verdicts() -> bool:
    good = True
    for pid, expected in FROZEN.items():
        vpath = ROOT / "deliverables" / "unit" / pid / "VERDICT.md"
        got = parse_verdict(vpath)
        if got != expected:
            fail(f"{pid} VERDICT.md: expected {expected}, got {got}")
            good = False
        else:
            ok(f"{pid} frozen L1 = {expected}")
    return good


def check_authority_doc() -> bool:
    if not AUTHORITY.is_file():
        fail("missing docs/FROZEN_VERDICT_AUTHORITY.md")
        return False
    text = AUTHORITY.read_text(encoding="utf-8")
    if "8 一致 · 1 非一致" not in text:
        fail("FROZEN_VERDICT_AUTHORITY.md missing frozen tally")
        return False
    ok("frozen authority doc present")
    return True


def check_json_parse() -> bool:
    good = True
    tables = list((ROOT / "results").glob("P*/tables/*.json"))
    for p in tables:
        try:
            json.loads(p.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            fail(f"invalid JSON {p.relative_to(ROOT)}: {e}")
            good = False
    ok(f"parsed {len(tables)} results JSON files")
    return good


def check_evidence_quality_gate() -> bool:
    """Warn-level: key evidence files should expose quality_gate."""
    required = [
        "results/P1/tables/ce_B3LYP_6-31gs.json",
        "results/P1/tables/ce_zpe_B3LYP_6-31gs.json",
        "results/P1/tables/gl2007_butadiene_B3LYP_6-31gs.json",
        "results/P5/tables/p5_v5c_objections_B3LYP_6-31gs.json",
        "results/P9/tables/p9_v2b.json",
        "results/P2/tables/p2_v2_aggregate.json",
    ]
    good = True
    for rel in required:
        p = ROOT / rel
        if not p.is_file():
            fail(f"missing evidence {rel}")
            good = False
            continue
        data = json.loads(p.read_text(encoding="utf-8"))
        qg = data.get("quality_gate")
        if not isinstance(qg, dict) or "passed" not in qg:
            fail(f"{rel}: missing quality_gate.passed")
            good = False
        else:
            ok(f"{rel} has quality_gate")
    return good


def main() -> int:
    checks = [
        check_authority_doc(),
        check_frozen_verdicts(),
        check_json_parse(),
        check_evidence_quality_gate(),
    ]
    if all(checks):
        print("\nvalidate_repo: ALL PASSED")
        return 0
    print("\nvalidate_repo: FAILED", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
