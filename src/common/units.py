"""Common constants and I/O helpers for AI4S-OrgChem."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HARTREE_TO_KCAL = 627.5094740631
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def ha_to_kcal(e_ha: float) -> float:
    return float(e_ha) * HARTREE_TO_KCAL


def kcal_to_ha(e_kcal: float) -> float:
    return float(e_kcal) / HARTREE_TO_KCAL


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, data: Any) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))
