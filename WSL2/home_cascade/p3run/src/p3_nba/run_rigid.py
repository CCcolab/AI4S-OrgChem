"""Rigid (non-relaxed) twist scan — fast tier for P3."""
from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pyscf import dft, gto, scf  # noqa: E402

from src.common.units import HARTREE_TO_KCAL, ensure_dir, ha_to_kcal, write_json  # noqa: E402
from src.p3_nba.geometry import (  # noqa: E402
    IDX_C_IMINE,
    IDX_CIPSO_N,
    IDX_CORTHO_N,
    IDX_N,
    build_nba,
    dihedral_deg,
    format_xyz,
)

SCAN_ANGLES = [0, 15, 30, 45, 60, 75, 90]


def energy_en(symbols: list[str], coords: np.ndarray, basis: str, method: str) -> tuple[float, float]:
    mol = gto.Mole()
    mol.atom = [(symbols[i], tuple(coords[i])) for i in range(len(symbols))]
    mol.basis = basis
    mol.unit = "Angstrom"
    mol.verbose = 0
    mol.build()
    if method.upper() == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
    else:
        mf = scf.RHF(mol)
    mf.kernel()
    e = float(mf.e_tot)
    return e, float(mol.energy_nuc())


def run_rigid(basis: str, method: str, out: Path, angles: list[float]) -> dict:
    tag = f"rigid_{method}_{basis.replace('*', 's').replace('+', 'p')}"
    raw = ensure_dir(out / "raw" / tag)
    tables = ensure_dir(out / "tables")

    symbols, _ = build_nba(0.0)
    rows = []
    for ang in angles:
        print(f"[P3 rigid] theta_target={ang:.0f}", flush=True)
        _, coords = build_nba(ang)
        th = dihedral_deg(coords, IDX_CORTHO_N, IDX_CIPSO_N, IDX_N, IDX_C_IMINE)
        e, en = energy_en(symbols, coords, basis, method)
        ee = e - en
        rows.append(
            {
                "theta_target": ang,
                "theta_deg": th,
                "E_ha": e,
                "Ee_ha": ee,
                "EN_ha": en,
            }
        )
        (raw / f"nba_rigid_{int(ang):03d}.xyz").write_text(
            f"{len(symbols)}\nrigid theta={th:.2f} {tag}\n{format_xyz(symbols, coords)}\n",
            encoding="utf-8",
        )
        print(f"      theta={th:.2f} E={e:.8f}", flush=True)

    ref = min(rows, key=lambda r: abs(r["theta_deg"]))
    e0, ee0, en0 = ref["E_ha"], ref["Ee_ha"], ref["EN_ha"]
    table = []
    for r in rows:
        table.append(
            {
                **r,
                "dE_kcal": ha_to_kcal(r["E_ha"] - e0),
                "dEe_kcal": ha_to_kcal(r["Ee_ha"] - ee0),
                "dEN_kcal": ha_to_kcal(r["EN_ha"] - en0),
            }
        )

    tmin = table[int(np.argmin([t["dE_kcal"] for t in table]))]
    abs_th = abs(tmin["theta_deg"])
    if abs_th > 90:
        abs_th = 180 - abs_th
    agree_geo = 25.0 <= abs_th <= 65.0
    agree_decomp = (tmin["dEN_kcal"] > 0.0) and (tmin["dEe_kcal"] < 0.0)
    planar = min(table, key=lambda t: abs(t["theta_deg"]))
    agree_not_planar = tmin["dE_kcal"] < planar["dE_kcal"] - 0.05
    agree = bool(agree_geo and agree_decomp and agree_not_planar)

    pack = {
        "proposition": "P3",
        "tag": tag,
        "method": method,
        "basis": basis,
        "protocol": "rigid twist scan (no relaxation); single-point E/Ee/EN",
        "note": "Book uses relaxed PES; full relaxation pending compute budget",
        "theta_definition": "C_ortho(N-Ph)-C_ipso(N-Ph)-N-C_imine",
        "book_theta_exp_deg": [36, 55],
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "reference_planar_theta_deg": ref["theta_deg"],
        "scan": table,
        "E_min_scan": {
            "theta_deg": tmin["theta_deg"],
            "dE_kcal": tmin["dE_kcal"],
            "dEe_kcal": tmin["dEe_kcal"],
            "dEN_kcal": tmin["dEN_kcal"],
        },
        "checks": {
            "agree_geo_30_60": agree_geo,
            "agree_decomp_EN_up_Ee_down": agree_decomp,
            "agree_not_planar_global_min": agree_not_planar,
        },
        "agree": agree,
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(tables / f"pes_{tag}.json", pack)

    with (tables / f"pes_{tag}.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["theta_target", "theta_deg", "dE_kcal", "dEe_kcal", "dEN_kcal"],
        )
        w.writeheader()
        for t in table:
            w.writerow(
                {
                    "theta_target": f"{t['theta_target']:.1f}",
                    "theta_deg": f"{t['theta_deg']:.3f}",
                    "dE_kcal": f"{t['dE_kcal']:.4f}",
                    "dEe_kcal": f"{t['dEe_kcal']:.4f}",
                    "dEN_kcal": f"{t['dEN_kcal']:.4f}",
                }
            )

    lines = [
        f"# P3 NBA rigid twist scan ({tag})",
        "",
        f"- Protocol: **rigid** theta scan + **{method}/{basis}** single-point",
        f"- E_min at theta = **{tmin['theta_deg']:.1f} deg** "
        f"(dE={tmin['dE_kcal']:+.3f}, dEe={tmin['dEe_kcal']:+.3f}, "
        f"dEN={tmin['dEN_kcal']:+.3f} kcal/mol vs near-planar)",
        "- Book theta_exp ~ 36-55 deg (relaxed PES in book)",
        f"- Checks: geo={agree_geo}, decomp={agree_decomp}, not_planar={agree_not_planar}",
        f"- Auto agree flag: **{agree}** (formal VERDICT only in deliverables/)",
        "",
    ]
    (out / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    write_json(
        out / "meta.json",
        {"proposition": "P3", "updated_utc": pack["timestamp_utc"], "latest_tag": tag},
    )
    print("\n".join(lines), flush=True)
    return pack


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP", choices=["B3LYP", "RHF"])
    ap.add_argument("--out", default=str(ROOT / "results" / "P3"))
    ap.add_argument("--angles", default=",".join(str(a) for a in SCAN_ANGLES))
    args = ap.parse_args()
    angles = [float(x) for x in args.angles.split(",") if x.strip()]
    run_rigid(args.basis, args.method, Path(args.out), angles)


if __name__ == "__main__":
    main()
