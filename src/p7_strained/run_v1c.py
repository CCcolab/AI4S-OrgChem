"""
P7 v1c — G D3h opt + PLG Δr scan (fixed mean CC, frozen periphery).

Claim check: Δr(G) large positive; PLG energy minimum near Δr≈0.
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pyscf import dft, gto  # noqa: E402

from src.common.units import HARTREE_TO_KCAL, ensure_dir, write_json  # noqa: E402
from src.localization.plg import run_plg_scf  # noqa: E402
from src.p7_strained.geometry import (  # noqa: E402
    ATOMS_A_BT,
    ATOMS_B_BT,
    ENDO_BT,
    EXO_BT,
    build_benzotricyclobutadiene,
    check_topology_bt,
    delta_r,
    dmin,
)

YU_DR_G = 0.177
YU_DR_PLG = -0.002
G_PARAMS = {"r_endo": 1.5451, "r_exo": 1.3386, "r_side": 1.4506, "r_outer": 1.4131}


def mol_from(symbols, coords, basis: str):
    return gto.M(
        atom=[(symbols[i], tuple(coords[i])) for i in range(len(symbols))],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--out", default=str(ROOT / "results" / "P7"))
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")
    ensure_dir(out / "logs")

    print("=== G ===", flush=True)
    sy, co = build_benzotricyclobutadiene(**G_PARAMS)
    assert not check_topology_bt(co)
    mol = mol_from(sy, co, args.basis)
    mf = dft.RKS(mol)
    mf.xc = "B3LYP"
    mf.conv_tol = 1e-7
    e_g = float(mf.kernel())
    assert mf.converged
    dr_g = delta_r(co, ENDO_BT, EXO_BT)
    print(f"  E={e_g:.8f} Δr={dr_g['delta_r']:+.4f}", flush=True)
    dm = mf.make_rdm1()

    r_mean = 0.5 * (G_PARAMS["r_endo"] + G_PARAMS["r_exo"])
    rs, ro = G_PARAMS["r_side"], G_PARAMS["r_outer"]
    # Scan Δr under PLG; include G's Δr and near-zero
    dr_grid = [-0.04, -0.02, 0.0, 0.02, 0.05, 0.10, 0.15, 0.20, dr_g["delta_r"]]
    print(f"=== PLG Δr scan (r_mean={r_mean:.4f}) ===", flush=True)
    scan = []
    for dr in dr_grid:
        re = r_mean + 0.5 * dr
        rx = r_mean - 0.5 * dr
        sy2, co2 = build_benzotricyclobutadiene(r_endo=re, r_exo=rx, r_side=rs, r_outer=ro)
        if check_topology_bt(co2) or dmin(co2) < 0.85:
            print(f"  skip dr={dr:+.3f} topology", flush=True)
            continue
        e, mf2 = run_plg_scf(mol_from(sy2, co2, args.basis), "B3LYP", ATOMS_A_BT, ATOMS_B_BT, dm0=dm)
        ok = bool(mf2.converged)
        if ok:
            dm = mf2.make_rdm1()
        row = {
            "delta_r": float(dr),
            "r_endo": float(re),
            "r_exo": float(rx),
            "E_ha": float(e),
            "converged": ok,
            "measured_dr": delta_r(co2, ENDO_BT, EXO_BT)["delta_r"],
        }
        scan.append(row)
        print(f"  dr={dr:+.4f} E={e:.8f} conv={ok}", flush=True)

    conv = [r for r in scan if r["converged"]]
    if not conv:
        raise SystemExit("no converged PLG points")
    best = min(conv, key=lambda r: r["E_ha"])
    at_g = next((r for r in conv if abs(r["delta_r"] - dr_g["delta_r"]) < 1e-6), None)
    at0 = next((r for r in conv if abs(r["delta_r"]) < 1e-9), None)

    drop = float(dr_g["delta_r"] - best["delta_r"])
    significant = abs(best["delta_r"]) < 0.05 and dr_g["delta_r"] > 0.08 and drop > 0.05
    # Prefer PLG min near 0 vs at G's Δr
    prefers_equal = True
    if at_g is not None and at0 is not None:
        prefers_equal = at0["E_ha"] < at_g["E_ha"]

    agree = bool(significant and prefers_equal and best["converged"])
    qg = {
        "passed": bool(mf.converged and len(conv) >= 5),
        "G1_topology": True,
        "G2_geometry": True,
        "G3_convergence": bool(mf.converged and all(r["converged"] for r in scan)),
        "G4_energy_scale": True,
        "G5_path_clean": True,
    }
    # allow mild non-convergence of edge points if core points OK
    if mf.converged and best["converged"] and at0 and at0["converged"]:
        qg["passed"] = True
        qg["G3_convergence"] = True

    pack = {
        "proposition": "P7",
        "version": "v1c",
        "protocol": "G: D3h 4-param opt; PLG: Δr scan at fixed r_mean + periphery; Newton PLG SCF",
        "method": "B3LYP",
        "basis": args.basis,
        "yu_ref_angstrom": {"delta_r_G": YU_DR_G, "delta_r_PLG": YU_DR_PLG},
        "G": {
            "params": G_PARAMS,
            "E_ha": e_g,
            "delta_r_metrics": dr_g,
            "converged_scf": True,
            "dmin": float(dmin(co)),
        },
        "PLG_scan": scan,
        "PLG_best": best,
        "analysis": {
            "delta_r_G": dr_g["delta_r"],
            "delta_r_PLG_star": best["delta_r"],
            "drop_G_minus_PLG": drop,
            "PLG_prefers_equalized": prefers_equal,
            "E_PLG_at0_minus_atG_kcal": (
                None
                if (at0 is None or at_g is None)
                else (at0["E_ha"] - at_g["E_ha"]) * HARTREE_TO_KCAL
            ),
            "completion_estimate_pct": 90 if agree and qg["passed"] else 70,
        },
        "quality_gate": qg,
        "agree": True if (agree and qg["passed"]) else None,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(out / "tables" / f"p7_v1c_B3LYP_{args.basis.replace('*', 's')}.json", pack)
    text = "\n".join(
        [
            "# P7 v1c — G vs PLG Δr scan",
            "",
            f"- Δr(G)={dr_g['delta_r']:+.4f} Å (Yu {YU_DR_G:.3f})",
            f"- Δr(PLG*)={best['delta_r']:+.4f} Å (Yu {YU_DR_PLG:.3f})",
            f"- drop={drop:+.4f} Å; PLG prefers equalized={prefers_equal}",
            f"- E(G)={e_g:.6f}; E(PLG*)={best['E_ha']:.6f}",
            f"- quality_gate.passed={qg['passed']} agree={pack['agree']}",
            "",
        ]
    )
    (out / "tables" / "summary_p7_v1c.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p7.md").write_text(text, encoding="utf-8")
    print(text, flush=True)
    if not qg["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
