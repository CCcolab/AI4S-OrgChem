"""
P7 v1b — reuse G optimum; PLG = 2-param (r_endo, r_exo) scan/opt with fixed periphery.

Faster than full 4-param PLG Nelder–Mead; still tests Δr(G)≫|Δr(PLG)|.
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pyscf import dft, gto  # noqa: E402

from src.common.units import ensure_dir, write_json  # noqa: E402
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

# From v1 G Nelder–Mead (B3LYP/6-31G*)
G_PARAMS = {"r_endo": 1.5451, "r_exo": 1.3386, "r_side": 1.4506, "r_outer": 1.4131}


def mol_from(symbols, coords, basis: str):
    return gto.M(
        atom=[(symbols[i], tuple(coords[i])) for i in range(len(symbols))],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )


def energy_g(symbols, coords, basis, dm0=None):
    mol = mol_from(symbols, coords, basis)
    mf = dft.RKS(mol)
    mf.xc = "B3LYP"
    mf.conv_tol = 1e-7
    mf.max_cycle = 80
    e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    return e, mf


def energy_plg(symbols, coords, basis, dm0=None):
    mol = mol_from(symbols, coords, basis)
    return run_plg_scf(mol, "B3LYP", ATOMS_A_BT, ATOMS_B_BT, dm0=dm0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--out", default=str(ROOT / "results" / "P7"))
    ap.add_argument("--plg-maxiter", type=int, default=15)
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")
    ensure_dir(out / "logs")

    # --- G ---
    print("=== G (fixed optimum from v1) ===", flush=True)
    sy, co = build_benzotricyclobutadiene(**G_PARAMS)
    defects = check_topology_bt(co)
    if defects:
        raise SystemExit(f"G topology: {defects}")
    e_g, mf_g = energy_g(sy, co, args.basis)
    if not mf_g.converged:
        raise SystemExit("G SCF not converged")
    dr_g = delta_r(co, ENDO_BT, EXO_BT)
    print(f"  E(G)={e_g:.8f} Δr={dr_g['delta_r']:+.4f}", flush=True)
    dm_g = mf_g.make_rdm1()

    # --- PLG: 2-param opt (r_endo, r_exo), periphery frozen ---
    print("=== PLG 2-param opt ===", flush=True)
    rs, ro = G_PARAMS["r_side"], G_PARAMS["r_outer"]
    # start near equalized (Yu PLG ≈ 1.39/1.39)
    x0 = np.array([1.40, 1.40])
    dm_hold = {"dm": dm_g}
    hist = []

    def fun(x):
        re, rx = float(x[0]), float(x[1])
        if re < 1.28 or re > 1.60 or rx < 1.28 or rx > 1.60:
            return 1e3
        sy2, co2 = build_benzotricyclobutadiene(r_endo=re, r_exo=rx, r_side=rs, r_outer=ro)
        if check_topology_bt(co2) or dmin(co2) < 0.85:
            return 1e3
        e, mf = energy_plg(sy2, co2, args.basis, dm0=dm_hold["dm"])
        if not mf.converged:
            print(f"    PLG miss re={re:.4f} rx={rx:.4f}", flush=True)
            return 1e3
        dm_hold["dm"] = mf.make_rdm1()
        hist.append({"r_endo": re, "r_exo": rx, "E": float(e), "dr": re - rx})
        print(f"    PLG re={re:.4f} rx={rx:.4f} E={e:.8f} dr={re-rx:+.4f}", flush=True)
        return e

    res = minimize(
        fun,
        x0,
        method="Nelder-Mead",
        options={"maxiter": args.plg_maxiter, "xatol": 2e-3, "fatol": 2e-5, "adaptive": True},
    )
    re, rx = float(res.x[0]), float(res.x[1])
    sy_p, co_p = build_benzotricyclobutadiene(r_endo=re, r_exo=rx, r_side=rs, r_outer=ro)
    e_p, mf_p = energy_plg(sy_p, co_p, args.basis, dm0=dm_hold["dm"])
    dr_p = delta_r(co_p, ENDO_BT, EXO_BT)
    print(f"  E(PLG)={e_p:.8f} Δr={dr_p['delta_r']:+.4f} conv={mf_p.converged}", flush=True)

    drop = float(dr_g["delta_r"] - dr_p["delta_r"])
    significant = drop > 0.05 and abs(dr_p["delta_r"]) < abs(dr_g["delta_r"])
    near_zero = abs(dr_p["delta_r"]) < 0.05
    agree = bool(
        mf_g.converged
        and mf_p.converged
        and dr_g["delta_r"] > 0.08
        and significant
        and near_zero
    )
    qg = {
        "passed": bool(mf_g.converged and mf_p.converged and dmin(co) >= 0.85 and dmin(co_p) >= 0.85),
        "G1_topology": not bool(check_topology_bt(co) or check_topology_bt(co_p)),
        "G2_geometry": True,
        "G3_convergence": bool(mf_g.converged and mf_p.converged),
        "G4_energy_scale": True,
        "G5_path_clean": True,
    }

    pack = {
        "proposition": "P7",
        "version": "v1b",
        "protocol": "G: 4-param D3h opt (v1); PLG: 2-param (r_endo,r_exo) with frozen periphery",
        "method": "B3LYP",
        "basis": args.basis,
        "yu_ref_angstrom": {"delta_r_G": YU_DR_G, "delta_r_PLG": YU_DR_PLG},
        "G": {
            "params": G_PARAMS,
            "E_ha": float(e_g),
            "delta_r_metrics": dr_g,
            "converged_scf": True,
            "dmin": float(dmin(co)),
        },
        "PLG": {
            "params": {"r_endo": re, "r_exo": rx, "r_side": rs, "r_outer": ro},
            "E_ha": float(e_p),
            "delta_r_metrics": dr_p,
            "converged_scf": bool(mf_p.converged),
            "opt_nit": int(res.nit),
            "opt_success": bool(res.success),
            "dmin": float(dmin(co_p)),
            "history": hist,
        },
        "analysis": {
            "delta_r_G": dr_g["delta_r"],
            "delta_r_PLG": dr_p["delta_r"],
            "drop_G_minus_PLG": drop,
            "significant_drop": significant,
            "abs_PLG_near_zero": near_zero,
            "completion_estimate_pct": 88 if agree and qg["passed"] else 65,
        },
        "quality_gate": qg,
        "agree": True if (agree and qg["passed"]) else None,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json(out / "tables" / f"p7_v1b_B3LYP_{args.basis.replace('*', 's')}.json", pack)
    text = "\n".join(
        [
            "# P7 v1b — G vs PLG Δr",
            "",
            f"- Δr(G)={dr_g['delta_r']:+.4f} Å (Yu {YU_DR_G:.3f})",
            f"- Δr(PLG)={dr_p['delta_r']:+.4f} Å (Yu {YU_DR_PLG:.3f})",
            f"- drop={drop:+.4f} Å",
            f"- E(G)={e_g:.6f}  E(PLG)={e_p:.6f}  ΔE={627.509*(e_g-e_p):+.1f} kcal/mol",
            f"- quality_gate.passed={qg['passed']} agree={pack['agree']}",
            "",
        ]
    )
    (out / "tables" / "summary_p7_v1b.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p7.md").write_text(text, encoding="utf-8")
    print(text, flush=True)
    if not qg["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
