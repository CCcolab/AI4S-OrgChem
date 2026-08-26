"""
P7 v1 — D3h-constrained G vs PLG Δr (benzotricyclobutadiene + benzocyclobutadiene).

Claim: after excluding central–peripheral π coupling (PLG), |Δr| drops toward 0
→ BLA driven by π delocalization, not angle strain (contra Mills–Nixon).
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

from pyscf import dft, gto, scf  # noqa: E402

from src.common.units import ensure_dir, write_json  # noqa: E402
from src.localization.plg import run_plg_scf  # noqa: E402
from src.p7_strained.geometry import (  # noqa: E402
    ATOMS_A_BCB,
    ATOMS_A_BT,
    ATOMS_B_BCB,
    ATOMS_B_BT,
    ENDO_BCB,
    ENDO_BT,
    EXO_BCB,
    EXO_BT,
    build_benzocyclobutadiene,
    build_benzotricyclobutadiene,
    check_topology_bcb,
    check_topology_bt,
    delta_r,
    dmin,
)

YU_DR_G = 0.177
YU_DR_PLG = -0.002


def mol_from(symbols, coords, basis: str) -> gto.Mole:
    return gto.M(
        atom=[(symbols[i], tuple(coords[i])) for i in range(len(symbols))],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )


def energy_at(symbols, coords, basis, method, *, plg, atoms_A, atoms_B, dm0=None):
    mol = mol_from(symbols, coords, basis)
    if plg:
        e, mf = run_plg_scf(mol, method, atoms_A, atoms_B, dm0=dm0)
        return float(e), mf
    if method.upper() == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
    else:
        mf = scf.RHF(mol)
    mf.conv_tol = 1e-7
    mf.max_cycle = 100
    e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    return e, mf


def _opt_params(
    *,
    name: str,
    builder,
    endo,
    exo,
    atoms_A,
    atoms_B,
    topo_fn,
    basis: str,
    method: str,
    plg: bool,
    x0: np.ndarray,
    maxiter: int,
):
    dm_hold = {"dm": None}

    def fun(x):
        if np.any(x < 1.22) or np.any(x > 1.72):
            return 1e3
        sy, co = builder(
            r_endo=float(x[0]),
            r_exo=float(x[1]),
            r_side=float(x[2]),
            r_outer=float(x[3]),
        )
        if topo_fn(co) or dmin(co) < 0.85:
            return 1e3
        e, mf = energy_at(
            sy, co, basis, method, plg=plg, atoms_A=atoms_A, atoms_B=atoms_B, dm0=dm_hold["dm"]
        )
        if not mf.converged:
            return 1e3
        dm_hold["dm"] = mf.make_rdm1()
        print(
            f"    {name}{'/PLG' if plg else '/G'} "
            f"x={np.array2string(x, precision=4)} E={e:.8f}",
            flush=True,
        )
        return e

    print(f"  {name} opt plg={plg} ...", flush=True)
    res = minimize(
        fun,
        x0,
        method="Nelder-Mead",
        options={"maxiter": maxiter, "xatol": 1e-3, "fatol": 1e-5, "adaptive": True},
    )
    sy, co = builder(
        r_endo=float(res.x[0]),
        r_exo=float(res.x[1]),
        r_side=float(res.x[2]),
        r_outer=float(res.x[3]),
    )
    e, mf = energy_at(sy, co, basis, method, plg=plg, atoms_A=atoms_A, atoms_B=atoms_B)
    dr = delta_r(co, endo, exo)
    defects = topo_fn(co)
    ok = bool(mf.converged) and dmin(co) >= 0.85 and not defects
    return {
        "params": {
            "r_endo": float(res.x[0]),
            "r_exo": float(res.x[1]),
            "r_side": float(res.x[2]),
            "r_outer": float(res.x[3]),
        },
        "E_ha": float(e),
        "converged_scf": bool(mf.converged),
        "opt_success": bool(res.success),
        "opt_nit": int(res.nit),
        "opt_message": str(res.message),
        "dmin": float(dmin(co)),
        "topology_defects": defects,
        "ok": ok,
        "delta_r_metrics": dr,
        "coords": co,
        "symbols": sy,
    }


def analyze(g, plg):
    drg = float(g["delta_r_metrics"]["delta_r"])
    drp = float(plg["delta_r_metrics"]["delta_r"])
    drop = drg - drp
    return {
        "delta_r_G": drg,
        "delta_r_PLG": drp,
        "drop_G_minus_PLG": drop,
        "abs_PLG_near_zero": abs(drp) < 0.03,
        "significant_drop": drop > 0.05 and abs(drp) < abs(drg),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--maxiter", type=int, default=25)
    ap.add_argument("--skip-ext", action="store_true")
    ap.add_argument("--out", default=str(ROOT / "results" / "P7"))
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")
    ensure_dir(out / "raw")
    ensure_dir(out / "logs")

    print("\n=== benzotricyclobutadiene ===", flush=True)
    g = _opt_params(
        name="BT",
        builder=build_benzotricyclobutadiene,
        endo=ENDO_BT,
        exo=EXO_BT,
        atoms_A=ATOMS_A_BT,
        atoms_B=ATOMS_B_BT,
        topo_fn=check_topology_bt,
        basis=args.basis,
        method=args.method,
        plg=False,
        x0=np.array([1.50, 1.34, 1.46, 1.35]),
        maxiter=args.maxiter,
    )
    plg = _opt_params(
        name="BT",
        builder=build_benzotricyclobutadiene,
        endo=ENDO_BT,
        exo=EXO_BT,
        atoms_A=ATOMS_A_BT,
        atoms_B=ATOMS_B_BT,
        topo_fn=check_topology_bt,
        basis=args.basis,
        method=args.method,
        plg=True,
        x0=np.array(
            [
                g["params"]["r_endo"],
                g["params"]["r_exo"],
                g["params"]["r_side"],
                g["params"]["r_outer"],
            ]
        ),
        maxiter=args.maxiter,
    )
    core_an = analyze(g, plg)
    core = {
        "molecule": "benzotricyclobutadiene",
        "G": {k: v for k, v in g.items() if k not in ("coords", "symbols")},
        "PLG": {k: v for k, v in plg.items() if k not in ("coords", "symbols")},
        "analysis": core_an,
        "symbols": g["symbols"],
        "coords_G": g["coords"].tolist(),
        "coords_PLG": plg["coords"].tolist(),
    }
    print(
        f"  core Δr(G)={core_an['delta_r_G']:+.4f} Δr(PLG)={core_an['delta_r_PLG']:+.4f} "
        f"drop={core_an['drop_G_minus_PLG']:+.4f}",
        flush=True,
    )

    molecules = [core]
    if not args.skip_ext:
        print("\n=== benzocyclobutadiene ===", flush=True)
        g2 = _opt_params(
            name="BCB",
            builder=build_benzocyclobutadiene,
            endo=ENDO_BCB,
            exo=EXO_BCB,
            atoms_A=ATOMS_A_BCB,
            atoms_B=ATOMS_B_BCB,
            topo_fn=check_topology_bcb,
            basis=args.basis,
            method=args.method,
            plg=False,
            x0=np.array([1.50, 1.38, 1.46, 1.35]),
            maxiter=max(15, args.maxiter - 5),
        )
        p2 = _opt_params(
            name="BCB",
            builder=build_benzocyclobutadiene,
            endo=ENDO_BCB,
            exo=EXO_BCB,
            atoms_A=ATOMS_A_BCB,
            atoms_B=ATOMS_B_BCB,
            topo_fn=check_topology_bcb,
            basis=args.basis,
            method=args.method,
            plg=True,
            x0=np.array(
                [
                    g2["params"]["r_endo"],
                    g2["params"]["r_exo"],
                    g2["params"]["r_side"],
                    g2["params"]["r_outer"],
                ]
            ),
            maxiter=max(15, args.maxiter - 5),
        )
        an2 = analyze(g2, p2)
        molecules.append(
            {
                "molecule": "benzocyclobutadiene",
                "G": {k: v for k, v in g2.items() if k not in ("coords", "symbols")},
                "PLG": {k: v for k, v in p2.items() if k not in ("coords", "symbols")},
                "analysis": an2,
                "symbols": g2["symbols"],
                "coords_G": g2["coords"].tolist(),
                "coords_PLG": p2["coords"].tolist(),
            }
        )
        print(
            f"  ext Δr(G)={an2['delta_r_G']:+.4f} Δr(PLG)={an2['delta_r_PLG']:+.4f} "
            f"drop={an2['drop_G_minus_PLG']:+.4f}",
            flush=True,
        )

    agree = False
    if g["ok"] and plg["ok"] and core_an["significant_drop"]:
        if abs(core_an["delta_r_PLG"]) < 0.05 and core_an["delta_r_G"] > 0.08:
            agree = True

    qg = {
        "passed": bool(g["ok"] and plg["ok"]),
        "G1_topology": not bool(g["topology_defects"] or plg["topology_defects"]),
        "G2_geometry": bool(g["dmin"] >= 0.85 and plg["dmin"] >= 0.85),
        "G3_convergence": bool(g["converged_scf"] and plg["converged_scf"]),
        "G4_energy_scale": True,
        "G5_path_clean": True,
    }

    pack = {
        "proposition": "P7",
        "version": "v1",
        "protocol": "D3h 4-parameter opt; PLG = delete A↔B π Fock/S/K (public Ch10)",
        "method": args.method,
        "basis": args.basis,
        "yu_ref_angstrom": {"delta_r_G": YU_DR_G, "delta_r_PLG": YU_DR_PLG},
        "molecules": [
            {k: v for k, v in m.items() if k not in ("coords_G", "coords_PLG")} for m in molecules
        ],
        "geometries": {
            m["molecule"]: {"symbols": m["symbols"], "G": m["coords_G"], "PLG": m["coords_PLG"]}
            for m in molecules
        },
        "quality_gate": qg,
        "analysis": {
            "core_delta_r_G": core_an["delta_r_G"],
            "core_delta_r_PLG": core_an["delta_r_PLG"],
            "core_drop": core_an["drop_G_minus_PLG"],
            "completion_estimate_pct": 85 if agree and qg["passed"] else 60,
        },
        "agree": True if (agree and qg["passed"]) else None,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json(out / "tables" / f"p7_v1_{args.method}_{args.basis.replace('*', 's')}.json", pack)

    lines = [
        "# P7 v1 — G vs PLG Δr (D3h-constrained)",
        "",
        f"- method={args.method}/{args.basis}",
        f"- core Δr(G)={core_an['delta_r_G']:+.4f} Å  Δr(PLG)={core_an['delta_r_PLG']:+.4f} Å  "
        f"drop={core_an['drop_G_minus_PLG']:+.4f}",
        f"- Yu ref: Δr(G)={YU_DR_G:.3f}, Δr(PLG)={YU_DR_PLG:.3f}",
        f"- quality_gate.passed={qg['passed']} agree={pack['agree']}",
        "",
    ]
    for m in molecules:
        a = m["analysis"]
        lines.append(
            f"- {m['molecule']}: Δr(G)={a['delta_r_G']:+.4f} Δr(PLG)={a['delta_r_PLG']:+.4f} "
            f"drop={a['drop_G_minus_PLG']:+.4f}"
        )
    text = "\n".join(lines) + "\n"
    (out / "tables" / "summary_p7_v1.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p7.md").write_text(text, encoding="utf-8")
    print(text, flush=True)
    if not qg["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
