"""
P3 — NBA (Ph–N=CH–Ph) relaxed twist PES + Ee/EN.

Default protocol (speed + DFT energies):
  RHF/6-31G* constrained geometry scan
  B3LYP/6-31G* single-point E, Ee, EN at each point

θ = C_ortho(N–Ph)–C_ipso(N–Ph)–N–C_imine
"""
from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

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
    dihedral_grad_deg,
    format_xyz,
    set_twist_deg,
)

BOHR = 0.52917721092


def make_mol(symbols: list[str], coords_ang: np.ndarray, basis: str) -> gto.Mole:
    mol = gto.Mole()
    mol.atom = [(symbols[i], tuple(coords_ang[i])) for i in range(len(symbols))]
    mol.basis = basis
    mol.unit = "Angstrom"
    mol.verbose = 0
    mol.build()
    return mol


def make_mf(mol: gto.Mole, method: str):
    if method.upper() == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
        return mf
    return scf.RHF(mol)


def energy_grad_en(
    symbols: list[str], coords_ang: np.ndarray, basis: str, method: str
) -> tuple[float, np.ndarray, float]:
    mol = make_mol(symbols, coords_ang, basis)
    mf = make_mf(mol, method)
    mf.kernel()
    e = float(mf.e_tot)
    g = np.asarray(mf.nuc_grad_method().kernel(), dtype=float)
    en = float(mol.energy_nuc())
    return e, g, en


def energy_en(
    symbols: list[str], coords_ang: np.ndarray, basis: str, method: str
) -> tuple[float, float]:
    mol = make_mol(symbols, coords_ang, basis)
    mf = make_mf(mol, method)
    mf.kernel()
    e = float(mf.e_tot)
    en = float(mol.energy_nuc())
    return e, en


def optimize_free(
    symbols: list[str], coords0: np.ndarray, basis: str, method: str, maxiter: int = 40
) -> dict:
    x0 = (coords0 / BOHR).ravel().copy()
    cache: dict[bytes, tuple[float, np.ndarray]] = {}

    def unpack(x: np.ndarray) -> np.ndarray:
        return x.reshape(-1, 3) * BOHR

    def fun(x: np.ndarray) -> float:
        key = x.tobytes()
        if key not in cache:
            e, g, _ = energy_grad_en(symbols, unpack(x), basis, method)
            cache[key] = (e, g.ravel())
        return cache[key][0]

    def jac(x: np.ndarray) -> np.ndarray:
        key = x.tobytes()
        if key not in cache:
            fun(x)
        return cache[key][1]

    def cb(xk):
        print(f"    free-opt iter cache={len(cache)}", flush=True)

    res = minimize(
        fun,
        x0,
        method="BFGS",
        jac=jac,
        callback=cb,
        options={"maxiter": maxiter, "gtol": 3e-4, "disp": False},
    )
    coords = unpack(res.x)
    e, _, en = energy_grad_en(symbols, coords, basis, method)
    return {
        "energy_ha": e,
        "en_ha": en,
        "ee_ha": e - en,
        "theta_deg": dihedral_deg(coords, IDX_CORTHO_N, IDX_CIPSO_N, IDX_N, IDX_C_IMINE),
        "coords": coords,
        "success": bool(res.success),
        "nit": int(res.nit),
    }


def optimize_constrained(
    symbols: list[str],
    coords0: np.ndarray,
    theta_target: float,
    basis: str,
    method: str,
    maxiter: int = 30,
    k_penalty: float = 0.10,
) -> dict:
    coords0 = set_twist_deg(coords0.copy(), theta_target)
    x0 = (coords0 / BOHR).ravel().copy()
    cache: dict[bytes, tuple[float, np.ndarray]] = {}

    def unpack(x: np.ndarray) -> np.ndarray:
        return x.reshape(-1, 3) * BOHR

    def eval_point(x: np.ndarray):
        key = x.tobytes()
        if key not in cache:
            coords = unpack(x)
            e, g_ha_bohr, _en = energy_grad_en(symbols, coords, basis, method)
            th = dihedral_deg(coords, IDX_CORTHO_N, IDX_CIPSO_N, IDX_N, IDX_C_IMINE)
            dth = th - theta_target
            while dth > 180:
                dth -= 360
            while dth < -180:
                dth += 360
            dth_rad = np.radians(dth)
            pen = 0.5 * k_penalty * dth_rad**2
            gθ = dihedral_grad_deg(coords, IDX_CORTHO_N, IDX_CIPSO_N, IDX_N, IDX_C_IMINE)
            g_pen = k_penalty * dth_rad * (np.pi / 180.0) * gθ * BOHR
            cache[key] = (e + pen, (g_ha_bohr + g_pen).ravel())
        return cache[key]

    def cb(xk):
        print(f"    constr θ={theta_target:.0f} cache={len(cache)}", flush=True)

    res = minimize(
        lambda x: eval_point(x)[0],
        x0,
        method="BFGS",
        jac=lambda x: eval_point(x)[1],
        callback=cb,
        options={"maxiter": maxiter, "gtol": 4e-4, "disp": False},
    )
    coords = unpack(res.x)
    e, _, en = energy_grad_en(symbols, coords, basis, method)
    th = dihedral_deg(coords, IDX_CORTHO_N, IDX_CIPSO_N, IDX_N, IDX_C_IMINE)
    err = abs(((th - theta_target + 180) % 360) - 180)
    return {
        "theta_target": theta_target,
        "theta_deg": th,
        "energy_ha": e,
        "en_ha": en,
        "ee_ha": e - en,
        "coords": coords,
        "success": bool(res.success),
        "nit": int(res.nit),
        "theta_error_deg": err,
    }


def annotate_energy(pt: dict, symbols: list[str], basis: str, energy_method: str, opt_method: str) -> dict:
    out = dict(pt)
    if energy_method != opt_method:
        e, en = energy_en(symbols, pt["coords"], basis, energy_method)
        out["energy_ha"] = e
        out["en_ha"] = en
        out["ee_ha"] = e - en
    return out


def run_scan(
    basis: str,
    opt_method: str,
    energy_method: str,
    out: Path,
    angles: list[float],
    maxiter: int,
    skip_free: bool = True,
) -> dict:
    if energy_method == opt_method:
        tag = f"{opt_method}_{basis.replace('*', 's').replace('+', 'p')}"
    else:
        tag = f"{energy_method}_on_{opt_method}_{basis.replace('*', 's').replace('+', 'p')}"
    raw = ensure_dir(out / "raw" / tag)
    tables = ensure_dir(out / "tables")

    symbols, coords = build_nba(40.0)
    free = None
    if skip_free:
        print(f"[P3] skip free opt; seed twist=40 deg; tag={tag}", flush=True)
        e0, en0 = energy_en(symbols, coords, basis, energy_method)
        free = {
            "theta_deg": dihedral_deg(coords, IDX_CORTHO_N, IDX_CIPSO_N, IDX_N, IDX_C_IMINE),
            "energy_ha": e0,
            "en_ha": en0,
            "ee_ha": e0 - en0,
            "coords": coords,
            "success": True,
            "nit": 0,
        }
    else:
        print(f"[P3] free opt ({opt_method}) tag={tag}", flush=True)
        free = annotate_energy(
            optimize_free(symbols, coords, basis, opt_method, maxiter=maxiter),
            symbols,
            basis,
            energy_method,
            opt_method,
        )
    (raw / "nba_free.xyz").write_text(
        f"{len(symbols)}\nNBA free {tag} theta={free['theta_deg']:.2f}\n"
        + format_xyz(symbols, free["coords"])
        + "\n",
        encoding="utf-8",
    )
    print(
        f"[P3] seed/free theta={free['theta_deg']:.2f} E={free['energy_ha']:.8f} nit={free['nit']}",
        flush=True,
    )

    rows = []
    coords_prev = free["coords"]
    for ang in angles:
        print(f"[P3] constrained -> {ang:.0f}", flush=True)
        pt = annotate_energy(
            optimize_constrained(symbols, coords_prev, ang, basis, opt_method, maxiter=maxiter),
            symbols,
            basis,
            energy_method,
            opt_method,
        )
        coords_prev = pt["coords"]
        (raw / f"nba_theta_{int(ang):03d}.xyz").write_text(
            f"{len(symbols)}\nNBA target={ang} theta={pt['theta_deg']:.2f} {tag}\n"
            + format_xyz(symbols, pt["coords"])
            + "\n",
            encoding="utf-8",
        )
        rows.append(pt)
        print(
            f"      theta={pt['theta_deg']:.2f} err={pt['theta_error_deg']:.2f} "
            f"E={pt['energy_ha']:.8f} nit={pt['nit']}",
            flush=True,
        )

    ref = min(rows, key=lambda r: abs(r["theta_deg"]))
    e0, ee0, en0 = ref["energy_ha"], ref["ee_ha"], ref["en_ha"]
    table = []
    for r in rows:
        table.append(
            {
                "theta_target": r["theta_target"],
                "theta_deg": r["theta_deg"],
                "theta_error_deg": r["theta_error_deg"],
                "E_ha": r["energy_ha"],
                "Ee_ha": r["ee_ha"],
                "EN_ha": r["en_ha"],
                "dE_kcal": ha_to_kcal(r["energy_ha"] - e0),
                "dEe_kcal": ha_to_kcal(r["ee_ha"] - ee0),
                "dEN_kcal": ha_to_kcal(r["en_ha"] - en0),
                "nit": r["nit"],
                "success": r["success"],
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
        "opt_method": opt_method,
        "energy_method": energy_method,
        "basis": basis,
        "protocol": (
            f"{opt_method} constrained relaxed scan (maxiter={maxiter}); "
            f"{energy_method} single-point energies; skip_free={skip_free}"
        ),
        "theta_definition": "C_ortho(N-Ph)-C_ipso(N-Ph)-N-C_imine",
        "book_theta_exp_deg": [36, 55],
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "free_opt": {
            "theta_deg": free["theta_deg"],
            "E_ha": free["energy_ha"],
            "Ee_ha": free["ee_ha"],
            "EN_ha": free["en_ha"],
            "success": free["success"],
            "nit": free["nit"],
        },
        "reference_planar_scan_theta_deg": ref["theta_deg"],
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
            fieldnames=[
                "theta_target",
                "theta_deg",
                "dE_kcal",
                "dEe_kcal",
                "dEN_kcal",
                "theta_error_deg",
            ],
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
                    "theta_error_deg": f"{t['theta_error_deg']:.3f}",
                }
            )

    lines = [
        f"# P3 NBA twist PES ({tag})",
        "",
        f"- Protocol: **{opt_method}** geometry + **{energy_method}** energy",
        f"- Free opt theta = **{free['theta_deg']:.1f} deg**",
        f"- Scan E_min at theta = **{tmin['theta_deg']:.1f} deg** "
        f"(dE={tmin['dE_kcal']:+.3f}, dEe={tmin['dEe_kcal']:+.3f}, "
        f"dEN={tmin['dEN_kcal']:+.3f} kcal/mol vs near-planar)",
        "- Book theta_exp ~ 36-55 deg",
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
    ap.add_argument("--opt-method", default="RHF", choices=["B3LYP", "RHF"])
    ap.add_argument("--energy-method", default="B3LYP", choices=["B3LYP", "RHF"])
    ap.add_argument("--out", default=str(ROOT / "results" / "P3"))
    ap.add_argument("--angles", default="15,30,45,60,90")
    ap.add_argument("--maxiter", type=int, default=12)
    ap.add_argument(
        "--skip-free",
        action="store_true",
        default=True,
        help="Skip free optimization (default on)",
    )
    ap.add_argument("--do-free", action="store_true", help="Run free opt before scan")
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out)
    angles = [float(x) for x in args.angles.split(",") if x.strip()]
    skip_free = not args.do_free
    run_scan(
        args.basis,
        args.opt_method,
        args.energy_method,
        out,
        angles,
        args.maxiter,
        skip_free=skip_free,
    )


if __name__ == "__main__":
    main()
