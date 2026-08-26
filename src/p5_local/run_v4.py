"""
P5 deepen v4 — Ci hexatriene (fix GE2); adiabatic 1D; benzene Kekulé GE-m.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pyscf import dft, gto  # noqa: E402

from src.common.units import HARTREE_TO_KCAL, ensure_dir, ha_to_kcal, write_json  # noqa: E402
from src.localization.gl_2007 import make_localized_mf  # noqa: E402
from src.localization.molecules import (  # noqa: E402
    build_benzene_kekule,
    build_hexatriene,
    check_polyene_topology,
    set_single_bond_keep_sides,
    to_xyz,
)


def mol_from(symbols, coords, basis: str) -> gto.Mole:
    return gto.M(
        atom=[(symbols[i], tuple(coords[i])) for i in range(len(symbols))],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )


def scf_energy(symbols, coords, basis, method, doubles, *, mode: str, allow_pair=None, dm0=None):
    mol = mol_from(symbols, coords, basis)
    if mode == "G":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
    else:
        mf = make_localized_mf(
            mol, method, doubles, allow_pair=allow_pair, zero_overlap=True
        ).newton()
    e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"SCF failed {mode} {allow_pair}")
    return e, mf


def bond(coords, i, j) -> float:
    return float(np.linalg.norm(coords[i] - coords[j]))


def scan_min(symbols, coords0, basis, method, doubles, n_c, *, ij, right, r_grid, mode, allow_pair=None, dm_ref=None):
    best = None
    dm = dm_ref
    for r in r_grid:
        c = set_single_bond_keep_sides(coords0, ij[0], ij[1], float(r), right)
        if check_polyene_topology(symbols, c, n_c, doubles):
            continue
        e, mf = scf_energy(symbols, c, basis, method, doubles, mode=mode, allow_pair=allow_pair, dm0=dm)
        dm = mf.make_rdm1()
        print(f"    {mode}{allow_pair or ''} r={r:.3f} E={e:.8f}", flush=True)
        if best is None or e < best["E_ha"]:
            best = {"r": float(r), "E_ha": e, "coords": c, "dm": dm}
    if best is None:
        raise RuntimeError(f"empty scan {mode} {allow_pair}")
    best["at_edge"] = abs(best["r"] - float(r_grid[0])) < 1e-12 or abs(best["r"] - float(r_grid[-1])) < 1e-12
    return best


def butadiene_from_p1() -> dict:
    print("\n=== Butadiene (P1) ===", flush=True)
    gl = json.loads((ROOT / "results/P1/tables/gl2007_butadiene_B3LYP_6-31gs.json").read_text(encoding="utf-8"))
    primary = next(p for p in gl["protocols"] if p.get("zero_overlap"))
    e_g = primary["G"]["E_ha"]
    r_g = primary["G"]["r23_ang"]
    e_gl = primary["GL"]["E_ha"]
    r_gl = primary["GL"]["r23_ang"]
    dE = float(primary.get("deltaE_kcal", ha_to_kcal(e_g - e_gl)))
    dr = float(primary.get("delta_r23_ang", r_g - r_gl))
    print(f"  ΔEAm={dE:+.3f}  Δr={dr:+.4f}", flush=True)
    return {
        "molecule": "butadiene",
        "source": "P1",
        "deltaEAm_kcal": dE,
        "delta_r_ang": dr,
        "deltaEAm_positive": dE > 0,
        "delta_r_positive": dr > 0,
        "G": {"E_ha": e_g, "r23": r_g},
        "GL": {"E_ha": e_gl, "r23": r_gl},
    }


def hexatriene_adiabatic(basis: str, method: str, out: Path) -> dict:
    print("\n=== Hexatriene (Ci + adiabatic 1D) ===", flush=True)
    symbols, coords, doubles, _ = build_hexatriene(1.340, 1.450)
    n_c = 6
    ci_err = float(np.max(np.abs(coords[:6] + coords[5::-1])))
    print(f"  Ci carbon residual max={ci_err:.2e}", flush=True)
    defects = check_polyene_topology(symbols, coords, n_c, doubles)
    if defects:
        raise SystemExit(defects)

    rights = [[2, 3, 4, 5, 9, 10, 11, 12, 13], [4, 5, 11, 12, 13]]
    grid = np.linspace(1.42, 1.50, 9)

    e_g, _ = scf_energy(symbols, coords, basis, method, doubles, mode="G")
    e_gl_v, mf_gl = scf_energy(symbols, coords, basis, method, doubles, mode="GL")
    dm0 = mf_gl.make_rdm1()
    e_ge1_v, _ = scf_energy(symbols, coords, basis, method, doubles, mode="GE", allow_pair=(0, 1), dm0=dm0)
    e_ge2_v, _ = scf_energy(symbols, coords, basis, method, doubles, mode="GE", allow_pair=(1, 2), dm0=dm0)
    print(
        f"  vertical: ΔEA={ha_to_kcal(e_g - e_gl_v):+.3f} "
        f"ΔEA1={ha_to_kcal(e_ge1_v - e_gl_v):+.3f} "
        f"ΔEA2={ha_to_kcal(e_ge2_v - e_gl_v):+.3f}",
        flush=True,
    )

    print("  GL scan r12 then r34", flush=True)
    gl12 = scan_min(symbols, coords, basis, method, doubles, n_c, ij=(1, 2), right=rights[0], r_grid=grid, mode="GL", dm_ref=dm0)
    gl34 = scan_min(symbols, gl12["coords"], basis, method, doubles, n_c, ij=(3, 4), right=rights[1], r_grid=grid, mode="GL", dm_ref=gl12["dm"])
    coords_gl = gl34["coords"]
    e_gl, mf_gl = scf_energy(symbols, coords_gl, basis, method, doubles, mode="GL")
    dm_gl = mf_gl.make_rdm1()
    r12_gl, r34_gl = bond(coords_gl, 1, 2), bond(coords_gl, 3, 4)

    print("  GE1 / GE2 scans", flush=True)
    ge1 = scan_min(symbols, coords_gl, basis, method, doubles, n_c, ij=(1, 2), right=rights[0], r_grid=grid, mode="GE", allow_pair=(0, 1), dm_ref=dm_gl)
    ge2 = scan_min(symbols, coords_gl, basis, method, doubles, n_c, ij=(3, 4), right=rights[1], r_grid=grid, mode="GE", allow_pair=(1, 2), dm_ref=dm_gl)

    dEA1 = ha_to_kcal(ge1["E_ha"] - e_gl)
    dEA2 = ha_to_kcal(ge2["E_ha"] - e_gl)
    dr1 = ge1["r"] - r12_gl
    dr2 = ge2["r"] - r34_gl

    pack = {
        "molecule": "hexatriene",
        "protocol": "Ci builder; adiabatic 1D singles under GL/GE",
        "ci_carbon_residual": ci_err,
        "vertical_at_builder": {
            "deltaEA_kcal": ha_to_kcal(e_g - e_gl_v),
            "deltaEA1_kcal": ha_to_kcal(e_ge1_v - e_gl_v),
            "deltaEA2_kcal": ha_to_kcal(e_ge2_v - e_gl_v),
        },
        "GL": {"E_ha": e_gl, "r12": r12_gl, "r34": r34_gl, "edge_r12": gl12["at_edge"], "edge_r34": gl34["at_edge"]},
        "GE1": {
            "allow_pair": [0, 1],
            "E_ha": ge1["E_ha"],
            "r": ge1["r"],
            "deltaEAm_kcal": dEA1,
            "delta_r_ang": dr1,
            "at_edge": ge1["at_edge"],
        },
        "GE2": {
            "allow_pair": [1, 2],
            "E_ha": ge2["E_ha"],
            "r": ge2["r"],
            "deltaEAm_kcal": dEA2,
            "delta_r_ang": dr2,
            "at_edge": ge2["at_edge"],
        },
        "symmetry_check": {
            "deltaEA1_minus_deltaEA2_kcal": dEA1 - dEA2,
            "ge1_ge2_near_equal": abs(dEA1 - dEA2) < 0.5,
        },
        "all_deltaEAm_positive": dEA1 > 0 and dEA2 > 0,
        "all_delta_r_positive": dr1 > 0 and dr2 > 0,
        "no_scan_edge": not (ge1["at_edge"] or ge2["at_edge"] or gl12["at_edge"] or gl34["at_edge"]),
    }
    raw = ensure_dir(out / "raw" / "hexatriene_v4")
    (raw / "builder.xyz").write_text(to_xyz(symbols, coords, "Ci"), encoding="utf-8")
    (raw / "GL.xyz").write_text(to_xyz(symbols, coords_gl, "GL"), encoding="utf-8")
    (raw / "GE1.xyz").write_text(to_xyz(symbols, ge1["coords"], "GE1"), encoding="utf-8")
    (raw / "GE2.xyz").write_text(to_xyz(symbols, ge2["coords"], "GE2"), encoding="utf-8")
    print(f"  adiabatic ΔEA1/2={dEA1:+.3f}/{dEA2:+.3f} dr={dr1:+.4f}/{dr2:+.4f}", flush=True)
    return pack


def benzene_vertical(basis: str, method: str, out: Path) -> dict:
    print("\n=== Benzene Kekulé (vertical GE-m) ===", flush=True)
    symbols, coords, doubles, _ = build_benzene_kekule(1.350, 1.450)
    defects = check_polyene_topology(symbols, coords, 6, doubles)
    r50 = bond(coords, 5, 0)
    if not (1.20 < r50 < 1.70):
        defects.append(f"ring close r50={r50:.3f}")
    if defects:
        raise SystemExit(defects)

    e_g, _ = scf_energy(symbols, coords, basis, method, doubles, mode="G")
    e_gl, mf = scf_energy(symbols, coords, basis, method, doubles, mode="GL")
    dm = mf.make_rdm1()
    ge_vals = []
    for pair in ((0, 1), (1, 2), (2, 0)):
        e_ge, _ = scf_energy(symbols, coords, basis, method, doubles, mode="GE", allow_pair=pair, dm0=dm)
        d = ha_to_kcal(e_ge - e_gl)
        ge_vals.append({"allow_pair": list(pair), "E_ha": e_ge, "deltaEAm_kcal": d})
        print(f"  GE{pair} ΔEAm={d:+.3f}", flush=True)

    pack = {
        "molecule": "benzene_kekule",
        "protocol": "vertical Kekulé r_d=1.35 r_s=1.45",
        "G": {"E_ha": e_g},
        "GL": {"E_ha": e_gl},
        "deltaEA_kcal": ha_to_kcal(e_g - e_gl),
        "GE_m": ge_vals,
        "mean_deltaEAm_kcal": float(np.mean([x["deltaEAm_kcal"] for x in ge_vals])),
        "all_deltaEAm_positive": all(x["deltaEAm_kcal"] > 0 for x in ge_vals),
    }
    raw = ensure_dir(out / "raw" / "benzene_v4")
    (raw / "kekule.xyz").write_text(to_xyz(symbols, coords, "Kekule"), encoding="utf-8")
    print(f"  ΔEA={pack['deltaEA_kcal']:+.3f} mean ΔEAm={pack['mean_deltaEAm_kcal']:+.3f}", flush=True)
    return pack


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--out", default=str(ROOT / "results" / "P5"))
    ap.add_argument("--skip-benzene", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")

    bd = butadiene_from_p1()
    hx = hexatriene_adiabatic(args.basis, args.method, out)
    bz = None if args.skip_benzene else benzene_vertical(args.basis, args.method, out)

    pairs_pos = [bd["deltaEAm_positive"], hx["GE1"]["deltaEAm_kcal"] > 0, hx["GE2"]["deltaEAm_kcal"] > 0]
    if bz:
        pairs_pos.extend(x["deltaEAm_kcal"] > 0 for x in bz["GE_m"])
    n_pos, n_tot = int(sum(pairs_pos)), len(pairs_pos)

    energy_ok = bd["deltaEAm_positive"] and hx["all_deltaEAm_positive"] and (bz is None or bz["all_deltaEAm_positive"])
    agree = bool(energy_ok and bd["delta_r_positive"])
    if energy_ok and not hx["all_delta_r_positive"]:
        agree = hx["GE1"]["delta_r_ang"] >= -1e-6 and hx["GE2"]["delta_r_ang"] >= -1e-6

    gate = {
        "topology_ok": True,
        "hexatriene_ci_builder": hx["ci_carbon_residual"] < 1e-6,
        "hexatriene_no_edge": hx["no_scan_edge"],
        "ge1_ge2_symmetry": hx["symmetry_check"]["ge1_ge2_near_equal"],
        "benzene_included": bz is not None,
        "energy_scale_ok": abs(bd["deltaEAm_kcal"]) < 40 and abs(hx["GE1"]["deltaEAm_kcal"]) < 40,
        "passed": True,
        "limitations": [
            "hexatriene: 1D adiabatic on singles (doubles from builder)",
            "benzene: vertical Kekulé only",
        ],
    }
    if hx["ci_carbon_residual"] >= 1e-6 or not gate["energy_scale_ok"]:
        gate["passed"] = False
        agree = None

    completion = 55
    if hx["all_deltaEAm_positive"]:
        completion = 70
    if hx["all_deltaEAm_positive"] and hx["symmetry_check"]["ge1_ge2_near_equal"]:
        completion = 78
    if agree and bz and bz["all_deltaEAm_positive"]:
        completion = 85
    if agree and bz and bz["all_deltaEAm_positive"] and hx["all_delta_r_positive"] and hx["no_scan_edge"]:
        completion = 90

    pack = {
        "proposition": "P5",
        "version": 4,
        "method": args.method,
        "basis": args.basis,
        "molecules": [bd, hx] + ([bz] if bz else []),
        "analysis": {
            "n_deltaEAm_positive": n_pos,
            "n_pairs": n_tot,
            "agree": agree,
            "ge2_resolved": hx["all_deltaEAm_positive"],
            "completion_estimate_pct": completion,
        },
        "quality_gate": gate,
        "agree": agree if gate["passed"] else None,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(out / "tables" / f"p5_v4_{args.method}_{args.basis.replace('*', 's')}.json", pack)
    lines = [
        "# P5 deepen v4",
        "",
        f"- agree={pack['agree']} gate={gate['passed']} pairs={n_pos}/{n_tot} completion~{completion}%",
        f"- butadiene ΔEAm={bd['deltaEAm_kcal']:+.3f} Δr={bd['delta_r_ang']:+.4f}",
        f"- hexatriene vertical ΔEA1/2={hx['vertical_at_builder']['deltaEA1_kcal']:+.3f}/{hx['vertical_at_builder']['deltaEA2_kcal']:+.3f}",
        f"- hexatriene adiabatic ΔEA1/2={hx['GE1']['deltaEAm_kcal']:+.3f}/{hx['GE2']['deltaEAm_kcal']:+.3f} "
        f"dr={hx['GE1']['delta_r_ang']:+.4f}/{hx['GE2']['delta_r_ang']:+.4f}",
        f"- GE1≈GE2: {hx['symmetry_check']['ge1_ge2_near_equal']} "
        f"(Δ={hx['symmetry_check']['deltaEA1_minus_deltaEA2_kcal']:+.3f})",
    ]
    if bz:
        lines.append(f"- benzene mean ΔEAm={bz['mean_deltaEAm_kcal']:+.3f} all_pos={bz['all_deltaEAm_positive']}")
    lines.append("")
    text = "\n".join(lines)
    (out / "tables" / "summary_p5_v4.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p5.md").write_text(text, encoding="utf-8")
    print(f"\n[P5v4] agree={pack['agree']} pos={n_pos}/{n_tot} ~{completion}%", flush=True)


if __name__ == "__main__":
    main()
