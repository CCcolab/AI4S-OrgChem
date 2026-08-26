"""
P1 — Butadiene conjugation energy CE1 vs CE2.

Uses SciPy BFGS + PySCF analytic gradients (no geometric/berny install).

  CE1 = ΔH_hyd(butadiene) - 2 ΔH_hyd(1-butene)
  CE2 = ΔH_hyd(butadiene) - 2 ΔH_hyd(trans-2-butene)

with ΔH_hyd defined from electronic energies to n-butane (sign: stabilizing < 0).
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

BOHR = 0.52917721092  # Angstrom per Bohr

GEOMS: dict[str, str] = {
    "h2": """
H 0.000000 0.000000 0.000000
H 0.000000 0.000000 0.740000
""",
    "butadiene": """
C  0.000000  0.000000  0.000000
C  1.340000  0.000000  0.000000
C  2.050000  1.260000  0.000000
C  3.390000  1.260000  0.000000
H -0.540000 -0.940000  0.000000
H -0.540000  0.940000  0.000000
H  1.880000 -0.940000  0.000000
H  1.510000  2.200000  0.000000
H  3.930000  2.200000  0.000000
H  3.930000  0.320000  0.000000
""",
    "1_butene": """
C  0.000000  0.000000  0.000000
C  1.340000  0.000000  0.000000
C  2.000000  1.400000  0.000000
C  3.500000  1.400000  0.000000
H -0.540000 -0.940000  0.000000
H -0.540000  0.940000  0.000000
H  1.880000 -0.940000  0.000000
H  1.460000  1.940000  0.890000
H  1.460000  1.940000 -0.890000
H  3.900000  2.300000  0.000000
H  3.900000  0.860000  0.890000
H  3.900000  0.860000 -0.890000
""",
    "trans_2_butene": """
C  0.000000  0.000000  0.000000
C  1.500000  0.000000  0.000000
C  2.200000  1.260000  0.000000
C  3.700000  1.260000  0.000000
H -0.400000  0.900000  0.890000
H -0.400000  0.900000 -0.890000
H -0.400000 -1.000000  0.000000
H  1.850000 -0.940000  0.000000
H  1.850000  2.200000  0.000000
H  4.100000  2.160000  0.890000
H  4.100000  2.160000 -0.890000
H  4.100000  0.260000  0.000000
""",
    "n_butane": """
C  0.000000  0.000000  0.000000
C  1.540000  0.000000  0.000000
C  2.090000  1.430000  0.000000
C  3.630000  1.430000  0.000000
H -0.360000  0.900000  0.890000
H -0.360000  0.900000 -0.890000
H -0.360000 -1.000000  0.000000
H  1.900000 -0.900000  0.890000
H  1.900000 -0.900000 -0.890000
H  1.730000  1.950000  0.890000
H  1.730000  1.950000 -0.890000
H  4.000000  2.330000  0.000000
H  4.000000  0.890000  0.890000
H  4.000000  0.890000 -0.890000
""",
}


def make_mol(xyz: str, basis: str) -> gto.Mole:
    mol = gto.Mole()
    mol.atom = xyz
    mol.basis = basis
    mol.unit = "Angstrom"
    mol.verbose = 0
    mol.build()
    return mol


def xyz_of(mol: gto.Mole) -> str:
    c = mol.atom_coords(unit="Angstrom")
    return "\n".join(
        f"{mol.atom_symbol(i)} {c[i,0]:.8f} {c[i,1]:.8f} {c[i,2]:.8f}"
        for i in range(mol.natm)
    )


def make_mf(mol: gto.Mole, method: str):
    if method == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
        return mf
    return scf.RHF(mol)


def energy_and_grad(mol: gto.Mole, method: str) -> tuple[float, np.ndarray]:
    """Energy (Ha) and gradient dE/dR in Ha/Bohr, shape (natm,3)."""
    mf = make_mf(mol, method)
    mf.kernel()
    # Prefer nuc_grad_method() so DFT/HF both resolve the correct Gradients class.
    g = mf.nuc_grad_method().kernel()
    return float(mf.e_tot), np.asarray(g, dtype=float)


def load_xyz_body(path: Path) -> str | None:
    """Parse XYZ file (natoms / comment / coords) → coord block only."""
    if not path.is_file():
        return None
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    if len(lines) < 3:
        return None
    try:
        n = int(lines[0].strip())
    except ValueError:
        return "\n".join(lines)  # bare coord block
    return "\n".join(lines[2 : 2 + n])


def starting_xyz(name: str, default: str, start_dir: Path | None) -> str:
    if start_dir is None:
        return default
    body = load_xyz_body(start_dir / f"{name}.xyz")
    return body if body else default


def optimize(xyz: str, basis: str, method: str, maxiter: int = 80) -> tuple[float, str]:
    """BFGS Cartesian optimization using SciPy + PySCF gradients."""
    mol0 = make_mol(xyz, basis)
    # Work in Bohr for consistency with PySCF gradient units
    x0 = mol0.atom_coords(unit="Bohr").ravel().copy()
    symbols = [mol0.atom_symbol(i) for i in range(mol0.natm)]

    def pack_mol(x: np.ndarray) -> gto.Mole:
        coords = x.reshape(-1, 3)
        atom = []
        for s, r in zip(symbols, coords):
            # PySCF default internal unit Bohr when unit='Bohr'
            atom.append([s, (float(r[0]), float(r[1]), float(r[2]))])
        mol = gto.Mole()
        mol.atom = atom
        mol.basis = basis
        mol.unit = "Bohr"
        mol.verbose = 0
        mol.build()
        return mol

    cache: dict[str, tuple[float, np.ndarray]] = {}

    def fun(x: np.ndarray) -> float:
        key = x.tobytes()
        if key not in cache:
            e, g = energy_and_grad(pack_mol(x), method if method != "MP2" else "RHF")
            cache[key] = (e, g.ravel())
        return cache[key][0]

    def jac(x: np.ndarray) -> np.ndarray:
        key = x.tobytes()
        if key not in cache:
            fun(x)
        return cache[key][1]

    # Looser gtol for DFT/HF cart opt of organics
    res = minimize(
        fun,
        x0,
        method="BFGS",
        jac=jac,
        options={"maxiter": maxiter, "gtol": 1.0e-4, "disp": False},
    )
    mol_eq = pack_mol(res.x)
    if method == "MP2":
        from pyscf import mp

        mf = scf.RHF(mol_eq).run()
        e = float(mp.MP2(mf).run().e_tot)
    else:
        e, _ = energy_and_grad(mol_eq, method)
    return e, xyz_of(mol_eq)


def run_level(
    basis: str,
    method: str,
    out: Path,
    *,
    start_dir: Path | None = None,
    maxiter: int = 80,
) -> dict:
    tag = f"{method}_{basis.replace('*', 's').replace('+', 'p')}"
    level_dir = ensure_dir(out / "raw" / tag)
    energies: dict[str, float] = {}

    for name, geom0 in GEOMS.items():
        geom = starting_xyz(name, geom0, start_dir)
        print(f"[P1] {tag} :: {name}", flush=True)
        if method == "MP2":
            # HF geometry, MP2 single-point energy
            _, xyz = optimize(geom, basis, "RHF", maxiter=maxiter)
            mol = make_mol(xyz, basis)
            from pyscf import mp

            mf = scf.RHF(mol).run()
            e = float(mp.MP2(mf).run().e_tot)
        else:
            e, xyz = optimize(geom, basis, method, maxiter=maxiter)
        energies[name] = e
        n = len(xyz.strip().splitlines())
        (level_dir / f"{name}.xyz").write_text(
            f"{n}\n{name} {tag}\n{xyz}\n", encoding="utf-8"
        )

    e_bd, e_h2 = energies["butadiene"], energies["h2"]
    e_1b, e_t2, e_bu = energies["1_butene"], energies["trans_2_butene"], energies["n_butane"]

    dh_bd = e_bu - e_bd - 2 * e_h2
    dh_1b = e_bu - e_1b - e_h2
    dh_t2 = e_bu - e_t2 - e_h2
    ce1 = dh_bd - 2 * dh_1b
    ce2 = dh_bd - 2 * dh_t2

    row = {
        "tag": tag,
        "method": method,
        "basis": basis,
        "optimizer": "scipy.BFGS + PySCF gradients (no geometric/berny)",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "energies_ha": energies,
        "deltaH_hyd_kcal": {
            "butadiene": ha_to_kcal(dh_bd),
            "1_butene": ha_to_kcal(dh_1b),
            "trans_2_butene": ha_to_kcal(dh_t2),
        },
        "CE1_kcal": ha_to_kcal(ce1),
        "CE2_kcal": ha_to_kcal(ce2),
        "sign_CE1": "stabilizing" if ce1 < 0 else "destabilizing",
        "sign_CE2": "stabilizing" if ce2 < 0 else "destabilizing",
        "hartree_to_kcal": HARTREE_TO_KCAL,
        "book_refs_kcal": {"classic_CE1": -3.9, "CE2_trans2": 1.9, "GL": 1.4},
        "agree_criteria": "CE1 < 0 and CE2 > 0",
        "sign_flip": bool(ce1 < 0 and ce2 > 0),
        "start_geom_dir": str(start_dir) if start_dir else None,
    }
    write_json(out / "tables" / f"ce_{tag}.json", row)
    return row


def rewrite_csv(path: Path, rows: list[dict]) -> None:
    ensure_dir(path.parent)
    fields = [
        "tag",
        "method",
        "basis",
        "CE1_kcal",
        "CE2_kcal",
        "sign_CE1",
        "sign_CE2",
        "sign_flip",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(
                {
                    "tag": row["tag"],
                    "method": row["method"],
                    "basis": row["basis"],
                    "CE1_kcal": f"{row['CE1_kcal']:.4f}",
                    "CE2_kcal": f"{row['CE2_kcal']:.4f}",
                    "sign_CE1": row["sign_CE1"],
                    "sign_CE2": row["sign_CE2"],
                    "sign_flip": "yes" if row.get("sign_flip") else "no",
                }
            )


def write_sensitivity_pack(out: Path, rows: list[dict], basis: str) -> None:
    """Cross-method comparison for L1/L2 robustness."""
    any_flip = any(r.get("sign_flip") for r in rows)
    all_ce1_pos = all(r["CE1_kcal"] > 0 for r in rows)
    all_ce2_pos = all(r["CE2_kcal"] > 0 for r in rows)
    pack = {
        "proposition": "P1",
        "kind": "method_sensitivity",
        "basis": basis,
        "agree_criteria": "CE1 < 0 and CE2 > 0",
        "methods": [
            {
                "method": r["method"],
                "tag": r["tag"],
                "CE1_kcal": r["CE1_kcal"],
                "CE2_kcal": r["CE2_kcal"],
                "sign_flip": r.get("sign_flip", False),
            }
            for r in rows
        ],
        "any_sign_flip": any_flip,
        "all_CE1_positive": all_ce1_pos,
        "all_CE2_positive": all_ce2_pos,
        "conclusion_zh": (
            "敏感性未改变主结论：各层次均未出现 CE1<0 且 CE2>0 的符号翻转。"
            if not any_flip
            else "至少一层出现符号翻转；须复核 VERDICT 是否改判。"
        ),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json(out / "tables" / "sensitivity_methods.json", pack)

    lines = [
        "# P1 method sensitivity",
        "",
        f"- Basis: `{basis}`",
        f"- Agree if CE1<0 and CE2>0",
        f"- Any sign flip across methods: **{'yes' if any_flip else 'no'}**",
        "",
        "| method | CE1 | CE2 | flip |",
        "|--------|-----|-----|------|",
    ]
    for r in rows:
        flip = "yes" if r.get("sign_flip") else "no"
        lines.append(
            f"| {r['method']} | {r['CE1_kcal']:+.3f} | {r['CE2_kcal']:+.3f} | {flip} |"
        )
    lines.extend(["", pack["conclusion_zh"], ""])
    (out / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines), flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP", choices=["B3LYP", "RHF", "MP2"])
    ap.add_argument("--out", default=str(ROOT / "results" / "P1"))
    ap.add_argument(
        "--sensitivity",
        action="store_true",
        help="Run RHF + MP2 (reuse B3LYP start geoms if present); merge with existing B3LYP JSON",
    )
    ap.add_argument(
        "--start-from",
        default="",
        help="Directory of .xyz start geometries (default: results/.../raw/B3LYP_* when sensitivity)",
    )
    ap.add_argument("--maxiter", type=int, default=80)
    args = ap.parse_args()

    out = Path(args.out)
    ensure_dir(out / "tables")
    ensure_dir(out / "raw")

    b3_tag = f"B3LYP_{args.basis.replace('*', 's').replace('+', 'p')}"
    default_start = out / "raw" / b3_tag
    start_dir = Path(args.start_from) if args.start_from else (
        default_start if args.sensitivity and default_start.is_dir() else None
    )

    if args.sensitivity:
        methods = ["RHF", "MP2"]
        new_rows = [
            run_level(
                args.basis, m, out, start_dir=start_dir, maxiter=args.maxiter
            )
            for m in methods
        ]
        # Merge with prior B3LYP if available
        b3_path = out / "tables" / f"ce_{b3_tag}.json"
        rows: list[dict] = []
        if b3_path.is_file():
            from src.common.units import read_json

            b3 = read_json(b3_path)
            b3.setdefault("sign_flip", bool(b3["CE1_kcal"] < 0 and b3["CE2_kcal"] > 0))
            rows.append(b3)
        rows.extend(new_rows)
        # stable order: RHF, B3LYP, MP2 when all present
        order = {"RHF": 0, "B3LYP": 1, "MP2": 2}
        rows.sort(key=lambda r: order.get(r["method"], 9))
        rewrite_csv(out / "tables" / "ce_summary.csv", rows)
        write_sensitivity_pack(out, rows, args.basis)
        write_json(
            out / "meta.json",
            {
                "proposition": "P1",
                "updated_utc": datetime.now(timezone.utc).isoformat(),
                "latest_tag": "sensitivity_methods",
                "methods": [r["tag"] for r in rows],
            },
        )
        return

    row = run_level(
        args.basis, args.method, out, start_dir=start_dir, maxiter=args.maxiter
    )
    # Keep multi-method CSV if present: update/replace this method only
    existing: list[dict] = []
    csv_path = out / "tables" / "ce_summary.csv"
    if csv_path.is_file():
        with csv_path.open(encoding="utf-8") as f:
            for rec in csv.DictReader(f):
                if rec["method"] == row["method"]:
                    continue
                existing.append(
                    {
                        "tag": rec["tag"],
                        "method": rec["method"],
                        "basis": rec["basis"],
                        "CE1_kcal": float(rec["CE1_kcal"]),
                        "CE2_kcal": float(rec["CE2_kcal"]),
                        "sign_CE1": rec["sign_CE1"],
                        "sign_CE2": rec["sign_CE2"],
                        "sign_flip": rec.get("sign_flip", "no") == "yes",
                    }
                )
    existing.append(row)
    rewrite_csv(csv_path, existing)
    write_json(
        out / "meta.json",
        {"proposition": "P1", "updated_utc": row["timestamp_utc"], "latest_tag": row["tag"]},
    )
    md = "\n".join(
        [
            f"# P1 summary ({row['tag']})",
            "",
            f"- CE1 (vs 1-butene) = **{row['CE1_kcal']:.3f} kcal/mol** ({row['sign_CE1']})",
            f"- CE2 (vs trans-2-butene) = **{row['CE2_kcal']:.3f} kcal/mol** ({row['sign_CE2']})",
            f"- Sign flip: **{'yes' if row['sign_flip'] else 'no'}**",
            "- Book: classic≈−3.9; CE2≈+1.9; GL≈+1.4 kcal/mol",
            "- Agree if CE1<0 and CE2>0",
            "",
        ]
    )
    (out / "summary.md").write_text(md, encoding="utf-8")
    print(md, flush=True)


if __name__ == "__main__":
    main()

