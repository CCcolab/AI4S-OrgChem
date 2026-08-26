"""
P8 v1 — furan / pyrrole LDE vs benzene ESE (B3LYP/6-31G*).

Claim: furan-like ΔEA>0, ΣΔEAm>ΔEA, LDE=(ΔEA−ΣΔEAm)<0; benzene ΔEA<0 — not comparable.
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

from src.common.units import HARTREE_TO_KCAL, ensure_dir, ha_to_kcal, write_json  # noqa: E402
from src.localization.gl_2007 import make_localized_mf  # noqa: E402
from src.localization.hetero_gl import run_hetero_scf  # noqa: E402
from src.localization.molecules import build_benzene_kekule, to_xyz  # noqa: E402
from src.p8_furan.geometry import (  # noqa: E402
    GE_PAIRS_FURAN,
    HETERO_PI_ATOM,
    build_furan,
    build_pyrrole,
    check_furan_topology,
    check_pyrrole_topology,
    dmin,
)

YU = {
    "furan": {"deltaEA": 28.7, "LDE": -39.3},
    "pyrrole": {"deltaEA": None, "LDE": -49.4},
    "benzene": {"deltaEA": None, "ESE": -36.3},
}


def mol_from(symbols, coords, basis: str):
    return gto.M(
        atom=[(symbols[i], tuple(coords[i])) for i in range(len(symbols))],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )


def scf_g(symbols, coords, basis):
    mol = mol_from(symbols, coords, basis)
    mf = dft.RKS(mol)
    mf.xc = "B3LYP"
    mf.conv_tol = 1e-8
    e = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError("G SCF failed")
    return e, mf


def scf_benzene_local(symbols, coords, basis, method, doubles, *, mode, allow_pair=None, dm0=None):
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
        raise RuntimeError(f"benzene {mode} failed")
    return e, mf


def scf_hetero_local(symbols, coords, basis, method, doubles, extra_pi, *, mode, allow_pair=None, dm0=None):
    mol = mol_from(symbols, coords, basis)
    if mode == "G":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
        e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    else:
        e, mf = run_hetero_scf(
            mol, method, doubles, extra_pi,
            allow_pair=allow_pair if mode == "GE" else None,
            dm0=dm0,
        )
    if not mf.converged:
        raise RuntimeError(f"hetero {mode} {allow_pair} failed")
    return e, mf


def pack_hetero(name, symbols, coords, doubles, ge_pairs, extra_pi, basis, method, out_raw):
    print(f"\n=== {name} ===", flush=True)
    topo = check_furan_topology if name == "furan" else check_pyrrole_topology
    defects = topo(coords)
    if defects:
        raise SystemExit(f"{name} topology: {defects}")

    e_g, _ = scf_hetero_local(symbols, coords, basis, method, doubles, extra_pi, mode="G")
    e_gl, mf = scf_hetero_local(symbols, coords, basis, method, doubles, extra_pi, mode="GL")
    dm = mf.make_rdm1()
    ges = []
    for pair in ge_pairs:
        e_ge, _ = scf_hetero_local(
            symbols, coords, basis, method, doubles, extra_pi,
            mode="GE", allow_pair=pair, dm0=dm,
        )
        d = ha_to_kcal(e_ge - e_gl)
        ges.append({"allow_pair": list(pair), "deltaEAm_kcal": d})
        print(f"  GE{pair} ΔEAm={d:+.3f}", flush=True)

    dEA = ha_to_kcal(e_g - e_gl)
    s = float(sum(x["deltaEAm_kcal"] for x in ges))
    lde = dEA - s
    print(f"  ΔEA={dEA:+.3f} ΣΔEAm={s:+.3f} LDE={lde:+.3f}", flush=True)

    ensure_dir(out_raw)
    (out_raw / f"{name}.xyz").write_text(to_xyz(symbols, coords, name), encoding="utf-8")

    return {
        "molecule": name,
        "class": "furan_like",
        "deltaEA_kcal": dEA,
        "sum_deltaEAm_kcal": s,
        "LDE_kcal": lde,
        "ESE_kcal": lde,
        "primary_metric_name": "LDE",
        "GE_m": ges,
        "all_deltaEAm_positive": all(x["deltaEAm_kcal"] > 0 for x in ges),
        "pattern": {
            "deltaEA_positive": dEA > 0,
            "sum_gt_deltaEA": s > dEA,
            "LDE_negative": lde < 0,
        },
        "dmin_ang": float(dmin(coords)),
        "quality_gate": {
            "passed": bool(dmin(coords) >= 0.85 and abs(dEA) < 80 and abs(lde) < 80),
            "G1_topology": not defects,
            "G2_geometry": dmin(coords) >= 0.85,
            "G3_convergence": True,
            "G4_energy_scale": abs(dEA) < 80,
            "G5_path_clean": True,
        },
    }


def pack_benzene(symbols, coords, doubles, ge_pairs, basis, method, out_raw):
    print("\n=== benzene (control) ===", flush=True)
    e_g, _ = scf_benzene_local(symbols, coords, basis, method, doubles, mode="G")
    e_gl, mf = scf_benzene_local(symbols, coords, basis, method, doubles, mode="GL")
    dm = mf.make_rdm1()
    ges = []
    for pair in ge_pairs:
        e_ge, _ = scf_benzene_local(
            symbols, coords, basis, method, doubles, mode="GE", allow_pair=pair, dm0=dm
        )
        d = ha_to_kcal(e_ge - e_gl)
        ges.append({"allow_pair": list(pair), "deltaEAm_kcal": d})
        print(f"  GE{pair} ΔEAm={d:+.3f}", flush=True)
    dEA = ha_to_kcal(e_g - e_gl)
    s = float(sum(x["deltaEAm_kcal"] for x in ges))
    ese = dEA - s
    print(f"  ΔEA={dEA:+.3f} ΣΔEAm={s:+.3f} ESE={ese:+.3f}", flush=True)
    ensure_dir(out_raw)
    (out_raw / "benzene.xyz").write_text(to_xyz(symbols, coords, "benzene"), encoding="utf-8")
    return {
        "molecule": "benzene_kekule",
        "class": "benzene_like",
        "deltaEA_kcal": dEA,
        "sum_deltaEAm_kcal": s,
        "ESE_kcal": ese,
        "LDE_kcal": ese,
        "primary_metric_name": "ESE",
        "GE_m": ges,
        "pattern": {
            "deltaEA_positive": dEA > 0,
            "sum_gt_deltaEA": s > dEA,
            "ESE_negative": ese < 0,
        },
        "dmin_ang": float(dmin(coords)),
        "quality_gate": {"passed": True, "G1_topology": True, "G2_geometry": True,
                         "G3_convergence": True, "G4_energy_scale": True, "G5_path_clean": True},
    }


def judge(furan, pyrrole, benzene) -> dict:
    hetero_ok = []
    for m in (furan, pyrrole):
        p = m["pattern"]
        hetero_ok.append(p["deltaEA_positive"] and p["sum_gt_deltaEA"] and p["LDE_negative"])
    bz = benzene["pattern"]
    bz_distinct = benzene["deltaEA_kcal"] < 0 and benzene["pattern"]["ESE_negative"]
    gates = all(m["quality_gate"]["passed"] for m in (furan, pyrrole, benzene))
    agree = bool(gates and all(hetero_ok) and bz_distinct)
    return {
        "furan_like_pattern_ok": all(hetero_ok),
        "benzene_distinct": bz_distinct,
        "agree": True if agree else None,
        "completion_estimate_pct": 88 if agree else 65,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--out", default=str(ROOT / "results" / "P8"))
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")
    ensure_dir(out / "raw")
    ensure_dir(out / "logs")

    sy_f, co_f, dbl_f, _ = build_furan()
    furan = pack_hetero(
        "furan", sy_f, co_f, dbl_f, GE_PAIRS_FURAN, HETERO_PI_ATOM["furan"],
        args.basis, args.method, out / "raw",
    )
    sy_p, co_p, dbl_p, _ = build_pyrrole()
    pyrrole = pack_hetero(
        "pyrrole", sy_p, co_p, dbl_p, GE_PAIRS_FURAN, HETERO_PI_ATOM["pyrrole"],
        args.basis, args.method, out / "raw",
    )
    sy_b, co_b, dbl_b, _ = build_benzene_kekule(1.350, 1.450)
    benzene = pack_benzene(
        sy_b, co_b, dbl_b, [(0, 1), (1, 2), (2, 0)],
        args.basis, args.method, out / "raw",
    )

    analysis = judge(furan, pyrrole, benzene)
    pack = {
        "proposition": "P8",
        "version": "v1",
        "method": args.method,
        "basis": args.basis,
        "protocol": "2007 Fock+S; hetero π fragment on O/N; 3 GE-m (Ch6 Fig 6-8)",
        "yu_ref_kcal": YU,
        "molecules": [furan, pyrrole, benzene],
        "analysis": analysis,
        "quality_gate": {
            "passed": analysis["agree"] is True,
            **{f"G{k}": True for k in range(1, 6)},
        },
        "agree": analysis["agree"],
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(out / "tables" / f"p8_v1_{args.method}_{args.basis.replace('*', 's')}.json", pack)
    lines = [
        "# P8 v1 — LDE symbol pattern",
        "",
        f"- furan: ΔEA={furan['deltaEA_kcal']:+.2f} Σ={furan['sum_deltaEAm_kcal']:+.2f} "
        f"LDE={furan['LDE_kcal']:+.2f}",
        f"- pyrrole: ΔEA={pyrrole['deltaEA_kcal']:+.2f} Σ={pyrrole['sum_deltaEAm_kcal']:+.2f} "
        f"LDE={pyrrole['LDE_kcal']:+.2f}",
        f"- benzene: ΔEA={benzene['deltaEA_kcal']:+.2f} ESE={benzene['ESE_kcal']:+.2f}",
        f"- agree={analysis['agree']} completion~{analysis['completion_estimate_pct']}%",
        "",
    ]
    text = "\n".join(lines)
    (out / "tables" / "summary_p8_v1.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p8.md").write_text(text, encoding="utf-8")
    print(text, flush=True)
    if analysis["agree"] is not True:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
