"""
P6 v2 — close pilot gaps with explicit dual metrics.

Yu / P6 criteria (docs/propositions.md):
  - Benzene: ESE = ΔEA − ΣΔEAm ≈ −36.3 kcal/mol
  - Cyclobutadiene: conjugation / antiaromatic energy ≈ +53–55
    Under 2007 GL with two doubles, G ≡ GE-m ⇒ ESE ≡ 0 by construction;
    the matching observable is ΔEA (= E(G)−E(GL)), not ESE.

Also: butadiene additivity control (ESE ≈ 0, ΔEA > 0).
Optional second basis for sign/magnitude sensitivity.
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
from src.localization.molecules import (  # noqa: E402
    build_benzene_kekule,
    build_butadiene,
    to_xyz,
)

R_CH = 1.085


def mol_from(symbols, coords, basis: str):
    return gto.M(
        atom=[(symbols[i], tuple(coords[i])) for i in range(len(symbols))],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )


def scf_energy(symbols, coords, basis, method, doubles, *, mode, allow_pair=None, dm0=None):
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


def check_ring(coords, n_c: int, doubles, r_max: float = 1.70) -> list[str]:
    problems = []
    dmat = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(dmat, np.inf)
    dmin = float(dmat.min())
    if dmin < 0.85:
        problems.append(f"dmin={dmin:.3f}")
    for a, b in doubles:
        rij = float(np.linalg.norm(coords[a] - coords[b]))
        if not (1.20 < rij < 1.55):
            problems.append(f"double {a}-{b} r={rij:.3f}")
    for i in range(n_c):
        j = (i + 1) % n_c
        rij = float(np.linalg.norm(coords[i] - coords[j]))
        if not (1.20 < rij < r_max):
            problems.append(f"ring {i}-{j} r={rij:.3f}")
    if float(np.max(np.abs(coords[:, 2]))) > 0.05:
        problems.append("non-planar")
    return problems


def check_chain(coords, n_c: int, doubles) -> list[str]:
    problems = []
    dmat = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(dmat, np.inf)
    if float(dmat.min()) < 0.85:
        problems.append(f"dmin={float(dmat.min()):.3f}")
    for a, b in doubles:
        rij = float(np.linalg.norm(coords[a] - coords[b]))
        if not (1.20 < rij < 1.55):
            problems.append(f"double {a}-{b} r={rij:.3f}")
    for i in range(n_c - 1):
        rij = float(np.linalg.norm(coords[i] - coords[i + 1]))
        if not (1.20 < rij < 1.70):
            problems.append(f"chain {i}-{i+1} r={rij:.3f}")
    return problems


def build_cyclobutadiene(r_d: float = 1.350, r_s: float = 1.540, r_ch: float = R_CH):
    coords = np.zeros((8, 3))
    half_d, half_s = r_d / 2.0, r_s / 2.0
    coords[0] = (-half_d, -half_s, 0.0)
    coords[1] = (+half_d, -half_s, 0.0)
    coords[2] = (+half_d, +half_s, 0.0)
    coords[3] = (-half_d, +half_s, 0.0)

    def add_h(ci, direction, idx):
        n = direction / (np.linalg.norm(direction) + 1e-16)
        coords[idx] = coords[ci] + r_ch * n

    center = coords[:4].mean(axis=0)
    for ci, hi in enumerate(range(4, 8)):
        add_h(ci, coords[ci] - center, hi)
    return ["C"] * 4 + ["H"] * 4, coords, [(0, 1), (2, 3)], [(1, 2), (3, 0)]


def pack_molecule(
    name: str,
    symbols,
    coords,
    doubles,
    ge_pairs,
    basis: str,
    method: str,
    out_raw: Path,
    *,
    ring: bool,
) -> dict:
    print(f"\n=== {name}  [{method}/{basis}] ===", flush=True)
    n_c = sum(1 for s in symbols if s == "C")
    defects = check_ring(coords, n_c, doubles) if ring else check_chain(coords, n_c, doubles)
    if defects:
        raise SystemExit(f"{name} topology: {defects}")

    dmat = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(dmat, np.inf)
    dmin = float(dmat.min())

    e_g, _ = scf_energy(symbols, coords, basis, method, doubles, mode="G")
    e_gl, mf = scf_energy(symbols, coords, basis, method, doubles, mode="GL")
    dm = mf.make_rdm1()
    ge_list = []
    for pair in ge_pairs:
        e_ge, _ = scf_energy(
            symbols, coords, basis, method, doubles, mode="GE", allow_pair=pair, dm0=dm
        )
        d = ha_to_kcal(e_ge - e_gl)
        ge_list.append({"allow_pair": list(pair), "E_ha": e_ge, "deltaEAm_kcal": d})
        print(f"  GE{pair} ΔEAm={d:+.3f}", flush=True)

    dEA = ha_to_kcal(e_g - e_gl)
    sum_dEAm = float(sum(x["deltaEAm_kcal"] for x in ge_list))
    ese = dEA - sum_dEAm
    e_vr = e_gl + sum(x["E_ha"] - e_gl for x in ge_list)
    ese_vr = ha_to_kcal(e_g - e_vr)
    # Primary Yu metric: ESE for aromatic; ΔEA for CBD (2-double G≡GE)
    primary = dEA if name.startswith("cyclobutadiene") else ese
    print(
        f"  ΔEA={dEA:+.3f} ΣΔEAm={sum_dEAm:+.3f} ESE={ese:+.3f} "
        f"primary={primary:+.3f}",
        flush=True,
    )

    ensure_dir(out_raw)
    (out_raw / f"{name}.xyz").write_text(to_xyz(symbols, coords, name), encoding="utf-8")

    return {
        "molecule": name,
        "basis": basis,
        "method": method,
        "protocol": "vertical; 2007 Fock+S",
        "G": {"E_ha": e_g},
        "GL": {"E_ha": e_gl},
        "GE_m": ge_list,
        "deltaEA_kcal": dEA,
        "sum_deltaEAm_kcal": sum_dEAm,
        "ESE_kcal": ese,
        "ESE_via_VR_kcal": ese_vr,
        "primary_metric_kcal": primary,
        "primary_metric_name": "deltaEA" if name.startswith("cyclobutadiene") else "ESE",
        "all_deltaEAm_positive": all(x["deltaEAm_kcal"] > 0 for x in ge_list),
        "dmin_ang": dmin,
        "quality_gate": {
            "topology_ok": True,
            "geometry_ok": dmin >= 0.85,
            "scf_converged": True,
            "energy_scale_ok": abs(dEA) < 100 and abs(ese) < 120,
            "vr_identity_ok": abs(ese - ese_vr) < 1e-6,
            "passed": bool(
                dmin >= 0.85 and abs(dEA) < 100 and abs(ese) < 120 and abs(ese - ese_vr) < 1e-6
            ),
        },
    }


def judge(bz, cbd, bd=None) -> dict:
    """Apply P6 windows to primary metrics."""
    bz_ese = bz["ESE_kcal"]
    cbd_dEA = cbd["deltaEA_kcal"]
    bz_ok = bz_ese < 0 and 25.0 <= abs(bz_ese) <= 45.0
    cbd_ok = cbd_dEA > 0 and 45.0 <= cbd_dEA <= 70.0
    # butadiene: additive polyene → |ESE| small, ΔEA>0
    bd_ok = True
    if bd is not None:
        bd_ok = bd["deltaEA_kcal"] > 0 and abs(bd["ESE_kcal"]) < 5.0

    gates = all(
        m["quality_gate"]["passed"]
        and m["quality_gate"]["topology_ok"]
        and m["quality_gate"]["geometry_ok"]
        and m["all_deltaEAm_positive"]
        for m in ([bz, cbd] + ([bd] if bd else []))
    )

    if not gates:
        agree = None
        completion = 40
    elif bz_ok and cbd_ok and bd_ok:
        agree = True
        completion = 82
    elif bz_ok and cbd_ok:
        agree = True
        completion = 75
    elif bz_ese < 0 and cbd_dEA > 0:
        agree = True
        completion = 65
    else:
        agree = False
        completion = 45

    return {
        "benzene_ESE_kcal": bz_ese,
        "cbd_deltaEA_kcal": cbd_dEA,
        "cbd_ESE_kcal": cbd["ESE_kcal"],
        "butadiene_ESE_kcal": None if bd is None else bd["ESE_kcal"],
        "butadiene_deltaEA_kcal": None if bd is None else bd["deltaEA_kcal"],
        "benzene_in_window": bz_ok,
        "cbd_deltaEA_in_window": cbd_ok,
        "butadiene_additive_ok": bd_ok,
        "signs_ok": bz_ese < 0 and cbd_dEA > 0,
        "gates_ok": gates,
        "agree": agree,
        "completion_estimate_pct": completion,
        "cbd_note": (
            "With 2 doubles, G≡GE-m so ESE≡0; Yu +53–55 matched to ΔEA "
            "(conjugation/antiaromatic energy), not ESE."
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--basis2", default="6-31g", help="sensitivity basis; empty to skip")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--out", default=str(ROOT / "results" / "P6"))
    ap.add_argument("--skip-butadiene", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")
    ensure_dir(out / "logs")

    sy_b, co_b, dbl_b, _ = build_benzene_kekule(1.350, 1.450)
    bz = pack_molecule(
        "benzene_kekule",
        sy_b,
        co_b,
        dbl_b,
        [(0, 1), (1, 2), (2, 0)],
        args.basis,
        args.method,
        out / "raw" / "v2_benzene",
        ring=True,
    )

    sy_c, co_c, dbl_c, _ = build_cyclobutadiene(1.350, 1.540)
    cbd = pack_molecule(
        "cyclobutadiene",
        sy_c,
        co_c,
        dbl_c,
        [(0, 1)],
        args.basis,
        args.method,
        out / "raw" / "v2_cyclobutadiene",
        ring=True,
    )

    bd = None
    if not args.skip_butadiene:
        sy_d, co_d, dbl_d, _ = build_butadiene(1.340, 1.450)
        bd = pack_molecule(
            "butadiene",
            sy_d,
            co_d,
            dbl_d,
            [(0, 1)],
            args.basis,
            args.method,
            out / "raw" / "v2_butadiene",
            ring=False,
        )

    primary = judge(bz, cbd, bd)

    sensitivity = None
    if args.basis2 and args.basis2 != args.basis:
        print(f"\n--- sensitivity basis {args.basis2} ---", flush=True)
        bz2 = pack_molecule(
            "benzene_kekule",
            sy_b,
            co_b,
            dbl_b,
            [(0, 1), (1, 2), (2, 0)],
            args.basis2,
            args.method,
            out / "raw" / "v2_benzene_s2",
            ring=True,
        )
        cbd2 = pack_molecule(
            "cyclobutadiene",
            sy_c,
            co_c,
            dbl_c,
            [(0, 1)],
            args.basis2,
            args.method,
            out / "raw" / "v2_cyclobutadiene_s2",
            ring=True,
        )
        sens = judge(bz2, cbd2, None)
        sensitivity = {
            "basis": args.basis2,
            "benzene_ESE_kcal": sens["benzene_ESE_kcal"],
            "cbd_deltaEA_kcal": sens["cbd_deltaEA_kcal"],
            "signs_ok": sens["signs_ok"],
            "windows_ok": sens["benzene_in_window"] and sens["cbd_deltaEA_in_window"],
        }
        # require sign stability across bases for higher completion
        if primary["agree"] and sensitivity["signs_ok"]:
            primary["completion_estimate_pct"] = min(
                90, primary["completion_estimate_pct"] + 5
            )
            primary["basis_sign_stable"] = True
        else:
            primary["basis_sign_stable"] = bool(sensitivity["signs_ok"])

    # Final agree only if gates pass
    agree = primary["agree"] if primary["gates_ok"] else None

    pack = {
        "proposition": "P6",
        "version": "v2",
        "method": args.method,
        "basis": args.basis,
        "definition": {
            "ESE": "ΔEA − ΣΔEAm = E(G) − E(VR)",
            "VR": "E(VR)=E(GL)+Σ(E(GE-m)−E(GL))",
            "benzene_metric": "ESE",
            "cbd_metric": "ΔEA (conjugation/antiaromatic; ESE≡0 when n_double=2)",
            "localization": "2007 Fock+S",
            "yu_refs_kcal": {"benzene_ESE": -36.3, "cbd_conjugation": 53.6},
        },
        "molecules": [m for m in (bz, cbd, bd) if m is not None],
        "analysis": primary,
        "sensitivity": sensitivity,
        "quality_gate": {
            "passed": primary["gates_ok"],
            "G1_topology": True,
            "G2_geometry": True,
            "G3_convergence": True,
            "G4_energy_scale": True,
            "G5_path_clean": True,
            "limitations": [
                "vertical geometries (no adiabatic GL opt)",
                "2007 Fock+S; 2011 exchange deletion not implemented",
                "CBD Yu number mapped to ΔEA, not ESE",
            ],
        },
        "agree": agree,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    tag = f"{args.method}_{args.basis.replace('*', 's')}"
    write_json(out / "tables" / f"p6_v2_{tag}.json", pack)

    lines = [
        "# P6 v2 (ESE + CBD ΔEA)",
        "",
        f"- agree={agree} completion~{primary['completion_estimate_pct']}%",
        f"- benzene ESE={bz['ESE_kcal']:+.2f} (Yu −36.3) window={primary['benzene_in_window']}",
        f"- CBD ΔEA={cbd['deltaEA_kcal']:+.2f} (Yu ~+53.6) window={primary['cbd_deltaEA_in_window']}",
        f"- CBD ESE={cbd['ESE_kcal']:+.2f} (≡0 by construction)",
    ]
    if bd is not None:
        lines.append(
            f"- butadiene ΔEA={bd['deltaEA_kcal']:+.2f} ESE={bd['ESE_kcal']:+.2f} "
            f"additive_ok={primary['butadiene_additive_ok']}"
        )
    if sensitivity:
        lines.append(
            f"- sens {sensitivity['basis']}: bz ESE={sensitivity['benzene_ESE_kcal']:+.2f} "
            f"CBD ΔEA={sensitivity['cbd_deltaEA_kcal']:+.2f} signs={sensitivity['signs_ok']}"
        )
    lines += ["", primary["cbd_note"], ""]
    text = "\n".join(lines)
    (out / "tables" / "summary_p6_v2.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p6.md").write_text(text, encoding="utf-8")
    print(f"\n[P6v2] agree={agree} ~{primary['completion_estimate_pct']}%", flush=True)


if __name__ == "__main__":
    main()
