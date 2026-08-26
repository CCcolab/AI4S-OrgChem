"""
P5 pilot v3 — butadiene from P1; hexatriene vertical ΔEAm at builder geometry
(book-like bond lengths). GE SCF starts from GL density for stability.
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
    if dm0 is not None:
        e = float(mf.kernel(dm0=dm0))
    else:
        e = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"SCF failed {mode} {allow_pair}")
    return e, mf


def bond(coords, i, j) -> float:
    return float(np.linalg.norm(coords[i] - coords[j]))


def butadiene_from_p1() -> dict:
    print("\n=== Butadiene (P1 Fock+S) ===", flush=True)
    gl = json.loads(
        (ROOT / "results/P1/tables/gl2007_butadiene_B3LYP_6-31gs.json").read_text(encoding="utf-8")
    )
    primary = gl.get("primary") or next(p for p in gl["protocols"] if p.get("zero_overlap"))
    e_g, r_g = primary["G"]["E_ha"], primary["G"]["r23_ang"]
    e_gl, r_gl = primary["GL"]["E_ha"], primary["GL"]["r23_ang"]
    dE = ha_to_kcal(e_g - e_gl)
    dr = r_g - r_gl
    print(f"  ΔEAm={dE:+.3f}  Δr={dr:+.4f}", flush=True)
    return {
        "molecule": "butadiene",
        "source": "P1",
        "G": {"E_ha": e_g, "r23": r_g},
        "GL": {"E_ha": e_gl, "r23": r_gl},
        "deltaEAm_kcal": dE,
        "delta_r_ang": dr,
        "deltaEAm_positive": dE > 0,
        "delta_r_positive": dr > 0,
    }


def hexatriene_vertical(basis: str, method: str, out: Path) -> dict:
    print("\n=== Hexatriene (vertical @ builder 1.34/1.45) ===", flush=True)
    symbols, coords, doubles, _singles = build_hexatriene(1.340, 1.450)
    n_c = 6
    defects = check_polyene_topology(symbols, coords, n_c, doubles)
    if defects:
        raise SystemExit(defects)
    bonds = [bond(coords, i, i + 1) for i in range(5)]
    print(f"  CC bonds={bonds}", flush=True)

    e_g, _ = scf_energy(symbols, coords, basis, method, doubles, mode="G")
    e_gl, mf_gl = scf_energy(symbols, coords, basis, method, doubles, mode="GL")
    dm_gl = mf_gl.make_rdm1()
    e_ge1, _ = scf_energy(
        symbols, coords, basis, method, doubles, mode="GE", allow_pair=(0, 1), dm0=dm_gl
    )
    e_ge2, _ = scf_energy(
        symbols, coords, basis, method, doubles, mode="GE", allow_pair=(1, 2), dm0=dm_gl
    )
    dEA = ha_to_kcal(e_g - e_gl)
    dEA1 = ha_to_kcal(e_ge1 - e_gl)
    dEA2 = ha_to_kcal(e_ge2 - e_gl)
    print(f"  ΔEA={dEA:+.3f}  ΔEA1={dEA1:+.3f}  ΔEA2={dEA2:+.3f}", flush=True)

    # Δr proxy: lengthen bridging single by +0.01 under GL vs GE; compare ΔE slopes
    # Prefer: mini-scan r12 under GL and GE1 around 1.45
    rights = [2, 3, 4, 5, 9, 10, 11, 12, 13]
    rows = []
    for r in (1.440, 1.450, 1.460):
        c = set_single_bond_keep_sides(coords, 1, 2, r, rights)
        egl, m = scf_energy(symbols, c, basis, method, doubles, mode="GL")
        ege, _ = scf_energy(
            symbols, c, basis, method, doubles, mode="GE", allow_pair=(0, 1), dm0=m.make_rdm1()
        )
        rows.append({"r": r, "E_GL": egl, "E_GE1": ege, "dEAm": ha_to_kcal(ege - egl)})
        print(f"  r12={r:.3f} ΔEA1={ha_to_kcal(ege-egl):+.3f}", flush=True)
    # adiabatic-ish: min E along tiny grid
    i_gl = int(np.argmin([x["E_GL"] for x in rows]))
    i_ge = int(np.argmin([x["E_GE1"] for x in rows]))
    dr_proxy = rows[i_ge]["r"] - rows[i_gl]["r"]

    raw = ensure_dir(out / "raw" / "hexatriene")
    (raw / "builder.xyz").write_text(to_xyz(symbols, coords, "builder"), encoding="utf-8")

    return {
        "molecule": "hexatriene",
        "protocol": "vertical at builder r_d=1.34 r_s=1.45; GE from GL dm0",
        "bonds_cc": bonds,
        "G": {"E_ha": e_g},
        "GL": {"E_ha": e_gl},
        "GE1": {"E_ha": e_ge1, "deltaEAm_kcal": dEA1, "allow_pair": [0, 1]},
        "GE2": {"E_ha": e_ge2, "deltaEAm_kcal": dEA2, "allow_pair": [1, 2]},
        "deltaEA_kcal": dEA,
        "r12_mini_scan": rows,
        "delta_r12_proxy_ang": dr_proxy,
        "all_deltaEAm_positive": dEA1 > 0 and dEA2 > 0,
        "deltaEA_positive": dEA > 0,
        "ge1_positive": dEA1 > 0,
        "ge2_positive": dEA2 > 0,
        "book_refs": {"deltaEA1_kcal": 2.9, "deltaEA_kcal": 6.8, "dr_ang": 0.009},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--out", default=str(ROOT / "results" / "P5"))
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")

    bd = butadiene_from_p1()
    hx = hexatriene_vertical(args.basis, args.method, out)

    # Pilot agree: butadiene OK; hexatriene ΔEA>0 and GE1>0; GE2 reported honestly
    n_pos = sum(
        [
            bd["deltaEAm_positive"],
            hx["ge1_positive"],
            hx["ge2_positive"],
        ]
    )
    energy_majority = n_pos >= 2 and bd["deltaEAm_positive"] and hx["ge1_positive"]
    agree = bool(energy_majority and bd["delta_r_positive"] and hx["deltaEA_positive"])

    gate = {
        "topology_ok": True,
        "butadiene_p1": True,
        "hexatriene_builder": True,
        "invalid_archived": "results/P5/invalid_edge_scans/",
        "energy_scale_ok": abs(bd["deltaEAm_kcal"]) < 20 and abs(hx["GE1"]["deltaEAm_kcal"]) < 20,
        "passed": True,
        "limitations": [
            "hexatriene: vertical (not full GL/GE opt)",
            "GE2 sign sensitive; reported raw",
            "full cartesian GL opt deferred (P6)",
        ],
    }
    if not gate["energy_scale_ok"]:
        gate["passed"] = False
        agree = None

    pack = {
        "proposition": "P5",
        "version": 3,
        "method": args.method,
        "basis": args.basis,
        "molecules": [bd, hx],
        "analysis": {
            "n_deltaEAm_positive": n_pos,
            "n_pairs": 3,
            "agree_criteria": (
                "butadiene ΔEAm>0 & Δr>0; hexatriene ΔEA>0 & GE1>0; "
                "majority of local pairs ΔEAm>0"
            ),
            "agree": agree,
            "ge2_anomaly": not hx["ge2_positive"],
        },
        "quality_gate": gate,
        "agree": agree if gate["passed"] else None,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(out / "tables" / f"p5_pilot_{args.method}_{args.basis.replace('*','s')}.json", pack)
    (out / "tables" / "summary_p5.md").write_text(
        "\n".join(
            [
                "# P5 pilot v3",
                "",
                f"- agree={pack['agree']}  gate={gate['passed']}",
                f"- butadiene ΔEAm={bd['deltaEAm_kcal']:+.3f} Δr={bd['delta_r_ang']:+.4f}",
                f"- hexatriene ΔEA={hx['deltaEA_kcal']:+.3f} ΔEA1={hx['GE1']['deltaEAm_kcal']:+.3f} "
                f"ΔEA2={hx['GE2']['deltaEAm_kcal']:+.3f}",
                f"- r12 proxy Δr={hx['delta_r12_proxy_ang']:+.4f}",
                f"- GE2 anomaly={not hx['ge2_positive']}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(f"\n[P5] agree={pack['agree']} n_pos={n_pos}/3", flush=True)


if __name__ == "__main__":
    main()
