"""
P6 pilot — parameter-free ESE = ΔEA − ΣΔEAm = E(G) − E(VR)
with VR: E(VR) = E(GL) + Σ_m (E(GE-m) − E(GL)).

2007 Fock+S localization (2011 exchange deletion deferred).
Benzene Kekulé vertical + rectangular cyclobutadiene vertical.
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
from src.localization.molecules import build_benzene_kekule, to_xyz  # noqa: E402

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


def check_ring(symbols, coords, n_c: int, doubles, r_max: float = 1.70) -> list[str]:
    problems = []
    dmat = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(dmat, np.inf)
    if float(dmat.min()) < 0.85:
        problems.append(f"dmin={float(dmat.min()):.3f}")
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


def build_cyclobutadiene(r_d: float = 1.350, r_s: float = 1.540, r_ch: float = R_CH):
    """Rectangular CBD: doubles (0,1)(2,3); singles (1,2)(3,0)."""
    coords = np.zeros((8, 3))
    half_d = r_d / 2.0
    half_s = r_s / 2.0
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

    symbols = ["C"] * 4 + ["H"] * 4
    doubles = [(0, 1), (2, 3)]
    singles = [(1, 2), (3, 0)]
    return symbols, coords, doubles, singles


def ese_pack(name, symbols, coords, doubles, ge_pairs, basis, method, out_raw: Path) -> dict:
    print(f"\n=== {name} ===", flush=True)
    n_c = sum(1 for s in symbols if s == "C")
    defects = check_ring(symbols, coords, n_c, doubles)
    if defects:
        raise SystemExit(f"{name} topology: {defects}")

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
    ese_from_vr = ha_to_kcal(e_g - e_vr)
    print(
        f"  ΔEA={dEA:+.3f} ΣΔEAm={sum_dEAm:+.3f} ESE={ese:+.3f} "
        f"(via VR {ese_from_vr:+.3f})",
        flush=True,
    )

    ensure_dir(out_raw)
    (out_raw / f"{name}.xyz").write_text(to_xyz(symbols, coords, name), encoding="utf-8")

    return {
        "molecule": name,
        "protocol": "vertical; 2007 Fock+S; ESE=ΔEA−ΣΔEAm",
        "G": {"E_ha": e_g},
        "GL": {"E_ha": e_gl},
        "GE_m": ge_list,
        "deltaEA_kcal": dEA,
        "sum_deltaEAm_kcal": sum_dEAm,
        "ESE_kcal": ese,
        "ESE_via_VR_kcal": ese_from_vr,
        "E_VR_ha": e_vr,
        "all_deltaEAm_positive": all(x["deltaEAm_kcal"] > 0 for x in ge_list),
        "quality_gate": {
            "topology_ok": True,
            "energy_scale_ok": abs(ese) < 120 and abs(dEA) < 80,
            "vr_identity_ok": abs(ese - ese_from_vr) < 1e-6,
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--out", default=str(ROOT / "results" / "P6"))
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")

    sy_b, co_b, dbl_b, _ = build_benzene_kekule(1.350, 1.450)
    bz = ese_pack(
        "benzene_kekule",
        sy_b,
        co_b,
        dbl_b,
        [(0, 1), (1, 2), (2, 0)],
        args.basis,
        args.method,
        out / "raw" / "benzene",
    )

    sy_c, co_c, dbl_c, _ = build_cyclobutadiene(1.350, 1.540)
    cbd = ese_pack(
        "cyclobutadiene",
        sy_c,
        co_c,
        dbl_c,
        [(0, 1)],
        args.basis,
        args.method,
        out / "raw" / "cyclobutadiene",
    )

    bz_ok = bz["ESE_kcal"] < 0 and 25 <= abs(bz["ESE_kcal"]) <= 45
    cbd_ok = cbd["ESE_kcal"] > 0 and 40 <= cbd["ESE_kcal"] <= 70
    gates_ok = (
        bz["quality_gate"]["topology_ok"]
        and cbd["quality_gate"]["topology_ok"]
        and bz["quality_gate"]["energy_scale_ok"]
        and cbd["quality_gate"]["energy_scale_ok"]
        and bz["all_deltaEAm_positive"]
        and cbd["all_deltaEAm_positive"]
    )

    if not gates_ok:
        agree = None
        completion = 25
    elif bz_ok and cbd_ok:
        agree = True
        completion = 70
    elif bz["ESE_kcal"] < 0 and cbd["ESE_kcal"] > 0:
        agree = True
        completion = 55
    else:
        agree = False
        completion = 40

    pack = {
        "proposition": "P6",
        "version": "pilot_v1",
        "method": args.method,
        "basis": args.basis,
        "definition": {
            "ESE": "ΔEA − ΣΔEAm",
            "VR": "E(VR)=E(GL)+Σ(E(GE-m)−E(GL))",
            "localization": "2007 Fock+S (exchange deletion not yet)",
            "yu_refs_kcal": {"benzene": -36.3, "cyclobutadiene": 53.0},
        },
        "molecules": [bz, cbd],
        "analysis": {
            "benzene_ESE_kcal": bz["ESE_kcal"],
            "cbd_ESE_kcal": cbd["ESE_kcal"],
            "benzene_in_window": bz_ok,
            "cbd_in_window": cbd_ok,
            "signs_ok": bz["ESE_kcal"] < 0 and cbd["ESE_kcal"] > 0,
            "agree": agree,
            "completion_estimate_pct": completion,
        },
        "quality_gate": {
            "passed": gates_ok,
            "limitations": [
                "vertical Kekulé / rectangle only",
                "2007 Fock+S; 2011 exchange deletion deferred",
                "CBD: single GE-m pair",
            ],
        },
        "agree": agree if gates_ok else None,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(out / "tables" / f"p6_pilot_{args.method}_{args.basis.replace('*', 's')}.json", pack)
    lines = [
        "# P6 pilot v1 (ESE)",
        "",
        f"- agree={pack['agree']} completion~{completion}%",
        f"- benzene ESE={bz['ESE_kcal']:+.2f} (Yu −36.3; ΔEA={bz['deltaEA_kcal']:+.2f} Σ={bz['sum_deltaEAm_kcal']:+.2f})",
        f"- CBD ESE={cbd['ESE_kcal']:+.2f} (Yu ~+53; ΔEA={cbd['deltaEA_kcal']:+.2f} Σ={cbd['sum_deltaEAm_kcal']:+.2f})",
        f"- signs_ok={pack['analysis']['signs_ok']} windows bz/cbd={bz_ok}/{cbd_ok}",
        "",
    ]
    text = "\n".join(lines)
    (out / "tables" / "summary_p6_pilot.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p6.md").write_text(text, encoding="utf-8")
    print(f"\n[P6pilot] agree={pack['agree']} ~{completion}%", flush=True)


if __name__ == "__main__":
    main()
