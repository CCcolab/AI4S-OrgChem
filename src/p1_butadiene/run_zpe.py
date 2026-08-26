"""
P1 — harmonic ZPE / thermal correction sensitivity on stored geometries.

Uses existing PySCF hessian + thermo (no new packages). Default: B3LYP/6-31G*
electronic geometries already in results/P1/raw/B3LYP_6-31gs/.

CE definitions unchanged; energies become E_el + ZPE (and optionally E_el+H_corr).
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from pyscf import dft, gto, hessian
from pyscf.hessian import thermo

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.common.units import HARTREE_TO_KCAL, ensure_dir, ha_to_kcal, write_json  # noqa: E402
from src.common.topology_p1 import check_p1_molecule  # noqa: E402

SPECIES = ["h2", "butadiene", "1_butene", "trans_2_butene", "n_butane"]


def load_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    n = int(lines[0].split()[0])
    syms, rows = [], []
    for ln in lines[2 : 2 + n]:
        p = ln.split()
        syms.append(p[0])
        rows.append([float(x) for x in p[1:4]])
    return syms, np.asarray(rows, dtype=float)


def make_mf(syms, coords, basis: str, method: str):
    mol = gto.Mole()
    mol.atom = [(syms[i], tuple(coords[i])) for i in range(len(syms))]
    mol.basis = basis
    mol.unit = "Angstrom"
    mol.verbose = 0
    mol.build()
    if method.upper() == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
    else:
        from pyscf import scf

        mf = scf.RHF(mol)
    mf.kernel()
    return mol, mf


def _scalar(x) -> float:
    if isinstance(x, (tuple, list)) and len(x) >= 1:
        x = x[0]
    return float(np.asarray(x, dtype=float).ravel()[0])


def zpe_of(mf, method: str) -> dict:
    if method.upper() == "B3LYP":
        hobj = hessian.rks.Hessian(mf)
    else:
        hobj = hessian.rhf.Hessian(mf)
    h = hobj.kernel()
    freq = thermo.harmonic_analysis(mf.mol, h)
    th = thermo.thermo(mf, freq["freq_au"], 298.15, 101325)
    wn = np.asarray(freq["freq_wavenumber"], dtype=float)
    return {
        "E_el_ha": float(mf.e_tot),
        "ZPE_ha": _scalar(th["ZPE"]),
        "E0K_ha": _scalar(th["E_0K"]),
        "H_298_ha": _scalar(th["H_tot"]),
        "n_imag": int(np.sum(wn < -10)),
        "scf_converged": bool(mf.converged),
    }


def ce_from(energies: dict[str, float]) -> tuple[float, float]:
    dh_bd = energies["n_butane"] - energies["butadiene"] - 2 * energies["h2"]
    dh_1b = energies["n_butane"] - energies["1_butene"] - energies["h2"]
    dh_t2 = energies["n_butane"] - energies["trans_2_butene"] - energies["h2"]
    return dh_bd - 2 * dh_1b, dh_bd - 2 * dh_t2


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP", choices=["B3LYP", "RHF"])
    ap.add_argument("--raw", default="")
    ap.add_argument("--out", default=str(ROOT / "results" / "P1"))
    args = ap.parse_args()

    out = Path(args.out)
    tag = f"{args.method}_{args.basis.replace('*', 's')}"
    raw = Path(args.raw) if args.raw else out / "raw" / tag
    ensure_dir(out / "tables")

    per: dict[str, dict] = {}
    for name in SPECIES:
        path = raw / f"{name}.xyz"
        print(f"[P1-ZPE] {name}", flush=True)
        syms, coords = load_xyz(path)
        defects = check_p1_molecule(
            name,
            "\n".join(
                f"{s} {coords[i,0]:.8f} {coords[i,1]:.8f} {coords[i,2]:.8f}"
                for i, s in enumerate(syms)
            ),
        )
        if defects:
            raise SystemExit(f"topology fail {name}: {defects}")
        mol, mf = make_mf(syms, coords, args.basis, args.method)
        info = zpe_of(mf, args.method)
        info["defects"] = defects
        per[name] = info
        print(
            f"  E={info['E_el_ha']:.8f} ZPE={info['ZPE_ha']*HARTREE_TO_KCAL:.3f} kcal "
            f"n_imag={info['n_imag']}",
            flush=True,
        )

    e_el = {k: v["E_el_ha"] for k, v in per.items()}
    e_zpe = {k: v["E0K_ha"] for k, v in per.items()}
    ce1_el, ce2_el = ce_from(e_el)
    ce1_z, ce2_z = ce_from(e_zpe)

    pack = {
        "proposition": "P1",
        "kind": "zpe_sensitivity",
        "method": args.method,
        "basis": args.basis,
        "protocol": "harmonic ZPE on stored electronic geometries (no re-opt)",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "species": per,
        "CE_electronic_kcal": {
            "CE1": ha_to_kcal(ce1_el),
            "CE2": ha_to_kcal(ce2_el),
            "sign_flip": bool(ce1_el < 0 and ce2_el > 0),
        },
        "CE_E0K_kcal": {
            "CE1": ha_to_kcal(ce1_z),
            "CE2": ha_to_kcal(ce2_z),
            "sign_flip": bool(ce1_z < 0 and ce2_z > 0),
        },
        "conclusion_zh": (
            "ZPE 校正后仍无 CE1<0 且 CE2>0 翻转"
            if not (ce1_z < 0 and ce2_z > 0)
            else "ZPE 校正后出现符号翻转——须复核 VERDICT"
        ),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(out / "tables" / f"ce_zpe_{tag}.json", pack)
    lines = [
        f"# P1 ZPE sensitivity ({tag})",
        "",
        f"- Electronic CE1/CE2 = {ha_to_kcal(ce1_el):+.3f} / {ha_to_kcal(ce2_el):+.3f}",
        f"- E+ZPE   CE1/CE2 = {ha_to_kcal(ce1_z):+.3f} / {ha_to_kcal(ce2_z):+.3f}",
        f"- Sign flip (E+ZPE): **{pack['CE_E0K_kcal']['sign_flip']}**",
        f"- {pack['conclusion_zh']}",
        "",
    ]
    (out / "summary_zpe.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
