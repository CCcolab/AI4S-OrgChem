"""
P1 L3 GL — Ch6 public 2007 definition (Fock + overlap π–π deletion between
double bonds), plus Fock-only diagnostic. Geometry: vertical at G, then 1D
r23 scan under Newton GL-SCF.
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
from src.localization.gl_2007 import make_gl_mf  # noqa: E402
from src.p1_butadiene.run import GEOMS, load_xyz_body  # noqa: E402

DBONDS = [(0, 1), (2, 3)]


def load_planar_butadiene(path: Path | None, basis: str) -> gto.Mole:
    xyz = load_xyz_body(path) if path and path.is_file() else GEOMS["butadiene"]
    mol = gto.M(atom=xyz, basis=basis, unit="Angstrom", verbose=0)
    c = mol.atom_coords(unit="Angstrom").copy()
    c[:, 2] = 0.0
    return gto.M(
        atom=[(mol.atom_symbol(i), tuple(c[i])) for i in range(mol.natm)],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )


def scf_energy(
    mol: gto.Mole, method: str, localized: bool, *, zero_overlap: bool
) -> tuple[float, object]:
    if localized:
        mf = make_gl_mf(mol, method, DBONDS, zero_overlap=zero_overlap).newton()
    elif method.upper() == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
    else:
        from pyscf import scf

        mf = scf.RHF(mol)
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"SCF not converged (localized={localized}, S0={zero_overlap})")
    return float(e), mf


def set_r23(mol: gto.Mole, r23: float) -> gto.Mole:
    """Scale C1–C2 vector to length r23; shift right-hand fragment."""
    c = mol.atom_coords(unit="Angstrom").copy()
    v = c[2] - c[1]
    v = v / (np.linalg.norm(v) + 1e-16) * r23
    shift = v - (c[2] - c[1])
    for i in (2, 3, 7, 8, 9):
        if i < len(c):
            c[i] = c[i] + shift
    c[:, 2] = 0.0
    return gto.M(
        atom=[(mol.atom_symbol(i), tuple(c[i])) for i in range(mol.natm)],
        basis=mol.basis,
        unit="Angstrom",
        verbose=0,
    )


def run_protocol(
    mol_g: gto.Mole,
    method: str,
    *,
    zero_overlap: bool,
    r_min: float,
    r_max: float,
    n_pts: int,
) -> dict:
    label = "Fock+S" if zero_overlap else "Fock-only"
    print(f"\n=== Protocol {label} ===", flush=True)

    e_g, _ = scf_energy(mol_g, method, False, zero_overlap=False)
    r_g = float(
        np.linalg.norm(
            mol_g.atom_coords(unit="Angstrom")[1] - mol_g.atom_coords(unit="Angstrom")[2]
        )
    )
    print(f"  E(G)={e_g:.8f} r23={r_g:.4f}", flush=True)

    e_gl_v, mf_v = scf_energy(mol_g, method, True, zero_overlap=zero_overlap)
    dE_v = ha_to_kcal(e_g - e_gl_v)
    print(f"  E(GL@G)={e_gl_v:.8f}  vertical ΔE={dE_v:+.3f} kcal", flush=True)

    scan = []
    best = None
    best_mol = mol_g
    for r in np.linspace(r_min, r_max, n_pts):
        mol_r = set_r23(mol_g, float(r))
        e_r, _ = scf_energy(mol_r, method, True, zero_overlap=zero_overlap)
        row = {"r23": float(r), "E_ha": e_r}
        scan.append(row)
        print(f"  r23={r:.3f} E={e_r:.8f}", flush=True)
        if best is None or e_r < best["E_ha"]:
            best = row
            best_mol = mol_r

    e_gl = best["E_ha"]
    r_gl = best["r23"]
    dE = ha_to_kcal(e_g - e_gl)
    at_edge = abs(r_gl - r_min) < 1e-9 or abs(r_gl - r_max) < 1e-9
    print(
        f"  ΔE(scan)={dE:+.3f} kcal  r23(GL)={r_gl:.4f}  edge={at_edge}",
        flush=True,
    )
    return {
        "zero_overlap": zero_overlap,
        "label": label,
        "G": {"E_ha": e_g, "r23_ang": r_g},
        "GL_vertical_at_G": {"E_ha": e_gl_v, "deltaE_kcal": dE_v},
        "GL": {"E_ha": e_gl, "r23_ang": r_gl, "scan": scan, "min_at_scan_edge": at_edge},
        "deltaE_kcal": dE,
        "delta_r23_ang": r_g - r_gl,
        "sign_destabilizing": dE > 0,
        "near_book_1p4": abs(dE - 1.4) < 1.0,
        "l3_agree_with_book_plus1p4": bool(dE > 0 and abs(dE - 1.4) < 1.5),
        "pi_aos": getattr(mf_v, "_gl_pi_aos", None),
        "fragments": getattr(mf_v, "_gl_fragments", None),
        "best_mol_xyz": (
            best_mol,
            "\n".join(
                f"{best_mol.atom_symbol(i)} "
                f"{best_mol.atom_coords(unit='Angstrom')[i,0]:.8f} "
                f"{best_mol.atom_coords(unit='Angstrom')[i,1]:.8f} "
                f"{best_mol.atom_coords(unit='Angstrom')[i,2]:.8f}"
                for i in range(best_mol.natm)
            ),
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", default="B3LYP", choices=["B3LYP", "RHF"])
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--out", default=str(ROOT / "results" / "P1"))
    ap.add_argument("--r-min", type=float, default=1.40)
    ap.add_argument("--r-max", type=float, default=1.52)
    ap.add_argument("--n-pts", type=int, default=13)
    ap.add_argument(
        "--protocols",
        default="both",
        choices=["both", "fock_s", "fock_only"],
        help="fock_s = Ch6 literal; fock_only = diagnostic; both = default",
    )
    args = ap.parse_args()

    out = Path(args.out)
    tag = f"{args.method}_{args.basis.replace('*', 's')}"
    raw = ensure_dir(out / "raw" / f"GL2007_{tag}")
    start = out / "raw" / f"B3LYP_{args.basis.replace('*', 's')}" / "butadiene.xyz"
    mol_g = load_planar_butadiene(start if start.is_file() else None, args.basis)

    do_fs = args.protocols in ("both", "fock_s")
    do_fo = args.protocols in ("both", "fock_only")

    protocols = []
    if do_fs:
        protocols.append(
            run_protocol(
                mol_g,
                args.method,
                zero_overlap=True,
                r_min=args.r_min,
                r_max=args.r_max,
                n_pts=args.n_pts,
            )
        )
    if do_fo:
        protocols.append(
            run_protocol(
                mol_g,
                args.method,
                zero_overlap=False,
                r_min=args.r_min,
                r_max=args.r_max,
                n_pts=args.n_pts,
            )
        )

    # Primary = Ch6 literal (Fock+S) if present, else first
    primary = next((p for p in protocols if p["zero_overlap"]), protocols[0])

    for p in protocols:
        mol_xyz, xyz = p.pop("best_mol_xyz")
        suffix = "FS" if p["zero_overlap"] else "Fonly"
        (raw / f"butadiene_GL_{suffix}_scanmin.xyz").write_text(
            f"{mol_xyz.natm}\nGL {p['label']} r23={p['GL']['r23_ang']}\n{xyz}\n",
            encoding="utf-8",
        )

    pack = {
        "proposition": "P1",
        "kind": "L3_GL_2007_butadiene",
        "method": args.method,
        "basis": args.basis,
        "definition": (
            "Ch6 2007: identify π-AOs via hcore; zero inter-double-bond π–π "
            "Fock and (primary) overlap blocks; Newton SCF; planar; "
            "ΔE=E(G)-E(GL) from r23 scan minimum under GL-SCF"
        ),
        "implementation_notes": {
            "primary_protocol": "Fock+S (Ch6 literal)",
            "diagnostic": "Fock-only retained for comparison",
            "scf_solver": "newton",
            "geometry": f"r23 1D scan {args.r_min}–{args.r_max} Å ({args.n_pts} pts)",
            "vs_book_2014": "exchange-integral deletion (2011/2014) not implemented",
            "book_plus1p4_source": "Ch12 cites 2014 method for butadiene +1.4",
        },
        "book_refs": {"deltaE_kcal": 1.4, "r23_GL_ang": 1.451, "r23_G_ang": 1.457},
        "protocols": protocols,
        "primary": {k: v for k, v in primary.items() if k != "best_mol_xyz"},
        "deltaE_kcal": primary["deltaE_kcal"],
        "sign_destabilizing": primary["sign_destabilizing"],
        "near_book_1p4": primary["near_book_1p4"],
        "l3_agree_with_book_plus1p4": primary["l3_agree_with_book_plus1p4"],
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(out / "tables" / f"gl2007_butadiene_{tag}.json", pack)
    print(
        f"\n[primary {primary['label']}] ΔE={primary['deltaE_kcal']:+.3f} kcal "
        f"(book +1.4); agree={primary['l3_agree_with_book_plus1p4']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
