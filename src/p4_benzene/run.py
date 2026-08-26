"""
P4 — Benzene BLA scan with E / Ee / EN decomposition (B3LYP/6-31G*).

Public claim (Ch9 / preface): along equalization with dra = −drb, nuclear
repulsion is minimized when ra = rb (D6h); electronic energy prefers
alternation. Agreement: E_tot and EN minima near δ=0; path δ_max→0 has
ΔEN < 0 and ΔEe > 0 (qualitative).
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

from pyscf import dft, gto  # noqa: E402

from src.common.units import HARTREE_TO_KCAL, ensure_dir, ha_to_kcal, write_json  # noqa: E402
from src.p4_benzene.geometry import (  # noqa: E402
    bla_delta,
    build_benzene_d3h,
    build_benzene_equal,
    cc_bonds,
    check_topology,
    to_xyz,
)


def mol_from(symbols: list[str], coords: np.ndarray, basis: str) -> gto.Mole:
    atom = [(symbols[i], tuple(coords[i])) for i in range(len(symbols))]
    return gto.M(atom=atom, basis=basis, unit="Angstrom", verbose=0)


def scf_b3lyp(mol: gto.Mole) -> tuple[float, float, object]:
    mf = dft.RKS(mol)
    mf.xc = "B3LYP"
    mf.conv_tol = 1e-9
    e = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError("SCF not converged")
    en = float(mol.energy_nuc())
    return e, en, mf


def optimize_d6h(basis: str, r0: float = 1.397) -> tuple[float, list[str], np.ndarray, dict]:
    """BFGS on all coords starting from equal-bond hexagon; enforce planarity soft."""

    symbols, c0 = build_benzene_equal(r0)
    x0 = c0.ravel().copy()

    def energy(x):
        coords = x.reshape(-1, 3)
        coords = coords.copy()
        coords[:, 2] = 0.0
        defects = check_topology(symbols, coords)
        if defects:
            return 1e6
        mol = mol_from(symbols, coords, basis)
        e, _, _ = scf_b3lyp(mol)
        return e

    def grad(x):
        coords = x.reshape(-1, 3).copy()
        coords[:, 2] = 0.0
        mol = mol_from(symbols, coords, basis)
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
        mf.kernel()
        g = mf.nuc_grad_method().kernel()
        g[:, 2] = 0.0
        return g.ravel()

    res = minimize(
        energy,
        x0,
        method="BFGS",
        jac=grad,
        options={"gtol": 1e-4, "maxiter": 80, "disp": False},
    )
    coords = res.x.reshape(-1, 3)
    coords[:, 2] = 0.0
    bonds = cc_bonds(coords)
    r_mean = float(np.mean(bonds))
    meta = {
        "opt_success": bool(res.success),
        "nit": int(res.nit),
        "fun_ha": float(res.fun) if res.success or np.isfinite(res.fun) else None,
        "r_cc_mean": r_mean,
        "r_cc_bonds": bonds,
        "bla": bla_delta(coords),
        "defects": check_topology(symbols, coords),
    }
    if meta["defects"]:
        raise SystemExit(f"D6h opt topology failed: {meta['defects']}")
    return r_mean, symbols, coords, meta


def scan_bla(
    r0: float,
    basis: str,
    deltas: np.ndarray,
) -> list[dict]:
    rows = []
    for d in deltas:
        # dra = −drb: r_a = r0 + δ/2, r_b = r0 − δ/2  (δ = r_a − r_b ≥ 0)
        r_a = r0 + 0.5 * float(d)
        r_b = r0 - 0.5 * float(d)
        if r_b < 1.20:
            raise SystemExit(f"r_b too short at δ={d}")
        symbols, coords = build_benzene_d3h(r_a, r_b)
        defects = check_topology(symbols, coords)
        if defects:
            raise SystemExit(f"topology at δ={d}: {defects}")
        mol = mol_from(symbols, coords, basis)
        e, en, mf = scf_b3lyp(mol)
        ee = e - en
        bonds = cc_bonds(coords)
        rows.append(
            {
                "delta_ang": float(d),
                "r_a": r_a,
                "r_b": r_b,
                "r_mean": 0.5 * (r_a + r_b),
                "bonds": bonds,
                "E_ha": e,
                "EN_ha": en,
                "Ee_ha": ee,
                "scf_converged": bool(mf.converged),
                "defects": defects,
                "dmin": float(
                    np.min(
                        [
                            np.linalg.norm(coords[i] - coords[j])
                            for i in range(len(coords))
                            for j in range(i + 1, len(coords))
                        ]
                    )
                ),
            }
        )
        print(
            f"  δ={d:.3f}  E={e:.8f}  EN={en:.8f}  Ee={ee:.8f}  "
            f"r=({r_a:.4f}/{r_b:.4f})",
            flush=True,
        )
    return rows


def analyze(rows: list[dict]) -> dict:
    e0 = rows[0]["E_ha"]
    en0 = rows[0]["EN_ha"]
    ee0 = rows[0]["Ee_ha"]
    for r in rows:
        r["dE_kcal"] = ha_to_kcal(r["E_ha"] - e0)
        r["dEN_kcal"] = ha_to_kcal(r["EN_ha"] - en0)
        r["dEe_kcal"] = ha_to_kcal(r["Ee_ha"] - ee0)

    i_e = int(np.argmin([r["E_ha"] for r in rows]))
    i_en = int(np.argmin([r["EN_ha"] for r in rows]))
    i_ee = int(np.argmin([r["Ee_ha"] for r in rows]))

    # path from max δ → δ=0 (equalization)
    r_max = rows[-1]
    r_min = rows[0]
    dE = ha_to_kcal(r_min["E_ha"] - r_max["E_ha"])
    dEN = ha_to_kcal(r_min["EN_ha"] - r_max["EN_ha"])
    dEe = ha_to_kcal(r_min["Ee_ha"] - r_max["Ee_ha"])

    e_min_at_zero = i_e == 0
    en_min_at_zero = i_en == 0
    # book qualitative: equalization EN favorable (<0), Ee unfavorable (>0)
    en_favors_equal = dEN < 0
    ee_opposes_equal = dEe > 0
    # |ΔEN| > |ΔEe| along this path (optional stronger)
    en_dominates = abs(dEN) > abs(dEe) and en_favors_equal

    agree = bool(
        e_min_at_zero and en_min_at_zero and en_favors_equal and ee_opposes_equal
    )

    span = max(r["dE_kcal"] for r in rows) - min(r["dE_kcal"] for r in rows)
    return {
        "idx_E_min": i_e,
        "idx_EN_min": i_en,
        "idx_Ee_min": i_ee,
        "delta_at_E_min": rows[i_e]["delta_ang"],
        "delta_at_EN_min": rows[i_en]["delta_ang"],
        "delta_at_Ee_min": rows[i_ee]["delta_ang"],
        "path_max_delta_to_zero": {
            "dE_kcal": dE,
            "dEN_kcal": dEN,
            "dEe_kcal": dEe,
            "EN_favors_equalization": en_favors_equal,
            "Ee_opposes_equalization": ee_opposes_equal,
            "EN_dominates_magnitude": en_dominates,
        },
        "e_min_at_delta0": e_min_at_zero,
        "en_min_at_delta0": en_min_at_zero,
        "deltaE_span_kcal": span,
        "agree_criteria": (
            "E_tot & EN minima at δ=0; path δ_max→0: ΔEN<0 and ΔEe>0"
        ),
        "agree": agree,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--delta-max", type=float, default=0.12)
    ap.add_argument("--n-pts", type=int, default=13)
    ap.add_argument("--skip-opt", action="store_true", help="use r0=1.397 without opt")
    ap.add_argument("--out", default=str(ROOT / "results" / "P4"))
    args = ap.parse_args()

    out = Path(args.out)
    raw = ensure_dir(out / "raw" / f"{args.method}_{args.basis.replace('*', 's')}")
    print("[P4] optimize D6h benzene", flush=True)
    if args.skip_opt:
        r0 = 1.397
        symbols, coords = build_benzene_equal(r0)
        opt_meta = {"opt_success": None, "skipped": True, "r_cc_mean": r0}
    else:
        r0, symbols, coords, opt_meta = optimize_d6h(args.basis)
    print(f"  r0={r0:.5f} Å  bla={bla_delta(coords):.5f}  opt={opt_meta.get('opt_success')}", flush=True)
    (raw / "benzene_d6h.xyz").write_text(to_xyz(symbols, coords), encoding="utf-8")

    e0, en0, _ = scf_b3lyp(mol_from(symbols, coords, args.basis))
    print(f"  E(D6h)={e0:.8f}  EN={en0:.8f}", flush=True)

    deltas = np.linspace(0.0, args.delta_max, args.n_pts)
    print("[P4] BLA scan (dra=−drb)", flush=True)
    # use mean CC from opt as r0 for scan so δ=0 matches opt length scale
    rows = scan_bla(r0, args.basis, deltas)

    analysis = analyze(rows)
    all_scf = all(r["scf_converged"] for r in rows)
    all_topo = all(not r["defects"] for r in rows)
    dmin_ok = all(r["dmin"] >= 0.85 for r in rows)
    span_ok = analysis["deltaE_span_kcal"] < 40.0
    gate = {
        "topology_ok": all_topo,
        "dmin_ok": dmin_ok,
        "scf_ok": all_scf,
        "opt_ok": bool(opt_meta.get("opt_success", True)),
        "energy_scale_ok": span_ok,
        "deltaE_span_kcal": analysis["deltaE_span_kcal"],
        "defects_any": [r["defects"] for r in rows if r["defects"]],
        "passed": bool(all_topo and dmin_ok and all_scf and span_ok),
    }
    if not gate["passed"]:
        analysis["agree"] = None

    pack = {
        "proposition": "P4",
        "method": args.method,
        "basis": args.basis,
        "protocol": (
            "D6h BFGS opt; BLA scan with r_a=r0+δ/2, r_b=r0−δ/2 (dra=−drb); "
            "B3LYP single points; Ee=E−EN"
        ),
        "r0_ang": r0,
        "opt_meta": opt_meta,
        "E_d6h_ha": e0,
        "EN_d6h_ha": en0,
        "scan": rows,
        "analysis": analysis,
        "quality_gate": gate,
        "agree": analysis["agree"] if gate["passed"] else None,
        "book_refs": {
            "claim": "EN minimized at ra=rb; Ee destabilizing along equalization; |ΔEN|>|ΔEe| in GL→G",
            "chapter": "Ch9.3 / preface",
        },
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    tag = f"{args.method}_{args.basis.replace('*', 's')}"
    write_json(out / "tables" / f"bla_scan_{tag}.json", pack)

    # summary md
    lines = [
        "# P4 BLA scan summary",
        "",
        f"- r0 = {r0:.5f} Å · {args.method}/{args.basis}",
        f"- quality_gate.passed = {gate['passed']}",
        f"- agree = {pack['agree']}",
        f"- E min at δ = {analysis['delta_at_E_min']:.3f} Å",
        f"- EN min at δ = {analysis['delta_at_EN_min']:.3f} Å",
        f"- Ee min at δ = {analysis['delta_at_Ee_min']:.3f} Å",
        f"- path δ_max→0: ΔE={analysis['path_max_delta_to_zero']['dE_kcal']:+.3f}, "
        f"ΔEN={analysis['path_max_delta_to_zero']['dEN_kcal']:+.3f}, "
        f"ΔEe={analysis['path_max_delta_to_zero']['dEe_kcal']:+.3f} kcal/mol",
        "",
        "| δ | ΔE | ΔEN | ΔEe |",
        "|---|-----|------|------|",
    ]
    for r in rows:
        lines.append(
            f"| {r['delta_ang']:.3f} | {r['dE_kcal']:+.3f} | {r['dEN_kcal']:+.3f} | {r['dEe_kcal']:+.3f} |"
        )
    (out / "tables" / "summary_bla.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n[analysis]", flush=True)
    print(f"  gate.passed={gate['passed']}  agree={pack['agree']}", flush=True)
    print(
        f"  E/EN/Ee min δ = {analysis['delta_at_E_min']:.3f} / "
        f"{analysis['delta_at_EN_min']:.3f} / {analysis['delta_at_Ee_min']:.3f}",
        flush=True,
    )
    p = analysis["path_max_delta_to_zero"]
    print(
        f"  δ_max→0: ΔE={p['dE_kcal']:+.3f} ΔEN={p['dEN_kcal']:+.3f} ΔEe={p['dEe_kcal']:+.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
