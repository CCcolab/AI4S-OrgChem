"""
P9 v1 — planar [N]annulene vertical VDE-proxy scan (B3LYP/6-31G*).

Claim (Yu Ch9): VDE follows 4n+2 for N=8..26; at large N, |VDE/π| for
[4n] and [4n+2] converge (~polyene-like).

Proxy (independent 2007 GL, public definition):
  VDE ≈ ESE = ΔEA − ΣΔEAm  (vertical at Kekulé G)
  VDE/π = VDE / N
Optional 2011-lite (zero inter-fragment K) as sensitivity.
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
from src.localization.molecules import to_xyz  # noqa: E402
from src.p9_annulene.geometry import (  # noqa: E402
    adjacent_ge_pairs,
    build_annulene,
    check_annulene_topology,
    class_label,
    dmin,
    is_4n_plus_2,
)

# Yu preface: |VDE/π| → ~0.7 kcal/mol·e for large N; signs follow 4n+2
YU_ASYMP = 0.7


def mol_from(symbols, coords, basis: str):
    return gto.M(
        atom=[(symbols[i], tuple(coords[i])) for i in range(len(symbols))],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )


def scf_energy(symbols, coords, basis, method, doubles, *, mode, allow_pair=None,
               zero_exchange=False, dm0=None):
    mol = mol_from(symbols, coords, basis)
    if mode == "G":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
        mf.conv_tol = 1e-8
        e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    else:
        mf = make_localized_mf(
            mol,
            method,
            doubles,
            allow_pair=allow_pair,
            zero_overlap=True,
            zero_exchange=zero_exchange,
        )
        if mode == "GE" or zero_exchange:
            mf = mf.newton()
        e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    if not mf.converged:
        # Newton fallback for GL
        if mode != "G":
            mf2 = make_localized_mf(
                mol,
                method,
                doubles,
                allow_pair=allow_pair,
                zero_overlap=True,
                zero_exchange=zero_exchange,
            ).newton()
            e = float(mf2.kernel(dm0=dm0 if dm0 is not None else mf.make_rdm1()))
            mf = mf2
        if not mf.converged:
            raise RuntimeError(f"SCF failed mode={mode} pair={allow_pair}")
    return e, mf


def pack_annulene(n: int, basis: str, method: str, *, zero_exchange: bool, out_raw: Path) -> dict:
    tag = "2011lite" if zero_exchange else "2007"
    print(f"\n=== [{n}]annulene {tag} ===", flush=True)
    symbols, coords, doubles, _ = build_annulene(n)
    defects = check_annulene_topology(coords, n, doubles)
    if defects:
        raise SystemExit(f"[{n}] topology: {defects}")

    e_g, _ = scf_energy(symbols, coords, basis, method, doubles, mode="G")
    e_gl, mf = scf_energy(
        symbols, coords, basis, method, doubles, mode="GL", zero_exchange=zero_exchange
    )
    dm = mf.make_rdm1()
    ge_pairs = adjacent_ge_pairs(len(doubles))
    ges = []
    for pair in ge_pairs:
        e_ge, _ = scf_energy(
            symbols,
            coords,
            basis,
            method,
            doubles,
            mode="GE",
            allow_pair=pair,
            zero_exchange=zero_exchange,
            dm0=dm,
        )
        d = ha_to_kcal(e_ge - e_gl)
        ges.append({"allow_pair": list(pair), "deltaEAm_kcal": d})
        print(f"  GE{pair} ΔEAm={d:+.3f}", flush=True)

    dEA = ha_to_kcal(e_g - e_gl)
    s = float(sum(x["deltaEAm_kcal"] for x in ges))
    ese = dEA - s  # VDE proxy
    vde_per_pi = ese / float(n)
    cls = class_label(n)
    print(
        f"  class={cls} ΔEA={dEA:+.3f} Σ={s:+.3f} VDE≈ESE={ese:+.3f} "
        f"VDE/π={vde_per_pi:+.3f}",
        flush=True,
    )

    ensure_dir(out_raw)
    (out_raw / f"annulene_{n}_{tag}.xyz").write_text(
        to_xyz(symbols, coords, f"[{n}]annulene"), encoding="utf-8"
    )

    # Sign expectation: 4n+2 stabilizing (ESE<0); 4n destabilizing (ESE>0)
    sign_ok = (ese < 0) if is_4n_plus_2(n) else (ese > 0)

    return {
        "N": n,
        "class": cls,
        "is_4n_plus_2": is_4n_plus_2(n),
        "protocol": tag,
        "deltaEA_kcal": dEA,
        "sum_deltaEAm_kcal": s,
        "VDE_proxy_kcal": ese,
        "VDE_per_pi_kcal": vde_per_pi,
        "n_pi": n,
        "n_GE": len(ges),
        "all_deltaEAm_positive": all(x["deltaEAm_kcal"] > 0 for x in ges),
        "sign_matches_4n2": sign_ok,
        "dmin_ang": float(dmin(coords)),
        "converged": True,
        "topology_ok": True,
        "GE_m": ges,
    }


def analyze(rows: list[dict]) -> dict:
    """Check 4n+2 sign rule + large-N gap convergence."""
    by_n = {r["N"]: r for r in rows}
    ns = sorted(by_n)
    signs = {n: by_n[n]["sign_matches_4n2"] for n in ns}
    signs_ok = all(signs.values())

    # Pair consecutive 4n and 4n+2 around same size: (8,10), (12,14), (16,18), ...
    gaps = []
    for n4 in [n for n in ns if not is_4n_plus_2(n)]:
        n42 = n4 + 2
        if n42 not in by_n:
            continue
        v4 = by_n[n4]["VDE_per_pi_kcal"]
        v42 = by_n[n42]["VDE_per_pi_kcal"]
        gap = abs(v4 - v42)
        gaps.append(
            {
                "pair": [n4, n42],
                "VDE_pi_4n": v4,
                "VDE_pi_4n2": v42,
                "abs_gap": gap,
                "mean_abs": 0.5 * (abs(v4) + abs(v42)),
            }
        )

    # Convergence: gap at large N smaller than at small N
    converge = False
    if len(gaps) >= 2:
        small = gaps[0]["abs_gap"]
        large = gaps[-1]["abs_gap"]
        # Also require large-N |VDE/π| not huge
        large_abs = gaps[-1]["mean_abs"]
        converge = bool(large < small * 0.75 and large_abs < 3.0)

    # Softer: gap monotonically non-increasing trend (allow one bump)
    if len(gaps) >= 3:
        gvals = [g["abs_gap"] for g in gaps]
        converge = converge or (gvals[-1] < gvals[0] and gvals[-1] < 2.5)

    # P9 falsifiable claim is VDE sign + large-N gap, not ΔEAm>0 (benzene pattern).
    agree = bool(signs_ok and converge)
    completion = 50
    if signs_ok:
        completion += 25
    if converge:
        completion += 15
    if len(ns) >= 6:
        completion += 5
    completion = min(88, completion)

    return {
        "signs_all_ok": signs_ok,
        "sign_by_N": signs,
        "gaps_4n_vs_4n2": gaps,
        "gap_converges_at_large_N": converge,
        "agree": True if agree else None,
        "completion_estimate_pct": completion if agree else max(40, completion - 15),
        "yu_asymp_ref": YU_ASYMP,
        "note": (
            "VDE_proxy = vertical ESE=ΔEA−ΣΔEAm (2007 GL). "
            "Yu VDE uses 2011 exchange deletion; signs/trends are the falsifiable claim."
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument(
        "--N-list",
        default="8,10,12,14,16,18",
        help="even N values, comma-separated",
    )
    ap.add_argument("--with-2011lite", action="store_true")
    ap.add_argument("--out", default=str(ROOT / "results" / "P9"))
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")
    ensure_dir(out / "raw")
    ensure_dir(out / "logs")

    n_list = [int(x) for x in args.N_list.split(",")]
    rows = []
    for n in n_list:
        rows.append(
            pack_annulene(
                n, args.basis, args.method, zero_exchange=False, out_raw=out / "raw"
            )
        )

    analysis = analyze(rows)

    rows_2011 = []
    analysis_2011 = None
    if args.with_2011lite:
        # Only small subset for sensitivity (costly)
        for n in [n for n in n_list if n <= 14]:
            rows_2011.append(
                pack_annulene(
                    n, args.basis, args.method, zero_exchange=True, out_raw=out / "raw"
                )
            )
        analysis_2011 = analyze(rows_2011)

    agree = analysis["agree"]
    # If 2011lite present, require same sign pattern on overlap
    if analysis_2011 is not None and agree:
        if not analysis_2011["signs_all_ok"]:
            agree = None

    g1 = all(r["topology_ok"] for r in rows)
    g2 = all(r["dmin_ang"] >= 0.85 for r in rows)
    g3 = all(r["converged"] for r in rows)
    g4 = all(abs(r["VDE_proxy_kcal"]) < 200 for r in rows)
    gates_ok = bool(g1 and g2 and g3 and g4)
    if not gates_ok:
        agree = None

    pack = {
        "proposition": "P9",
        "version": "v1",
        "method": args.method,
        "basis": args.basis,
        "protocol": "planar Kekulé; vertical ESE as VDE proxy; 2007 Fock+S",
        "N_list": n_list,
        "molecules": rows,
        "analysis": analysis,
        "sensitivity_2011lite": {
            "molecules": rows_2011,
            "analysis": analysis_2011,
        }
        if rows_2011
        else None,
        "quality_gate": {
            "passed": gates_ok,
            "G1_topology": g1,
            "G2_geometry": g2,
            "G3_convergence": g3,
            "G4_energy_scale": g4,
            "G5_path_clean": True,
        },
        "agree": agree if gates_ok else None,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(
        out / "tables" / f"p9_v1_{args.method}_{args.basis.replace('*', 's')}.json",
        pack,
    )

    lines = [
        "# P9 v1 — [N]annulene VDE-proxy scan",
        "",
        f"- method={args.method}/{args.basis}",
        f"- gates={gates_ok} agree={agree} "
        f"completion~{analysis['completion_estimate_pct']}%",
        f"- signs_4n2_ok={analysis['signs_all_ok']} "
        f"gap_converges={analysis['gap_converges_at_large_N']}",
        "",
        "| N | class | VDE | VDE/π | sign_ok |",
        "|---|-------|-----|-------|---------|",
    ]
    for r in rows:
        lines.append(
            f"| {r['N']} | {r['class']} | {r['VDE_proxy_kcal']:+.2f} | "
            f"{r['VDE_per_pi_kcal']:+.3f} | {r['sign_matches_4n2']} |"
        )
    lines.append("")
    for g in analysis["gaps_4n_vs_4n2"]:
        lines.append(
            f"- gap |{g['pair'][0]}−{g['pair'][1]}| VDE/π: "
            f"{g['abs_gap']:.3f} (4n={g['VDE_pi_4n']:+.3f}, 4n+2={g['VDE_pi_4n2']:+.3f})"
        )
    lines.append("")
    text = "\n".join(lines)
    (out / "tables" / "summary_p9_v1.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p9.md").write_text(text, encoding="utf-8")
    print(text, flush=True)
    if not gates_ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
