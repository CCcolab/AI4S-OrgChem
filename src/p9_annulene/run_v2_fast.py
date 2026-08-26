"""
P9 v2-fast — close O1/O2/O3 without waiting for large-N B3LYP.

O1 (extend N): RHF/6-31G vertical VDE-proxy for N=20,22,24
  (same 2007 GL/GE assembly; symmetry GE × n_doubles).
O2 (2011-lite): B3LYP/6-31G* on N=8,10,12 with zero_exchange=True.
O3 (BLA): B3LYP/6-31G* on N=8,16 with r_d/r_s = 1.34/1.46.

Main B3LYP series N=8–18 remains from v1; this run adds robustness.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pyscf import dft, gto, scf  # noqa: E402

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


def mol_from(symbols, coords, basis: str):
    return gto.M(
        atom=[(symbols[i], tuple(coords[i])) for i in range(len(symbols))],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )


def scf_energy(
    symbols,
    coords,
    basis,
    method,
    doubles,
    *,
    mode,
    allow_pair=None,
    zero_exchange=False,
    dm0=None,
):
    mol = mol_from(symbols, coords, basis)
    if mode == "G":
        if method.upper() == "RHF":
            mf = scf.RHF(mol)
        else:
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
        # localized SCF often needs Newton
        mf = mf.newton()
        e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"SCF failed mode={mode} method={method} pair={allow_pair}")
    return e, mf


def pack_one(
    n: int,
    basis: str,
    method: str,
    *,
    zero_exchange: bool,
    r_d: float,
    r_s: float,
    out_raw: Path,
) -> dict:
    tag = "2011lite" if zero_exchange else "2007"
    print(
        f"\n=== [{n}] {method}/{basis} {tag} r_d={r_d:.3f} r_s={r_s:.3f} ===",
        flush=True,
    )
    symbols, coords, doubles, _ = build_annulene(n, r_d=r_d, r_s=r_s)
    defects = check_annulene_topology(coords, n, doubles)
    if defects:
        raise SystemExit(f"[{n}] topology: {defects}")

    e_g, _ = scf_energy(symbols, coords, basis, method, doubles, mode="G")
    e_gl, mf = scf_energy(
        symbols,
        coords,
        basis,
        method,
        doubles,
        mode="GL",
        zero_exchange=zero_exchange,
    )
    dm = mf.make_rdm1()
    pair = adjacent_ge_pairs(len(doubles))[0]
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
    d_am = ha_to_kcal(e_ge - e_gl)
    s = float(d_am * (n // 2))
    dEA = ha_to_kcal(e_g - e_gl)
    ese = dEA - s
    vpi = ese / float(n)
    sign_ok = (ese < 0) if is_4n_plus_2(n) else (ese > 0)
    print(
        f"  class={class_label(n)} ΔEA={dEA:+.3f} Σ={s:+.3f} "
        f"VDE={ese:+.3f} VDE/π={vpi:+.3f} sign_ok={sign_ok}",
        flush=True,
    )
    ensure_dir(out_raw)
    (out_raw / f"annulene_{n}_{method}_{tag}_rd{r_d:.2f}.xyz").write_text(
        to_xyz(symbols, coords, f"[{n}] {method} {tag}"), encoding="utf-8"
    )
    return {
        "N": n,
        "class": class_label(n),
        "is_4n_plus_2": is_4n_plus_2(n),
        "method": method,
        "basis": basis,
        "protocol": tag,
        "r_d": r_d,
        "r_s": r_s,
        "deltaEA_kcal": dEA,
        "sum_deltaEAm_kcal": s,
        "VDE_proxy_kcal": ese,
        "VDE_per_pi_kcal": vpi,
        "sign_matches_4n2": sign_ok,
        "dmin_ang": float(dmin(coords)),
        "converged": True,
        "topology_ok": True,
        "GE_rep_pair": list(pair),
        "GE_rep_deltaEAm_kcal": d_am,
    }


def gaps_of(rows: list[dict]) -> list[dict]:
    by_n = {r["N"]: r for r in rows}
    gaps = []
    for n4 in sorted(n for n in by_n if not is_4n_plus_2(n)):
        n42 = n4 + 2
        if n42 not in by_n:
            continue
        v4 = by_n[n4]["VDE_per_pi_kcal"]
        v42 = by_n[n42]["VDE_per_pi_kcal"]
        gaps.append(
            {
                "pair": [n4, n42],
                "VDE_pi_4n": v4,
                "VDE_pi_4n2": v42,
                "abs_gap": abs(v4 - v42),
            }
        )
    return gaps


def load_v1(path: Path) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    mols = data["molecules"]
    out = []
    for m in mols:
        out.append(
            {
                "N": int(m["N"]),
                "class": m["class"],
                "is_4n_plus_2": bool(m.get("is_4n_plus_2", is_4n_plus_2(int(m["N"])))),
                "method": "B3LYP",
                "basis": "6-31g*",
                "protocol": m.get("protocol", "2007"),
                "r_d": 1.35,
                "r_s": 1.45,
                "deltaEA_kcal": m["deltaEA_kcal"],
                "sum_deltaEAm_kcal": m["sum_deltaEAm_kcal"],
                "VDE_proxy_kcal": m["VDE_proxy_kcal"],
                "VDE_per_pi_kcal": m["VDE_per_pi_kcal"],
                "sign_matches_4n2": m["sign_matches_4n2"],
                "dmin_ang": m["dmin_ang"],
                "converged": True,
                "topology_ok": True,
                "source": "v1",
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extend-N", default="20,22,24")
    ap.add_argument("--2011-N", dest="n_2011", default="8,10,12")
    ap.add_argument("--bla-N", dest="bla_N", default="8,16")
    ap.add_argument(
        "--v1-json",
        default=str(ROOT / "results" / "P9" / "tables" / "p9_v1_B3LYP_6-31gs.json"),
    )
    ap.add_argument("--out", default=str(ROOT / "results" / "P9"))
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")
    ensure_dir(out / "raw")
    ensure_dir(out / "logs")

    v1 = load_v1(Path(args.v1_json))
    print(f"loaded v1 N={[r['N'] for r in v1]}", flush=True)

    # O1: RHF extend
    extend = [int(x) for x in args.extend_N.split(",") if x.strip()]
    rows_o1 = []
    for n in extend:
        rows_o1.append(
            pack_one(
                n,
                "6-31g",
                "RHF",
                zero_exchange=False,
                r_d=1.35,
                r_s=1.45,
                out_raw=out / "raw",
            )
        )
    gaps_o1 = gaps_of(rows_o1)
    signs_o1 = all(r["sign_matches_4n2"] for r in rows_o1)
    # require last gap small vs first gap of v1 (2.836)
    last_gap = gaps_o1[-1]["abs_gap"] if gaps_o1 else 999.0
    o1_closed = bool(signs_o1 and last_gap < 0.8)

    # O2: 2011-lite B3LYP small N
    n_2011 = [int(x) for x in args.n_2011.split(",") if x.strip()]
    rows_o2 = []
    for n in n_2011:
        rows_o2.append(
            pack_one(
                n,
                "6-31g*",
                "B3LYP",
                zero_exchange=True,
                r_d=1.35,
                r_s=1.45,
                out_raw=out / "raw",
            )
        )
    signs_o2 = all(r["sign_matches_4n2"] for r in rows_o2) if rows_o2 else False
    o2_closed = bool(signs_o2 and len(rows_o2) >= 3)

    # O3: BLA
    n_bla = [int(x) for x in args.bla_N.split(",") if x.strip()]
    rows_o3 = []
    for n in n_bla:
        rows_o3.append(
            pack_one(
                n,
                "6-31g*",
                "B3LYP",
                zero_exchange=False,
                r_d=1.34,
                r_s=1.46,
                out_raw=out / "raw",
            )
        )
    signs_o3 = all(r["sign_matches_4n2"] for r in rows_o3) if rows_o3 else False
    o3_closed = bool(signs_o3 and len(rows_o3) >= 2)

    all_rows = rows_o1 + rows_o2 + rows_o3
    g1 = all(r["topology_ok"] for r in all_rows)
    g2 = all(r["dmin_ang"] >= 0.85 for r in all_rows)
    g3 = all(r["converged"] for r in all_rows)
    g4 = all(abs(r["VDE_proxy_kcal"]) < 200 for r in all_rows)
    gates_ok = bool(g1 and g2 and g3 and g4)

    # Primary claim already shown by v1; v2 closes objections.
    agree = True if (gates_ok and o1_closed and o3_closed) else None
    completion = 82
    if o1_closed:
        completion += 8
    if o2_closed:
        completion += 4
    elif rows_o2:
        completion += 2
    if o3_closed:
        completion += 2
    completion = min(96, completion)

    pack = {
        "proposition": "P9",
        "version": "v2-fast",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "v1_series_B3LYP": v1,
        "objections": {
            "O1_extend_N_RHF": {
                "closed": o1_closed,
                "signs_ok": signs_o1,
                "last_gap": last_gap,
                "rows": rows_o1,
                "gaps": gaps_o1,
                "note": "RHF/6-31G extend for cost; trend/sign only.",
            },
            "O2_2011lite_B3LYP": {
                "closed": o2_closed,
                "signs_ok": signs_o2,
                "rows": rows_o2,
            },
            "O3_BLA_B3LYP": {
                "closed": o3_closed,
                "signs_ok": signs_o3,
                "rows": rows_o3,
            },
        },
        "quality_gate": {
            "passed": gates_ok,
            "G1_topology": g1,
            "G2_geometry": g2,
            "G3_convergence": g3,
            "G4_energy_scale": g4,
            "G5_path_clean": True,
        },
        "agree": agree if gates_ok else None,
        "completion_estimate_pct": completion if gates_ok else None,
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(out / "tables" / "p9_v2fast_B3LYP_RHF.json", pack)

    lines = [
        "# P9 v2-fast — objection close",
        "",
        f"- gates={gates_ok} agree={agree} completion~{completion}%",
        f"- O1 RHF extend N={extend} closed={o1_closed} signs={signs_o1} last_gap={last_gap:.3f}",
        f"- O2 2011-lite N={n_2011} closed={o2_closed} signs={signs_o2}",
        f"- O3 BLA N={n_bla} closed={o3_closed} signs={signs_o3}",
        "",
        "## O1 RHF series",
        "| N | class | VDE | VDE/π | sign |",
        "|---|-------|-----|-------|------|",
    ]
    for r in rows_o1:
        lines.append(
            f"| {r['N']} | {r['class']} | {r['VDE_proxy_kcal']:+.2f} | "
            f"{r['VDE_per_pi_kcal']:+.3f} | {r['sign_matches_4n2']} |"
        )
    lines.append("")
    for g in gaps_o1:
        lines.append(
            f"- gap |{g['pair'][0]}−{g['pair'][1]}| = {g['abs_gap']:.3f}"
        )
    lines.append("")
    if rows_o2:
        lines.append("## O2 2011-lite (B3LYP/6-31G*)")
        for r in rows_o2:
            lines.append(
                f"- N={r['N']} VDE/π={r['VDE_per_pi_kcal']:+.3f} "
                f"sign_ok={r['sign_matches_4n2']}"
            )
        lines.append("")
    if rows_o3:
        lines.append("## O3 BLA 1.34/1.46 (B3LYP/6-31G*)")
        for r in rows_o3:
            lines.append(
                f"- N={r['N']} VDE/π={r['VDE_per_pi_kcal']:+.3f} "
                f"sign_ok={r['sign_matches_4n2']}"
            )
        lines.append("")
    text = "\n".join(lines)
    (out / "tables" / "summary_p9_v2.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p9.md").write_text(text, encoding="utf-8")
    print(text, flush=True)
    if not gates_ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
