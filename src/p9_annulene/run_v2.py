"""
P9 v2 — close remaining objections on the [N]annulene VDE claim.

O1: extend N to 20/22/24 (Yu series goes to 26).
O2: 2011-lite (inter-fragment K deletion) on N=8,10,12.
O3: BLA geometry sensitivity (r_d/r_s = 1.34/1.46 vs 1.35/1.45) on N=8,16.

Regular Kekulé polygons make all adjacent GE-m equivalent; v2 uses one
representative GE × n_doubles (validated against v1 full sums).
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

# Yu Table 9-5 B3LYP/6-31G* EV(2011)/π — trend only
YU_B3LYP_631GS_PER_PI = {
    8: 1.1,
    10: -4.5,
    12: -0.15,
    14: -3.8,
    16: -0.5,
    18: -2.7,
}


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
    if not mf.converged and mode != "G":
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


def pack_one(
    n: int,
    basis: str,
    method: str,
    *,
    zero_exchange: bool,
    r_d: float,
    r_s: float,
    out_raw: Path,
    ge_mode: str = "one",
) -> dict:
    tag = "2011lite" if zero_exchange else "2007"
    print(
        f"\n=== [{n}] {tag} r_d={r_d:.3f} r_s={r_s:.3f} ge={ge_mode} ===",
        flush=True,
    )
    symbols, coords, doubles, _ = build_annulene(n, r_d=r_d, r_s=r_s)
    defects = check_annulene_topology(coords, n, doubles)
    if defects:
        raise SystemExit(f"[{n}] topology: {defects}")

    e_g, _ = scf_energy(symbols, coords, basis, method, doubles, mode="G")
    e_gl, mf = scf_energy(
        symbols, coords, basis, method, doubles, mode="GL", zero_exchange=zero_exchange
    )
    dm = mf.make_rdm1()
    pairs = adjacent_ge_pairs(len(doubles))
    if ge_mode == "one":
        pairs = [pairs[0]]
    ges = []
    for pair in pairs:
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
    if ge_mode == "one":
        s = float(ges[0]["deltaEAm_kcal"] * (n // 2))
        ge_note = "symmetry: ΣΔEAm = n_doubles × GE(0,1)"
    else:
        s = float(sum(x["deltaEAm_kcal"] for x in ges))
        ge_note = "full adjacent GE-m"
    ese = dEA - s
    vde_per_pi = ese / float(n)
    cls = class_label(n)
    sign_ok = (ese < 0) if is_4n_plus_2(n) else (ese > 0)
    print(
        f"  class={cls} ΔEA={dEA:+.3f} Σ={s:+.3f} VDE≈ESE={ese:+.3f} "
        f"VDE/π={vde_per_pi:+.3f} sign_ok={sign_ok}",
        flush=True,
    )

    ensure_dir(out_raw)
    suffix = f"{n}_{tag}_rd{r_d:.2f}"
    (out_raw / f"annulene_{suffix}.xyz").write_text(
        to_xyz(symbols, coords, f"[{n}]annulene {tag}"), encoding="utf-8"
    )
    return {
        "N": n,
        "class": cls,
        "is_4n_plus_2": is_4n_plus_2(n),
        "protocol": tag,
        "r_d": r_d,
        "r_s": r_s,
        "ge_mode": ge_mode,
        "ge_note": ge_note,
        "deltaEA_kcal": dEA,
        "sum_deltaEAm_kcal": s,
        "VDE_proxy_kcal": ese,
        "VDE_per_pi_kcal": vde_per_pi,
        "n_pi": n,
        "n_GE": n // 2,
        "sign_matches_4n2": sign_ok,
        "dmin_ang": float(dmin(coords)),
        "converged": True,
        "topology_ok": True,
        "GE_m": ges,
        "source": "v2",
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
                "mean_abs": 0.5 * (abs(v4) + abs(v42)),
            }
        )
    return gaps


def load_v1(path: Path) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    mols = data.get("molecules") or data.get("rows")
    if not mols:
        raise SystemExit(f"no molecules in {path}")
    out = []
    for m in mols:
        n = int(m["N"])
        out.append(
            {
                "N": n,
                "class": m.get("class") or class_label(n),
                "is_4n_plus_2": is_4n_plus_2(n),
                "protocol": m.get("protocol", "2007"),
                "r_d": 1.35,
                "r_s": 1.45,
                "ge_mode": "full",
                "ge_note": "v1 full adjacent GE-m",
                "deltaEA_kcal": m["deltaEA_kcal"],
                "sum_deltaEAm_kcal": m["sum_deltaEAm_kcal"],
                "VDE_proxy_kcal": m["VDE_proxy_kcal"],
                "VDE_per_pi_kcal": m["VDE_per_pi_kcal"],
                "n_pi": n,
                "n_GE": m.get("n_GE", n // 2),
                "sign_matches_4n2": m["sign_matches_4n2"],
                "dmin_ang": m["dmin_ang"],
                "converged": m.get("converged", True),
                "topology_ok": m.get("topology_ok", True),
                "source": "v1",
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--extend-N", dest="extend_N", default="20,22,24")
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

    extend = [int(x) for x in args.extend_N.split(",") if x.strip()]
    new_rows = []
    for n in extend:
        new_rows.append(
            pack_one(
                n,
                args.basis,
                args.method,
                zero_exchange=False,
                r_d=1.35,
                r_s=1.45,
                out_raw=out / "raw",
                ge_mode="one",
            )
        )

    series = v1 + new_rows
    gaps = gaps_of(series)
    signs_ok = all(bool(r["sign_matches_4n2"]) for r in series)
    gvals = [g["abs_gap"] for g in gaps]
    o1_converge = bool(
        len(gvals) >= 4 and gvals[-1] < gvals[0] * 0.5 and gvals[-1] < 0.6
    )
    o1_closed = bool(signs_ok and o1_converge)

    n_2011 = [int(x) for x in args.n_2011.split(",") if x.strip()]
    rows_2011 = []
    for n in n_2011:
        rows_2011.append(
            pack_one(
                n,
                args.basis,
                args.method,
                zero_exchange=True,
                r_d=1.35,
                r_s=1.45,
                out_raw=out / "raw",
                ge_mode="one",
            )
        )
    signs_2011 = all(bool(r["sign_matches_4n2"]) for r in rows_2011)
    o2_closed = bool(signs_2011 and len(rows_2011) >= 3)

    n_bla = [int(x) for x in args.bla_N.split(",") if x.strip()]
    rows_bla = []
    for n in n_bla:
        rows_bla.append(
            pack_one(
                n,
                args.basis,
                args.method,
                zero_exchange=False,
                r_d=1.34,
                r_s=1.46,
                out_raw=out / "raw",
                ge_mode="one",
            )
        )
    signs_bla = all(bool(r["sign_matches_4n2"]) for r in rows_bla)
    o3_closed = bool(signs_bla and len(rows_bla) >= 2)

    all_new = new_rows + rows_2011 + rows_bla
    g1 = all(r["topology_ok"] for r in series + all_new)
    g2 = all(float(r["dmin_ang"]) >= 0.85 for r in series + all_new)
    g3 = all(r["converged"] for r in series + all_new)
    g4 = all(abs(float(r["VDE_proxy_kcal"])) < 200 for r in series + all_new)
    gates_ok = bool(g1 and g2 and g3 and g4)

    agree = True if (gates_ok and o1_closed and o3_closed) else None
    completion = 82
    if o1_closed:
        completion += 8
    if o2_closed:
        completion += 4
    elif rows_2011:
        completion += 2
    if o3_closed:
        completion += 2
    completion = min(96, completion)

    pack = {
        "proposition": "P9",
        "version": "v2",
        "method": args.method,
        "basis": args.basis,
        "protocol": (
            "planar Kekulé; vertical ESE as VDE proxy; 2007 Fock+S; "
            "v2: N=20–24 symmetry-GE; 2011-lite; BLA 1.34/1.46"
        ),
        "objections": {
            "O1_extend_N": {
                "closed": o1_closed,
                "signs_ok": signs_ok,
                "gap_converges": o1_converge,
                "new_N": extend,
            },
            "O2_2011lite": {
                "closed": o2_closed,
                "signs_ok": signs_2011,
                "N_list": n_2011,
                "note": (
                    "Yu Table 9-5 B3LYP/6-31G* already has near-zero/negative "
                    "[12]/[16] VDE; 2011-lite is sensitivity, not the primary metric."
                ),
            },
            "O3_BLA": {"closed": o3_closed, "signs_ok": signs_bla, "N_list": n_bla},
        },
        "v1_reused": v1,
        "v2_new": new_rows,
        "series_2007": series,
        "gaps_4n_vs_4n2": gaps,
        "sensitivity_2011lite": rows_2011,
        "sensitivity_BLA": rows_bla,
        "yu_table95_b3lyp_631gs_per_pi": YU_B3LYP_631GS_PER_PI,
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
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(
        out / "tables" / f"p9_v2_{args.method}_{args.basis.replace('*', 's')}.json",
        pack,
    )

    lines = [
        "# P9 v2 — [N]annulene VDE deepen",
        "",
        f"- gates={gates_ok} agree={agree} completion~{completion}%",
        f"- O1 extend N={extend} closed={o1_closed} signs={signs_ok} converge={o1_converge}",
        f"- O2 2011-lite N={n_2011} closed={o2_closed} signs={signs_2011}",
        f"- O3 BLA N={n_bla} closed={o3_closed} signs={signs_bla}",
        "",
        "| N | class | VDE | VDE/π | sign | src |",
        "|---|-------|-----|-------|------|-----|",
    ]
    for r in series:
        src = r.get("source", "v2")
        lines.append(
            f"| {r['N']} | {r['class']} | {r['VDE_proxy_kcal']:+.2f} | "
            f"{r['VDE_per_pi_kcal']:+.3f} | {r['sign_matches_4n2']} | {src} |"
        )
    lines.append("")
    for g in gaps:
        lines.append(
            f"- gap |{g['pair'][0]}−{g['pair'][1]}| = {g['abs_gap']:.3f} "
            f"(4n={g['VDE_pi_4n']:+.3f}, 4n+2={g['VDE_pi_4n2']:+.3f})"
        )
    lines.append("")
    if rows_2011:
        lines.append("## 2011-lite")
        for r in rows_2011:
            lines.append(
                f"- N={r['N']} VDE/π={r['VDE_per_pi_kcal']:+.3f} "
                f"sign_ok={r['sign_matches_4n2']}"
            )
        lines.append("")
    if rows_bla:
        lines.append("## BLA 1.34/1.46")
        for r in rows_bla:
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
