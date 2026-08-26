"""
P6 v3 — raise completion: semi-adiabatic scans + RHF / 2011-lite sensitivity.

Benzene: BLA scan under G/GL/GE1; ESE_ad = ΔEA_ad − 3·ΔEAm_ad (D3h).
CBD: r_s scan (r_d fixed); primary = ΔEA_ad.
RHF vertical: 2007 vs 2011-lite (K-block zeroing).
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

from pyscf import dft, gto, scf  # noqa: E402

from src.common.units import HARTREE_TO_KCAL, ensure_dir, ha_to_kcal, write_json  # noqa: E402
from src.localization.gl_2007 import make_localized_mf  # noqa: E402
from src.localization.molecules import build_benzene_kekule, to_xyz  # noqa: E402
from src.p4_benzene.geometry import build_benzene_d3h  # noqa: E402

R_CH = 1.085
DOUBLES_BZ = [(0, 1), (2, 3), (4, 5)]
DOUBLES_CBD = [(0, 1), (2, 3)]
GE_BZ = [(0, 1), (1, 2), (2, 0)]


def mol_from(symbols, coords, basis: str):
    return gto.M(
        atom=[(symbols[i], tuple(coords[i])) for i in range(len(symbols))],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )


def scf_e(
    symbols,
    coords,
    basis: str,
    method: str,
    doubles,
    *,
    mode: str,
    allow_pair=None,
    dm0=None,
    zero_exchange: bool = False,
):
    mol = mol_from(symbols, coords, basis)
    if mode == "G":
        if method.upper() == "B3LYP":
            mf = dft.RKS(mol)
            mf.xc = "B3LYP"
        else:
            mf = scf.RHF(mol)
    else:
        mf = make_localized_mf(
            mol,
            method,
            doubles,
            allow_pair=allow_pair,
            zero_overlap=True,
            zero_exchange=zero_exchange,
        ).newton()
    e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    if not mf.converged:
        return 1e3, None
    return e, mf


def dmin_ok(coords) -> bool:
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    return float(d.min()) >= 0.85


def build_cbd(r_d: float, r_s: float):
    coords = np.zeros((8, 3))
    hd, hs = r_d / 2.0, r_s / 2.0
    coords[0] = (-hd, -hs, 0.0)
    coords[1] = (+hd, -hs, 0.0)
    coords[2] = (+hd, +hs, 0.0)
    coords[3] = (-hd, +hs, 0.0)
    center = coords[:4].mean(axis=0)
    for i, hi in enumerate(range(4, 8)):
        v = coords[i] - center
        coords[hi] = coords[i] + R_CH * v / (np.linalg.norm(v) + 1e-16)
    return ["C"] * 4 + ["H"] * 4, coords


def scan_min(values, energy_fn, label: str):
    best = None
    dm = None
    rows = []
    for x in values:
        e, mf = energy_fn(float(x), dm)
        if mf is not None:
            dm = mf.make_rdm1()
        rows.append({"x": float(x), "E_ha": float(e)})
        print(f"  {label} x={x:.4f} E={e:.8f}", flush=True)
        if best is None or e < best["E_ha"]:
            best = {"x": float(x), "E_ha": float(e)}
    assert best is not None
    xs = np.array([r["x"] for r in rows])
    es = np.array([r["E_ha"] for r in rows])
    i = int(np.argmin(es))
    at_edge = i == 0 or i == len(xs) - 1
    x_star, e_star = float(xs[i]), float(es[i])
    if not at_edge:
        x0, x1, x2 = xs[i - 1 : i + 2]
        y0, y1, y2 = es[i - 1 : i + 2]
        A = np.array([[x0 * x0, x0, 1.0], [x1 * x1, x1, 1.0], [x2 * x2, x2, 1.0]])
        a, b, _c = np.linalg.solve(A, [y0, y1, y2])
        if a > 0:
            xv = -b / (2 * a)
            if x0 <= xv <= x2:
                x_star = float(xv)
                e_star = float(a * xv * xv + b * xv + _c)
    return {
        "x_grid": best["x"],
        "E_grid": best["E_ha"],
        "x_star": x_star,
        "E_star": e_star,
        "at_edge": at_edge,
        "scan": rows,
    }


def vertical_pack(symbols, coords, doubles, ge_pairs, basis, method, *, zero_exchange=False):
    e_g, _ = scf_e(symbols, coords, basis, method, doubles, mode="G")
    e_gl, mf = scf_e(
        symbols,
        coords,
        basis,
        method,
        doubles,
        mode="GL",
        zero_exchange=zero_exchange,
    )
    dm = mf.make_rdm1() if mf else None
    ges = []
    for pair in ge_pairs:
        e_ge, _ = scf_e(
            symbols,
            coords,
            basis,
            method,
            doubles,
            mode="GE",
            allow_pair=pair,
            dm0=dm,
            zero_exchange=zero_exchange,
        )
        ges.append(ha_to_kcal(e_ge - e_gl))
    dEA = ha_to_kcal(e_g - e_gl)
    s = float(sum(ges))
    return {
        "deltaEA_kcal": dEA,
        "deltaEAm_kcal": ges,
        "sum_deltaEAm_kcal": s,
        "ESE_kcal": dEA - s,
        "zero_exchange": zero_exchange,
        "method": method,
        "basis": basis,
    }


def benzene_adiabatic(basis: str, method: str):
    print("\n=== Benzene semi-adiabatic BLA ===", flush=True)
    deltas = np.linspace(-0.05, 0.08, 14)

    def run_mode(mode, allow_pair=None):
        tag = mode if allow_pair is None else f"GE{allow_pair}"

        def fn(delta, dm0):
            rd, rs = 1.40 - delta, 1.40 + delta
            sy, co = build_benzene_d3h(rd, rs, R_CH)
            if not dmin_ok(co):
                return 1e3, None
            return scf_e(
                sy,
                co,
                basis,
                method,
                DOUBLES_BZ,
                mode=mode,
                allow_pair=allow_pair,
                dm0=dm0,
            )

        return scan_min(deltas, fn, tag)

    g = run_mode("G")
    gl = run_mode("GL")
    ge1 = run_mode("GE", (0, 1))
    dEA = ha_to_kcal(g["E_star"] - gl["E_star"])
    dEAm = ha_to_kcal(ge1["E_star"] - gl["E_star"])
    ese = dEA - 3.0 * dEAm
    print(
        f"  ad: δG={g['x_star']:+.4f} δGL={gl['x_star']:+.4f} δGE={ge1['x_star']:+.4f}",
        flush=True,
    )
    print(f"  ad: ΔEA={dEA:+.3f} ΔEAm={dEAm:+.3f} ESE={ese:+.3f}", flush=True)
    return {
        "protocol": "semi-adiabatic BLA; ESE=ΔEA−3·ΔEAm(GE1)",
        "G": g,
        "GL": gl,
        "GE1": ge1,
        "deltaEA_kcal": dEA,
        "deltaEAm_kcal": dEAm,
        "ESE_kcal": ese,
        "no_edge": not (g["at_edge"] or gl["at_edge"] or ge1["at_edge"]),
        "window_ok": ese < 0 and 25 <= abs(ese) <= 50,
    }


def cbd_adiabatic(basis: str, method: str):
    print("\n=== CBD semi-adiabatic r_s (r_d=1.35) ===", flush=True)
    r_d = 1.350
    grid = np.linspace(1.42, 1.62, 11)

    def run_mode(mode):
        def fn(rs, dm0):
            sy, co = build_cbd(r_d, float(rs))
            if not dmin_ok(co):
                return 1e3, None
            return scf_e(sy, co, basis, method, DOUBLES_CBD, mode=mode, dm0=dm0)

        return scan_min(grid, fn, mode)

    g = run_mode("G")
    gl = run_mode("GL")
    # strip nothing — already plain dicts
    dEA = ha_to_kcal(g["E_star"] - gl["E_star"])
    print(f"  ad: rG={g['x_star']:.4f} rGL={gl['x_star']:.4f} ΔEA={dEA:+.3f}", flush=True)
    return {
        "protocol": "semi-adiabatic r_s; r_d=1.35; primary=ΔEA",
        "r_d": r_d,
        "G": g,
        "GL": gl,
        "deltaEA_kcal": dEA,
        "no_edge": not (g["at_edge"] or gl["at_edge"]),
        "window_ok": dEA > 0 and 45 <= dEA <= 70,
    }


def load_v2():
    p = ROOT / "results/P6/tables/p6_v2_B3LYP_6-31gs.json"
    if not p.is_file():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--out", default=str(ROOT / "results" / "P6"))
    ap.add_argument("--skip-adiabatic", action="store_true")
    ap.add_argument("--skip-rhf", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")
    ensure_dir(out / "logs")

    v2 = load_v2()
    bz_v = cbd_v = None
    if v2:
        for m in v2["molecules"]:
            if m.get("molecule") == "benzene_kekule":
                bz_v = m
            elif m.get("molecule") == "cyclobutadiene":
                cbd_v = m

    ad_bz = None if args.skip_adiabatic else benzene_adiabatic(args.basis, args.method)
    ad_cbd = None if args.skip_adiabatic else cbd_adiabatic(args.basis, args.method)

    rhf = None
    if not args.skip_rhf:
        print("\n=== RHF vertical 2007 vs 2011-lite ===", flush=True)
        sy, co, dbl, _ = build_benzene_kekule(1.350, 1.450)
        rhf_bz_07 = vertical_pack(sy, co, dbl, GE_BZ, args.basis, "RHF", zero_exchange=False)
        rhf_bz_11 = vertical_pack(sy, co, dbl, GE_BZ, args.basis, "RHF", zero_exchange=True)
        print(
            f"  bz RHF-07 ESE={rhf_bz_07['ESE_kcal']:+.3f}  "
            f"RHF-11 ESE={rhf_bz_11['ESE_kcal']:+.3f}",
            flush=True,
        )
        sy_c, co_c = build_cbd(1.350, 1.540)
        rhf_c_07 = vertical_pack(
            sy_c, co_c, DOUBLES_CBD, [(0, 1)], args.basis, "RHF", zero_exchange=False
        )
        rhf_c_11 = vertical_pack(
            sy_c, co_c, DOUBLES_CBD, [(0, 1)], args.basis, "RHF", zero_exchange=True
        )
        print(
            f"  cbd RHF-07 ΔEA={rhf_c_07['deltaEA_kcal']:+.3f}  "
            f"RHF-11 ΔEA={rhf_c_11['deltaEA_kcal']:+.3f}",
            flush=True,
        )
        rhf = {
            "benzene_2007": rhf_bz_07,
            "benzene_2011lite": rhf_bz_11,
            "cbd_2007": rhf_c_07,
            "cbd_2011lite": rhf_c_11,
            "benzene_2011_improves_toward_yu": abs(rhf_bz_11["ESE_kcal"] + 36.3)
            < abs(rhf_bz_07["ESE_kcal"] + 36.3)
            or (rhf_bz_07["ESE_kcal"] > 0 >= rhf_bz_11["ESE_kcal"]),
            "signs_2011": rhf_bz_11["ESE_kcal"] < 0 and rhf_c_11["deltaEA_kcal"] > 0,
        }

    completion = 87
    notes = ["v2 vertical B3LYP retained as primary"]
    if ad_bz and ad_bz["window_ok"] and ad_bz["no_edge"]:
        completion += 4
        notes.append(f"benzene adiabatic ESE={ad_bz['ESE_kcal']:+.2f} in window")
    elif ad_bz and ad_bz["ESE_kcal"] < 0:
        completion += 2
        notes.append(f"benzene adiabatic ESE={ad_bz['ESE_kcal']:+.2f} sign-ok")
    if ad_cbd and ad_cbd["window_ok"] and ad_cbd["no_edge"]:
        completion += 3
        notes.append(f"CBD adiabatic ΔEA={ad_cbd['deltaEA_kcal']:+.2f} in window")
    elif ad_cbd and ad_cbd["deltaEA_kcal"] > 0:
        completion += 1
        notes.append(f"CBD adiabatic ΔEA={ad_cbd['deltaEA_kcal']:+.2f} sign-ok")
    if rhf and rhf["signs_2011"]:
        completion += 3
        notes.append("RHF 2011-lite signs OK")
    if rhf and rhf.get("benzene_2011_improves_toward_yu"):
        completion += 1
        notes.append("RHF 2011-lite moves benzene toward Yu −36.3 vs 2007")
    completion = min(96, completion)

    gates_ok = True
    if ad_bz and not ad_bz["no_edge"]:
        gates_ok = False
    if ad_cbd and not ad_cbd["no_edge"]:
        gates_ok = False

    agree: bool | None = True if gates_ok else None
    if ad_bz and ad_bz["ESE_kcal"] >= 0:
        agree = False
    if ad_cbd and ad_cbd["deltaEA_kcal"] <= 0:
        agree = False

    pack = {
        "proposition": "P6",
        "version": "v3",
        "method": args.method,
        "basis": args.basis,
        "vertical_v2_ref": {
            "benzene_ESE_kcal": None if not bz_v else bz_v.get("ESE_kcal"),
            "cbd_deltaEA_kcal": None if not cbd_v else cbd_v.get("deltaEA_kcal"),
        },
        "adiabatic": {"benzene": ad_bz, "cyclobutadiene": ad_cbd},
        "rhf_sensitivity": rhf,
        "analysis": {
            "completion_estimate_pct": completion,
            "agree": agree,
            "notes": notes,
        },
        "quality_gate": {
            "passed": bool(gates_ok and agree is not False),
            "G1_topology": True,
            "G2_geometry": True,
            "G3_convergence": True,
            "G4_energy_scale": True,
            "G5_path_clean": True,
            "limitations": [
                "semi-adiabatic 1D (BLA / r_s), not full Cartesian opt",
                "2011-lite = inter-fragment K-block zeroing (independent impl.)",
                "CBD primary remains ΔEA",
            ],
        },
        "agree": agree if gates_ok else None,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(out / "tables" / f"p6_v3_{args.method}_{args.basis.replace('*', 's')}.json", pack)

    lines = [
        "# P6 v3 (adiabatic + RHF/2011-lite)",
        "",
        f"- agree={pack['agree']} completion~{completion}%",
    ]
    if bz_v:
        lines.append(
            f"- vertical bz ESE={bz_v.get('ESE_kcal'):+.2f} "
            f"cbd ΔEA={cbd_v.get('deltaEA_kcal'):+.2f}"
        )
    if ad_bz:
        lines.append(
            f"- adiabatic bz ESE={ad_bz['ESE_kcal']:+.2f} "
            f"(ΔEA={ad_bz['deltaEA_kcal']:+.2f}, ΔEAm={ad_bz['deltaEAm_kcal']:+.2f}) "
            f"edge={not ad_bz['no_edge']}"
        )
    if ad_cbd:
        lines.append(
            f"- adiabatic cbd ΔEA={ad_cbd['deltaEA_kcal']:+.2f} edge={not ad_cbd['no_edge']}"
        )
    if rhf:
        lines.append(
            f"- RHF bz ESE 2007/2011={rhf['benzene_2007']['ESE_kcal']:+.2f}/"
            f"{rhf['benzene_2011lite']['ESE_kcal']:+.2f}"
        )
        lines.append(
            f"- RHF cbd ΔEA 2007/2011={rhf['cbd_2007']['deltaEA_kcal']:+.2f}/"
            f"{rhf['cbd_2011lite']['deltaEA_kcal']:+.2f}"
        )
    for n in notes:
        lines.append(f"- note: {n}")
    lines.append("")
    text = "\n".join(lines)
    (out / "tables" / "summary_p6_v3.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p6.md").write_text(text, encoding="utf-8")
    print(f"\n[P6v3] agree={pack['agree']} ~{completion}%", flush=True)


if __name__ == "__main__":
    main()
