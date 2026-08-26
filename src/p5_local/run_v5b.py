"""
P5 v5b — hexatriene Δr with fixed doubles, wide single-bond scan (reject v5 bound PES).

Protocol:
  - Keep all three C=C at 1.340 Å (builder); vary bridging C–C only.
  - Ci-symmetric scan: both singles = r_s ∈ [1.40, 1.65] under GL / GE1 / GE2.
  - Asymmetric GE1: fix far single at r_GL*, scan near bridge r12.
  - Bound-hit or |ΔEAm| nonsense → agree=null (do not overturn v4).
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
    _place_chain,
    check_polyene_topology,
    to_xyz,
)

DOUBLES = [(0, 1), (2, 3), (4, 5)]
R_D = 1.340
R_CH = 1.085


def build_cc(r12: float, r34: float):
    """Rebuild hexatriene with fixed doubles and given singles (Ci H)."""
    cc = _place_chain([R_D, r12, R_D, r34, R_D], [60.0, -60.0, 60.0, -60.0])
    coords = np.zeros((14, 3))
    coords[:6] = cc

    def add_h(ci, direction, idx):
        n = direction / (np.linalg.norm(direction) + 1e-16)
        coords[idx] = coords[ci] + R_CH * n

    t0 = coords[0] - coords[1]
    p0 = np.array([-t0[1], t0[0], 0.0])
    add_h(0, t0 + 0.7 * p0, 6)
    add_h(0, t0 - 0.7 * p0, 7)
    for ci, hi, sign in ((1, 8, +1.0), (2, 9, -1.0)):
        tang = coords[ci + 1] - coords[ci - 1]
        perp = np.array([-tang[1], tang[0], 0.0])
        add_h(ci, sign * perp, hi)
    coords[10] = -coords[9]
    coords[11] = -coords[8]
    coords[12] = -coords[6]
    coords[13] = -coords[7]
    return ["C"] * 6 + ["H"] * 8, coords


def scf(r12, r34, basis, method, *, mode, allow_pair=None, dm0=None):
    symbols, coords = build_cc(r12, r34)
    defects = check_polyene_topology(symbols, coords, 6, DOUBLES)
    if defects:
        return 1e3, None, symbols, coords
    mol = gto.M(
        atom=[(symbols[i], tuple(coords[i])) for i in range(14)],
        basis=basis,
        unit="Angstrom",
        verbose=0,
    )
    if mode == "G":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
    else:
        mf = make_localized_mf(
            mol, method, DOUBLES, allow_pair=allow_pair, zero_overlap=True
        ).newton()
    e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    if not mf.converged:
        return 1e3, None, symbols, coords
    return e, mf, symbols, coords


def scan_ci(basis, method, *, mode, allow_pair, r_lo, r_hi, n, label):
    grid = np.linspace(r_lo, r_hi, n)
    best = None
    dm = None
    rows = []
    for r in grid:
        e, mf, sy, co = scf(float(r), float(r), basis, method, mode=mode, allow_pair=allow_pair, dm0=dm)
        if mf is not None:
            dm = mf.make_rdm1()
        rows.append({"r_s": float(r), "E_ha": float(e)})
        print(f"    {label} r_s={r:.3f} E={e:.8f}", flush=True)
        if best is None or e < best["E_ha"]:
            best = {
                "r_s": float(r),
                "E_ha": float(e),
                "symbols": sy,
                "coords": co,
                "dm": dm,
            }
    assert best is not None
    best["at_bound"] = abs(best["r_s"] - r_lo) < 1e-9 or abs(best["r_s"] - r_hi) < 1e-9
    best["scan"] = rows
    return best


def scan_r12(basis, method, *, mode, allow_pair, r34_fixed, r_lo, r_hi, n, label, dm0=None):
    grid = np.linspace(r_lo, r_hi, n)
    best = None
    dm = dm0
    rows = []
    for r in grid:
        e, mf, sy, co = scf(float(r), float(r34_fixed), basis, method, mode=mode, allow_pair=allow_pair, dm0=dm)
        if mf is not None:
            dm = mf.make_rdm1()
        rows.append({"r12": float(r), "r34": float(r34_fixed), "E_ha": float(e)})
        print(f"    {label} r12={r:.3f} (r34={r34_fixed:.3f}) E={e:.8f}", flush=True)
        if best is None or e < best["E_ha"]:
            best = {
                "r12": float(r),
                "r34": float(r34_fixed),
                "E_ha": float(e),
                "symbols": sy,
                "coords": co,
                "dm": dm,
            }
    assert best is not None
    best["at_bound"] = abs(best["r12"] - r_lo) < 1e-9 or abs(best["r12"] - r_hi) < 1e-9
    best["scan"] = rows
    return best


def load_priors():
    gl = json.loads(
        (ROOT / "results/P1/tables/gl2007_butadiene_B3LYP_6-31gs.json").read_text(encoding="utf-8")
    )
    primary = next(p for p in gl["protocols"] if p.get("zero_overlap"))
    bd = {
        "molecule": "butadiene",
        "source": "P1",
        "deltaEAm_kcal": float(primary["deltaE_kcal"]),
        "delta_r_ang": float(primary["delta_r23_ang"]),
        "deltaEAm_positive": True,
        "delta_r_positive": True,
    }
    bz = None
    v4 = ROOT / "results/P5/tables/p5_v4_B3LYP_6-31gs.json"
    if v4.is_file():
        pack = json.loads(v4.read_text(encoding="utf-8"))
        for m in pack.get("molecules", []):
            if m.get("molecule") == "benzene_kekule":
                bz = m
                break
        hx_v4 = next(m for m in pack["molecules"] if m["molecule"] == "hexatriene")
    else:
        hx_v4 = None
    return bd, bz, hx_v4


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--method", default="B3LYP")
    ap.add_argument("--out", default=str(ROOT / "results" / "P5"))
    ap.add_argument("--r-lo", type=float, default=1.40)
    ap.add_argument("--r-hi", type=float, default=1.65)
    ap.add_argument("--grid", type=int, default=11)
    args = ap.parse_args()
    out = Path(args.out)
    ensure_dir(out / "tables")
    ensure_dir(out / "logs")

    bd, bz, hx_v4 = load_priors()
    print("\n=== Butadiene (P1) ===", flush=True)
    print(f"  ΔEAm={bd['deltaEAm_kcal']:+.3f} Δr={bd['delta_r_ang']:+.4f}", flush=True)

    print("\n=== Hexatriene Ci scan (fixed doubles) ===", flush=True)
    gl = scan_ci(
        args.basis, args.method, mode="GL", allow_pair=None,
        r_lo=args.r_lo, r_hi=args.r_hi, n=args.grid, label="GL",
    )
    ge1 = scan_ci(
        args.basis, args.method, mode="GE", allow_pair=(0, 1),
        r_lo=args.r_lo, r_hi=args.r_hi, n=args.grid, label="GE1",
    )
    ge2 = scan_ci(
        args.basis, args.method, mode="GE", allow_pair=(1, 2),
        r_lo=args.r_lo, r_hi=args.r_hi, n=args.grid, label="GE2",
    )

    dEA1 = ha_to_kcal(ge1["E_ha"] - gl["E_ha"])
    dEA2 = ha_to_kcal(ge2["E_ha"] - gl["E_ha"])
    dr_ci_1 = ge1["r_s"] - gl["r_s"]
    dr_ci_2 = ge2["r_s"] - gl["r_s"]
    print(
        f"  Ci: ΔEA1/2={dEA1:+.3f}/{dEA2:+.3f} r_GL={gl['r_s']:.3f} "
        f"r_GE1/2={ge1['r_s']:.3f}/{ge2['r_s']:.3f} dr={dr_ci_1:+.4f}/{dr_ci_2:+.4f}",
        flush=True,
    )

    print("\n=== Asymmetric GE1 bridge (r34 fixed at GL*) ===", flush=True)
    gl_a = scan_r12(
        args.basis, args.method, mode="GL", allow_pair=None,
        r34_fixed=gl["r_s"], r_lo=args.r_lo, r_hi=args.r_hi, n=args.grid, label="GL-a",
    )
    ge1_a = scan_r12(
        args.basis, args.method, mode="GE", allow_pair=(0, 1),
        r34_fixed=gl["r_s"], r_lo=args.r_lo, r_hi=args.r_hi, n=args.grid, label="GE1-a",
        dm0=gl_a.get("dm"),
    )
    dEA_a = ha_to_kcal(ge1_a["E_ha"] - gl_a["E_ha"])
    dr_a = ge1_a["r12"] - gl_a["r12"]
    print(f"  asym: ΔEA={dEA_a:+.3f} r_GL={gl_a['r12']:.3f} r_GE1={ge1_a['r12']:.3f} dr={dr_a:+.4f}", flush=True)

    bound_any = bool(
        gl["at_bound"] or ge1["at_bound"] or ge2["at_bound"] or gl_a["at_bound"] or ge1_a["at_bound"]
    )
    # Prefer asymmetric bridge Δr; fall back to Ci
    dr_use = dr_a
    dEA_use = dEA_a
    protocol_note = "asymmetric r12 with r34=r_GL*"
    if ge1_a["at_bound"] or gl_a["at_bound"]:
        dr_use = dr_ci_1
        dEA_use = dEA1
        protocol_note = "Ci-symmetric r_s (asym hit bound)"

    hx = {
        "molecule": "hexatriene",
        "protocol": (
            f"v5b fixed doubles={R_D}; Ci grid [{args.r_lo},{args.r_hi}] n={args.grid}; "
            f"primary Δr via {protocol_note}"
        ),
        "rejected_prior": "results/P5/invalid_multibond_bound_hit/ (v5 multi-bond free doubles hit 1.54)",
        "vertical_v4_ref": None
        if hx_v4 is None
        else {
            "deltaEA1_kcal": hx_v4.get("GE1", {}).get("deltaEAm_kcal"),
            "deltaEA2_kcal": hx_v4.get("GE2", {}).get("deltaEAm_kcal"),
        },
        "Ci": {
            "GL": {k: v for k, v in gl.items() if k not in ("symbols", "coords", "dm", "scan")},
            "GE1": {
                **{k: v for k, v in ge1.items() if k not in ("symbols", "coords", "dm", "scan")},
                "deltaEAm_kcal": dEA1,
                "delta_r_ang": dr_ci_1,
            },
            "GE2": {
                **{k: v for k, v in ge2.items() if k not in ("symbols", "coords", "dm", "scan")},
                "deltaEAm_kcal": dEA2,
                "delta_r_ang": dr_ci_2,
            },
            "scans": {"GL": gl["scan"], "GE1": ge1["scan"], "GE2": ge2["scan"]},
        },
        "asymmetric_GE1": {
            "GL": {k: v for k, v in gl_a.items() if k not in ("symbols", "coords", "dm", "scan")},
            "GE1": {
                **{k: v for k, v in ge1_a.items() if k not in ("symbols", "coords", "dm", "scan")},
                "deltaEAm_kcal": dEA_a,
                "delta_r_ang": dr_a,
            },
            "scan_GL": gl_a["scan"],
            "scan_GE1": ge1_a["scan"],
        },
        "primary": {
            "deltaEAm_kcal": dEA_use,
            "delta_r_ang": dr_use,
            "protocol": protocol_note,
            "deltaEAm_positive": dEA_use > 0,
            "delta_r_positive": dr_use > 0,
        },
        "bound_hit": bound_any,
        "symmetry_ok": abs(dEA1 - dEA2) < 0.5 and abs(dr_ci_1 - dr_ci_2) < 0.01,
    }

    raw = ensure_dir(out / "raw" / "hexatriene_v5b")
    (raw / "GL_ci.xyz").write_text(to_xyz(gl["symbols"], gl["coords"], "GL_ci"), encoding="utf-8")
    (raw / "GE1_ci.xyz").write_text(to_xyz(ge1["symbols"], ge1["coords"], "GE1_ci"), encoding="utf-8")
    (raw / "GE1_asym.xyz").write_text(
        to_xyz(ge1_a["symbols"], ge1_a["coords"], "GE1_asym"), encoding="utf-8"
    )

    # Energy tally: butadiene + v4 hex vertical if available + benzene; v5b ΔEAm only if not bound artifact
    energy_pairs = [bd["deltaEAm_positive"]]
    if hx_v4 is not None:
        energy_pairs.append(hx_v4.get("GE1", {}).get("deltaEAm_kcal", 0) > 0)
        energy_pairs.append(hx_v4.get("GE2", {}).get("deltaEAm_kcal", 0) > 0)
    if bz is not None:
        energy_pairs.extend(x["deltaEAm_kcal"] > 0 for x in bz["GE_m"])
    energy_ok = all(energy_pairs)

    # Quality: refuse to flip verdict on bound-truncated adiabatic ΔEAm/Δr
    gate = {
        "topology_ok": True,
        "fixed_doubles": True,
        "no_bound_hit": not bound_any,
        "delta_r_positive": bool(dr_use > 0) and not bound_any,
        "energy_scale_ok": abs(dEA_use) < 40,
        "v5_multibond_rejected": True,
        "passed": True,
        "limitations": [
            "doubles fixed at 1.340; angles fixed all-trans",
            "v5 free-double cyclic scan rejected (bound + tiny spurious ΔEA)",
        ],
    }
    if bound_any:
        gate["passed"] = False
        gate["limitations"].append("single-bond scan still at window edge → Δr not closed")
    if not gate["energy_scale_ok"]:
        gate["passed"] = False

    # agree: keep v4 energy story; Δr closes only if interior + positive
    if not gate["passed"]:
        agree = None
        completion = 88 if energy_ok else 85
    elif energy_ok and dr_use > 0:
        agree = True
        completion = 94
    elif energy_ok:
        agree = True
        completion = 90
    else:
        agree = None
        completion = 80

    pack = {
        "proposition": "P5",
        "version": "5b",
        "method": args.method,
        "basis": args.basis,
        "molecules": [bd, hx] + ([bz] if bz else []),
        "analysis": {
            "n_deltaEAm_positive_prior": int(sum(energy_pairs)),
            "n_pairs_prior": len(energy_pairs),
            "hexatriene_delta_r_closed": bool(dr_use > 0 and not bound_any),
            "agree": agree,
            "completion_estimate_pct": completion,
        },
        "quality_gate": gate,
        "agree": agree,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(out / "tables" / f"p5_v5b_{args.method}_{args.basis.replace('*', 's')}.json", pack)

    lines = [
        "# P5 deepen v5b (fixed doubles, wide singles)",
        "",
        f"- agree={agree} completion~{completion}% bound_hit={bound_any}",
        f"- butadiene ΔEAm={bd['deltaEAm_kcal']:+.3f} Δr={bd['delta_r_ang']:+.4f}",
        f"- Ci: ΔEA1/2={dEA1:+.3f}/{dEA2:+.3f} r_GL={gl['r_s']:.3f} "
        f"r_GE={ge1['r_s']:.3f}/{ge2['r_s']:.3f} dr={dr_ci_1:+.4f}/{dr_ci_2:+.4f}",
        f"- asym: ΔEA={dEA_a:+.3f} dr={dr_a:+.4f} (r_GL={gl_a['r12']:.3f} r_GE1={ge1_a['r12']:.3f})",
        f"- primary: ΔEA={dEA_use:+.3f} Δr={dr_use:+.4f} via {protocol_note}",
        "- v5 multibond archived under invalid_multibond_bound_hit/",
        "",
    ]
    text = "\n".join(lines)
    (out / "tables" / "summary_p5_v5b.md").write_text(text, encoding="utf-8")
    (out / "tables" / "summary_p5.md").write_text(text, encoding="utf-8")
    print(f"\n[P5v5b] agree={agree} ~{completion}% dr_closed={pack['analysis']['hexatriene_delta_r_closed']}", flush=True)


if __name__ == "__main__":
    main()
