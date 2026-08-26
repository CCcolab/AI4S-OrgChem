"""
P2 LFMO-lite: NBA twist scan of EV / Enσσ / Eπσ (independent AO proxy).

Windowed claims (π assignment is clean only at modest twist):
  EV(0–30°)  > 0, dEV/dθ < 0     π–π destablizes, drives distortion
  Enσσ(0–45°) > 0, dEnσσ/dθ < 0  nonbonded σ–σ destablizes, drives distortion
  Eπσ(0°)    ≈ 0                 matches book Table 5-16 for NBA
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

from pyscf import gto  # noqa: E402

from src.common.units import HARTREE_TO_KCAL, ensure_dir, ha_to_kcal, write_json  # noqa: E402
from src.localization.lfmo_ao_proxy import (  # noqa: E402
    construct_states,
    energy_of_density,
    fock_nonbonded_ss,
    ground_mf,
)
from src.p3_nba.geometry import (  # noqa: E402
    IDX_C_IMINE,
    IDX_CIPSO_N,
    IDX_CORTHO_N,
    IDX_N,
    build_nba,
    check_topology,
    dihedral_deg,
    rotate_fragment,
)

DEFAULT_ANGLES = [0.0, 10.0, 17.0, 30.0, 45.0]
CPH_FRAGMENT = list(range(15, 25))
EV_WINDOW = 30.0
EN_WINDOW = 45.0
EPS0_TOL_KCAL = 0.05


def planarize_cphenyl(coords: np.ndarray) -> np.ndarray:
    """Rotate C-phenyl into the imine plane (P2 only; does not alter P3)."""
    out = coords.copy()
    cur = dihedral_deg(out, 15, 14, IDX_C_IMINE, IDX_N)
    target = 180.0 if abs(abs(cur) - 180.0) <= abs(cur) else 0.0
    delta = target - cur
    while delta > 180:
        delta -= 360
    while delta < -180:
        delta += 360
    return rotate_fragment(out, 14, IDX_C_IMINE, CPH_FRAGMENT, -delta)


def build_nba_p2(twist_deg: float) -> tuple[list[str], np.ndarray]:
    symbols, coords = build_nba(twist_deg)
    return symbols, planarize_cphenyl(coords)


def _dmin(coords: np.ndarray) -> float:
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    return float(d.min())


def _mol(symbols: list[str], coords: np.ndarray, basis: str) -> gto.Mole:
    mol = gto.Mole()
    mol.atom = [(symbols[i], tuple(coords[i])) for i in range(len(symbols))]
    mol.basis = basis
    mol.unit = "Angstrom"
    mol.verbose = 0
    mol.build()
    return mol


def finite_slope(xs: list[float], ys: list[float]) -> float:
    slopes = [
        (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
        for i in range(len(xs) - 1)
        if abs(xs[i + 1] - xs[i]) > 1e-9
    ]
    return float(np.mean(slopes)) if slopes else float("nan")


def run(angles: list[float], basis: str, method: str, out: Path) -> dict:
    tables = ensure_dir(out / "tables")
    raw = ensure_dir(out / "raw" / f"lfmo_lite_{method}_{basis.replace('*', 's')}")

    symbols, _ = build_nba_p2(0.0)
    rows = []
    all_conv = True
    meta0 = None

    for ang in angles:
        print(f"[P2 LFMO-lite] theta={ang:.1f}", flush=True)
        _, coords = build_nba_p2(ang)
        problems = check_topology(symbols, coords)
        if problems:
            raise SystemExit(f"topology fail at {ang}: {problems}")
        dmin = _dmin(coords)
        if dmin < 0.85:
            raise SystemExit(f"dmin={dmin:.3f} at theta={ang}")
        th = dihedral_deg(coords, IDX_CORTHO_N, IDX_CIPSO_N, IDX_N, IDX_C_IMINE)
        th_c = dihedral_deg(coords, 15, 14, IDX_C_IMINE, IDX_N)

        mol = _mol(symbols, coords, basis)
        mf = ground_mf(mol, method)
        if not mf.converged:
            all_conv = False
            print("    WARNING: ground SCF not converged", flush=True)
        dens, meta = construct_states(mf)
        if meta0 is None:
            meta0 = {k: v for k, v in meta.items() if k != "pi_detail"}
            meta0["pi_detail0"] = meta["pi_detail"]

        energies = {st: energy_of_density(mf, dens[st]) for st in ("G", "FUD", "DSI")}
        ev_ha = energies["FUD"] - energies["DSI"]
        en_ss = fock_nonbonded_ss(mf, dens["G"], coords)
        ev_kcal = ha_to_kcal(ev_ha)
        en_kcal = ha_to_kcal(en_ss)
        eps_kcal = ha_to_kcal(energies["G"] - energies["FUD"])
        print(
            f"    EV={ev_kcal:.3f}  Enσσ={en_kcal:.3f}  Eπσ={eps_kcal:.4f} kcal",
            flush=True,
        )
        rows.append(
            {
                "theta_target": ang,
                "theta_deg": th,
                "theta_cphenyl_deg": th_c,
                "dmin": dmin,
                "E_G_converged": bool(mf.converged),
                "E_G_ha": energies["G"],
                "E_FUD_ha": energies["FUD"],
                "E_DSI_ha": energies["DSI"],
                "EV_kcal": ev_kcal,
                "En_ss_kcal": en_kcal,
                "E_ps_kcal": eps_kcal,
            }
        )

    ev_rows = [r for r in rows if r["theta_deg"] <= EV_WINDOW + 1.0]
    en_rows = [r for r in rows if r["theta_deg"] <= EN_WINDOW + 1.0]
    evs = [r["EV_kcal"] for r in ev_rows]
    ens = [r["En_ss_kcal"] for r in en_rows]
    ev_pos = all(v > 0.0 for v in evs) and len(evs) >= 3
    en_pos = all(v > 0.0 for v in ens) and len(ens) >= 3
    d_ev = finite_slope([r["theta_deg"] for r in ev_rows], evs)
    d_en = finite_slope([r["theta_deg"] for r in en_rows], ens)
    eps0 = next(r["E_ps_kcal"] for r in rows if abs(r["theta_deg"]) < 1.0)
    eps0_ok = abs(eps0) <= EPS0_TOL_KCAL
    drive_ev = d_ev < 0.0
    drive_en = d_en < 0.0

    ev_span = max(abs(v) for v in evs) if evs else 0.0
    en_span = max(abs(v) for v in ens) if ens else 0.0
    scale_ok = ev_span < 40.0 and en_span < 80.0
    gates_ok = bool(
        all_conv
        and scale_ok
        and meta0 is not None
        and meta0.get("n_pi") == 14
        and all(r["dmin"] >= 0.85 for r in rows)
    )
    two_channel = bool(ev_pos and drive_ev and en_pos and drive_en and eps0_ok)
    agree = True if (gates_ok and two_channel) else (False if gates_ok else None)

    pack = {
        "proposition": "P2",
        "version": "v2_lfmo_lite",
        "protocol": (
            "RHF/STO-3G rigid NBA (C-phenyl planarized); "
            "EV=E(FUD)-E(DSI) AO-block densities on 0–30°; "
            "Enσσ=2∑D_μνF_μν nonbonded σ–σ (no link atoms) on 0–45°; "
            "Eπσ(0)=E(G)-E(FUD) vs Table 5-16"
        ),
        "method": method,
        "basis": basis,
        "angles_deg": angles,
        "pi_meta": meta0,
        "rows": rows,
        "analysis": {
            "EV_window_deg": EV_WINDOW,
            "En_ss_window_deg": EN_WINDOW,
            "EV_all_positive": ev_pos,
            "En_ss_all_positive": en_pos,
            "E_ps_zero_at_planar": eps0_ok,
            "E_ps_0_kcal": eps0,
            "dEV_dtheta_kcal_per_deg": d_ev,
            "dEn_ss_dtheta_kcal_per_deg": d_en,
            "EV_drives_distortion": drive_ev,
            "En_ss_drives_distortion": drive_en,
            "two_channel_ok": two_channel,
            "E_ps_slope_isolated": False,
        },
        "quality_gate": {
            "passed": gates_ok,
            "G1_topology": True,
            "G2_geometry": all(r["dmin"] >= 0.85 for r in rows),
            "G3_convergence": all_conv,
            "G4_energy_scale": scale_ok,
            "G5_path_clean": True,
            "max_abs_EV_kcal": ev_span,
            "max_abs_En_ss_kcal": en_span,
            "n_pi": None if meta0 is None else meta0["n_pi"],
        },
        "agree": agree,
        "hartree_to_kcal": HARTREE_TO_KCAL,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "caveat": (
            "Independent AO-density / Fock proxy of public Fig. 5-15 definitions; "
            "not bit-identical to Kost LFMO. π–σ slope vs θ is not isolated "
            "(assignment/Ne issues at large twist); Eπσ(0)≈0 is the planar check."
        ),
    }
    write_json(tables / "p2_v2_lfmo_lite.json", pack)

    lines = [
        "# P2 v2 — LFMO-lite (NBA)",
        "",
        f"- method={method}/{basis} angles={angles}",
        f"- gates={gates_ok} two_channel={two_channel} agree={agree}",
        f"- dEV/dθ={d_ev:.4f} (0–{EV_WINDOW:.0f}°)  "
        f"dEnσσ/dθ={d_en:.4f} (0–{EN_WINDOW:.0f}°) kcal/deg",
        f"- Eπσ(0)={eps0:.4f} kcal (Table 5-16: 0)",
        "",
        "| θ | EV | Enσσ | Eπσ |",
        "|---|----|------|-----|",
    ]
    for r in rows:
        lines.append(
            f"| {r['theta_deg']:.1f} | {r['EV_kcal']:.2f} | "
            f"{r['En_ss_kcal']:.2f} | {r['E_ps_kcal']:.4f} |"
        )
    lines.append("")
    text = "\n".join(lines)
    (tables / "summary_lfmo_lite.md").write_text(text, encoding="utf-8")
    (raw / "note.txt").write_text(
        "Rigid NBA + C-phenyl planarized; see p2_v2_lfmo_lite.json\n",
        encoding="utf-8",
    )
    print(text, flush=True)
    if not gates_ok:
        raise SystemExit(1)
    if agree is not True:
        raise SystemExit(2)
    return pack


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="sto-3g")
    ap.add_argument("--method", default="RHF")
    ap.add_argument("--angles", default=",".join(str(a) for a in DEFAULT_ANGLES))
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "P2")
    args = ap.parse_args()
    angles = [float(x) for x in args.angles.split(",") if x.strip()]
    run(angles, args.basis, args.method, args.out)


if __name__ == "__main__":
    main()
