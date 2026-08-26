"""
P3 — tightened NBA twist PES.

Three stages, addressing the initial run's failure mode (unrelaxed bond
lengths made each theta point descend by a different amount, inflating
conformer dE to ~10^2 kcal/mol):

  A. free pre-relaxation at `geom_basis` (bond lengths/angles become sane)
  B. constrained theta scan from that single relaxed seed, stiff penalty
  C. `energy_method`/`basis` single points at each scan geometry

A quality gate then decides whether the scan may feed a formal verdict:
conformer dE span must be chemically sane and theta errors small.
"""
from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.common.units import HARTREE_TO_KCAL, ensure_dir, ha_to_kcal, write_json  # noqa: E402
from src.p3_nba.geometry import (  # noqa: E402
    IDX_C_IMINE,
    IDX_CIPSO_N,
    IDX_CORTHO_N,
    IDX_N,
    build_nba,
    check_topology,
    dihedral_deg,
    format_xyz,
)
from src.p3_nba.run import (  # noqa: E402
    energy_en,
    optimize_constrained,
    optimize_free,
)

# Reject a geometry outright if any two atoms are closer than this (Angstrom).
MIN_CONTACT_ANG = 0.85
# Conformer energy differences above this are not chemically credible for NBA.
DE_SPAN_LIMIT_KCAL = 25.0
# Constraint is considered honoured below this deviation.
THETA_ERROR_LIMIT_DEG = 3.0
# Up/down sweeps must agree to this much for the profile to count as converged.
HYSTERESIS_LIMIT_KCAL = 1.0


def min_contact(coords: np.ndarray) -> float:
    n = len(coords)
    return min(
        float(np.linalg.norm(coords[i] - coords[j]))
        for i in range(n)
        for j in range(i + 1, n)
    )


def theta_of(coords: np.ndarray) -> float:
    return dihedral_deg(coords, IDX_CORTHO_N, IDX_CIPSO_N, IDX_N, IDX_C_IMINE)


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    n = int(lines[0].split()[0])
    symbols: list[str] = []
    rows: list[list[float]] = []
    for ln in lines[2 : 2 + n]:
        parts = ln.split()
        symbols.append(parts[0])
        rows.append([float(x) for x in parts[1:4]])
    return symbols, np.array(rows, dtype=float)


def constrained_tight(
    symbols: list[str],
    seed: np.ndarray,
    ang: float,
    geom_basis: str,
    geom_method: str,
    maxiter: int,
    k_penalty: float,
    retries: int = 2,
) -> dict:
    """Constrained relax, retried with a stiffer penalty while theta drifts.

    Cartesian BFGS with a quadratic theta penalty cannot enforce the
    constraint exactly; escalating k recovers the target angle without
    making the first (cheapest) pass needlessly stiff.
    """
    k = k_penalty
    pt = optimize_constrained(
        symbols, seed, ang, geom_basis, geom_method, maxiter=maxiter, k_penalty=k
    )
    pt["k_penalty"] = k
    pt["retries"] = 0
    for attempt in range(1, retries + 1):
        if pt["theta_error_deg"] <= THETA_ERROR_LIMIT_DEG:
            break
        k *= 4.0
        print(
            f"      theta err={pt['theta_error_deg']:.2f} > {THETA_ERROR_LIMIT_DEG}; "
            f"retry {attempt} with k={k}",
            flush=True,
        )
        nxt = optimize_constrained(
            symbols, pt["coords"], ang, geom_basis, geom_method,
            maxiter=maxiter, k_penalty=k,
        )
        nxt["k_penalty"] = k
        nxt["retries"] = attempt
        pt = nxt
    return pt


def run_tight(
    geom_basis: str,
    geom_method: str,
    energy_basis: str,
    energy_method: str,
    out: Path,
    angles: list[float],
    cascade: list[tuple[str, int]],
    scan_maxiter: int,
    k_penalty: float,
    seed_xyz: Path | None = None,
    seed_only: bool = False,
    direction: str = "both",
    refine_seed: bool = False,
) -> dict | None:
    tag = (
        f"tight_{energy_method}_{energy_basis.replace('*', 's').replace('+', 'p')}"
        f"_on_{geom_method}_{geom_basis.replace('*', 's').replace('+', 'p')}"
        f"_{direction}"
    )
    raw = ensure_dir(out / "raw" / tag)
    tables = ensure_dir(out / "tables")

    symbols, coords0 = build_nba(40.0)
    d0 = min_contact(coords0)
    print(f"[P3T] built seed: theta={theta_of(coords0):.2f} dmin={d0:.3f} A", flush=True)
    if d0 < MIN_CONTACT_ANG:
        raise SystemExit(f"[P3T] seed geometry rejected: dmin={d0:.3f} < {MIN_CONTACT_ANG}")
    built_defects = check_topology(symbols, coords0)
    if built_defects:
        raise SystemExit(f"[P3T] built geometry is not NBA: {built_defects}")

    # --- Stage A: cascaded free pre-relaxation ----------------------------
    # A crude hand-built geometry is far from equilibrium; relaxing it at a
    # small basis first is cheap and removes most of the bond-length error,
    # so the expensive level only has to polish.
    seed = coords0
    pre_stages: list[dict] = []
    pre = {"theta_deg": theta_of(coords0), "energy_ha": None, "success": False}
    stages = [] if (seed_xyz and not refine_seed) else cascade
    if seed_xyz:
        sym_s, seed = read_xyz(seed_xyz)
        if sym_s != symbols:
            raise SystemExit(f"[P3T] seed xyz atom order mismatch: {seed_xyz}")
        pre = {"theta_deg": theta_of(seed), "energy_ha": None, "success": True}
        print(
            f"[P3T] loaded pre-relaxed seed {seed_xyz} theta={pre['theta_deg']:.2f} "
            f"dmin={min_contact(seed):.3f}",
            flush=True,
        )
    for i, (bas, mit) in enumerate(stages, start=1):
        print(f"[P3T] stage A{i}: free relax {geom_method}/{bas} (maxiter={mit})", flush=True)
        pre = optimize_free(symbols, seed, bas, geom_method, maxiter=mit)
        seed = pre["coords"]
        dmin_i = min_contact(seed)
        print(
            f"[P3T] stage A{i} done: theta={pre['theta_deg']:.2f} E={pre['energy_ha']:.8f} "
            f"nit={pre['nit']} converged={pre['success']} dmin={dmin_i:.3f}",
            flush=True,
        )
        pre_stages.append(
            {
                "basis": bas,
                "maxiter": mit,
                "theta_deg": pre["theta_deg"],
                "E_ha": pre["energy_ha"],
                "nit": pre["nit"],
                "converged": pre["success"],
                "min_contact_ang": dmin_i,
            }
        )
        # Write after every stage so a crash later does not discard the work.
        (out / f"seed_A{i}_{bas.replace('*', 's')}.xyz").write_text(
            f"{len(symbols)}\nstage A{i} {geom_method}/{bas} theta={pre['theta_deg']:.2f}\n"
            + format_xyz(symbols, seed)
            + "\n",
            encoding="utf-8",
        )
        if dmin_i < MIN_CONTACT_ANG:
            raise SystemExit(f"[P3T] stage A{i} produced a clash: dmin={dmin_i:.3f}")
        defects_i = check_topology(symbols, seed)
        if defects_i:
            raise SystemExit(f"[P3T] stage A{i} changed the molecule: {defects_i}")
    dseed = min_contact(seed)
    seed_path = out / "seed_prerelaxed.xyz"
    if stages:
        cascade_desc = " ".join(f"{b}:{m}" for b, m in cascade)
        seed_path.write_text(
            f"{len(symbols)}\nprerelaxed {geom_method} cascade={cascade_desc} "
            f"theta={pre['theta_deg']:.2f}\n" + format_xyz(symbols, seed) + "\n",
            encoding="utf-8",
        )
        write_json(out / "prerelax_stages.json", {"stages": pre_stages})
        print(f"[P3T] seed written to {seed_path}", flush=True)
    if seed_only:
        print("[P3T] --seed-only: stopping before the theta scan", flush=True)
        return None

    # --- Stage B/C: theta sweeps ------------------------------------------
    # Seeding every point from one geometry biases the profile towards the
    # seed's own theta (points far from it stay under-relaxed). Sweeping
    # sequentially, and in both directions, exposes that as hysteresis.
    def sweep(label: str, order: list[float]) -> list[dict]:
        cur = seed
        acc: list[dict] = []
        for ang in order:
            print(f"[P3T] sweep {label}: theta -> {ang:.0f}", flush=True)
            pt = constrained_tight(
                symbols, cur, ang, geom_basis, geom_method, scan_maxiter, k_penalty
            )
            cur = pt["coords"]
            e_sp, en_sp = energy_en(symbols, pt["coords"], energy_basis, energy_method)
            pt["E_geom_ha"] = pt["energy_ha"]
            pt["energy_ha"] = e_sp
            pt["en_ha"] = en_sp
            pt["ee_ha"] = e_sp - en_sp
            pt["min_contact_ang"] = min_contact(pt["coords"])
            pt["defects"] = check_topology(symbols, pt["coords"])
            pt["sweep"] = label
            acc.append(pt)
            (raw / f"nba_{label}_theta_{int(ang):03d}.xyz").write_text(
                f"{len(symbols)}\nsweep={label} target={ang} "
                f"theta={pt['theta_deg']:.2f} {tag}\n"
                + format_xyz(symbols, pt["coords"])
                + "\n",
                encoding="utf-8",
            )
            print(
                f"      theta={pt['theta_deg']:.2f} err={pt['theta_error_deg']:.2f} "
                f"nit={pt['nit']} conv={pt['success']} "
                f"dmin={pt['min_contact_ang']:.3f} defects={len(pt['defects'])} "
                f"E({energy_method})={e_sp:.8f}",
                flush=True,
            )
            if pt["defects"]:
                print(f"      DEFECTS: {pt['defects'][:3]}", flush=True)
        return acc

    asc = sorted(angles)
    sweeps: dict[str, list[dict]] = {}
    if direction in ("up", "both"):
        sweeps["up"] = sweep("up", asc)
    if direction in ("down", "both"):
        sweeps["down"] = sweep("down", list(reversed(asc)))

    by_target = {
        lbl: {r["theta_target"]: r for r in rws} for lbl, rws in sweeps.items()
    }
    hysteresis: list[dict] = []
    if len(sweeps) == 2:
        for ang in asc:
            hysteresis.append(
                {
                    "theta_target": ang,
                    "dE_up_vs_down_kcal": ha_to_kcal(
                        by_target["up"][ang]["energy_ha"]
                        - by_target["down"][ang]["energy_ha"]
                    ),
                }
            )
    # Where both directions exist, keep the lower-energy geometry per angle.
    rows = [
        min((by_target[lbl][ang] for lbl in sweeps), key=lambda r: r["energy_ha"])
        for ang in asc
    ]

    ref = min(rows, key=lambda r: abs(r["theta_deg"]))
    e0, ee0, en0 = ref["energy_ha"], ref["ee_ha"], ref["en_ha"]
    table = [
        {
            "theta_target": r["theta_target"],
            "theta_deg": r["theta_deg"],
            "theta_error_deg": r["theta_error_deg"],
            "E_ha": r["energy_ha"],
            "Ee_ha": r["ee_ha"],
            "EN_ha": r["en_ha"],
            "dE_kcal": ha_to_kcal(r["energy_ha"] - e0),
            "dEe_kcal": ha_to_kcal(r["ee_ha"] - ee0),
            "dEN_kcal": ha_to_kcal(r["en_ha"] - en0),
            "nit": r["nit"],
            "converged": r["success"],
            "min_contact_ang": r["min_contact_ang"],
            "defects": r["defects"],
            "k_penalty": r["k_penalty"],
            "retries": r["retries"],
            "sweep": r["sweep"],
        }
        for r in rows
    ]

    de_values = [t["dE_kcal"] for t in table]
    de_span = float(max(de_values) - min(de_values))
    max_theta_err = float(max(t["theta_error_deg"] for t in table))
    gate = {
        "dE_span_kcal": de_span,
        "dE_span_limit_kcal": DE_SPAN_LIMIT_KCAL,
        "dE_span_ok": de_span <= DE_SPAN_LIMIT_KCAL,
        "max_theta_error_deg": max_theta_err,
        "theta_error_limit_deg": THETA_ERROR_LIMIT_DEG,
        "theta_error_ok": max_theta_err <= THETA_ERROR_LIMIT_DEG,
        "min_contact_ang": float(min(t["min_contact_ang"] for t in table)),
        "min_contact_ok": min(t["min_contact_ang"] for t in table) >= MIN_CONTACT_ANG,
        "topology_ok": all(not t["defects"] for t in table),
        "n_converged": sum(1 for t in table if t["converged"]),
        "n_points": len(table),
        "all_converged": all(t["converged"] for t in table),
    }
    # A profile where every point merely ran out of iterations proves nothing.
    # Either each point met the gradient threshold, or the two sweep
    # directions must agree, which is the practical convergence test.
    if hysteresis:
        max_hyst = max(abs(h["dE_up_vs_down_kcal"]) for h in hysteresis)
        gate["max_hysteresis_kcal"] = max_hyst
        gate["hysteresis_limit_kcal"] = HYSTERESIS_LIMIT_KCAL
        gate["hysteresis_ok"] = max_hyst <= HYSTERESIS_LIMIT_KCAL
    else:
        gate["max_hysteresis_kcal"] = None
        gate["hysteresis_ok"] = False
    gate["convergence_ok"] = bool(gate["all_converged"] or gate["hysteresis_ok"])
    gate["passed"] = bool(
        gate["dE_span_ok"]
        and gate["theta_error_ok"]
        and gate["min_contact_ok"]
        and gate["convergence_ok"]
        and gate["topology_ok"]
    )

    tmin = table[int(np.argmin(de_values))]
    abs_th = abs(tmin["theta_deg"])
    if abs_th > 90:
        abs_th = 180 - abs_th
    planar = [t["dE_kcal"] for t in table if abs(t["theta_deg"]) < 20.0]
    checks = {
        "E_min_in_30_60": bool(25.0 <= abs_th <= 65.0),
        "at_min_EN_up_Ee_down": bool(tmin["dEN_kcal"] > 0.0 and tmin["dEe_kcal"] < 0.0),
        "near_planar_not_global_min": bool(
            tmin["dE_kcal"] < min(planar) - 0.05 if planar else True
        ),
    }
    # Only a gate-passing scan may carry an agree/disagree signal.
    agree = (
        bool(checks["E_min_in_30_60"] and checks["at_min_EN_up_Ee_down"])
        if gate["passed"]
        else None
    )

    pack = {
        "proposition": "P3",
        "tag": tag,
        "protocol": (
            "A: cascaded free relax "
            + " -> ".join(f"{geom_method}/{b}(maxiter={m})" for b, m in cascade)
            + f"; B: sequential theta sweep ({direction}) at {geom_method}/{geom_basis} "
            f"(maxiter={scan_maxiter}, k0={k_penalty} Ha/rad^2, "
            f"stiffening retries until |dtheta|<={THETA_ERROR_LIMIT_DEG} deg, "
            "each point seeded from the previous one); "
            f"C: {energy_method}/{energy_basis} single points; "
            "per-angle energy = lower of the two sweeps"
        ),
        "direction": direction,
        "hysteresis": hysteresis,
        "geom_method": geom_method,
        "geom_basis": geom_basis,
        "energy_method": energy_method,
        "energy_basis": energy_basis,
        "theta_definition": "C_ortho(N-Ph)-C_ipso(N-Ph)-N-C_imine",
        "book_theta_exp_deg": [36, 55],
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "prerelax_stages": pre_stages,
        "prerelax_seed": {
            "theta_deg": pre["theta_deg"],
            "E_ha": pre["energy_ha"],
            "converged": pre["success"],
            "min_contact_ang": dseed,
        },
        "reference_theta_deg": ref["theta_deg"],
        "scan": table,
        "E_min_scan": {
            "theta_deg": tmin["theta_deg"],
            "dE_kcal": tmin["dE_kcal"],
            "dEe_kcal": tmin["dEe_kcal"],
            "dEN_kcal": tmin["dEN_kcal"],
        },
        "quality_gate": gate,
        "checks": checks,
        "agree": agree,
        "verdict_eligible": gate["passed"],
        "hartree_to_kcal": HARTREE_TO_KCAL,
    }
    write_json(tables / f"pes_{tag}.json", pack)

    with (tables / f"pes_{tag}.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "theta_target",
                "theta_deg",
                "theta_error_deg",
                "dE_kcal",
                "dEe_kcal",
                "dEN_kcal",
                "nit",
                "converged",
                "k_penalty",
                "retries",
                "sweep",
            ],
        )
        w.writeheader()
        for t in table:
            w.writerow(
                {
                    "theta_target": f"{t['theta_target']:.1f}",
                    "theta_deg": f"{t['theta_deg']:.3f}",
                    "theta_error_deg": f"{t['theta_error_deg']:.3f}",
                    "dE_kcal": f"{t['dE_kcal']:.4f}",
                    "dEe_kcal": f"{t['dEe_kcal']:.4f}",
                    "dEN_kcal": f"{t['dEN_kcal']:.4f}",
                    "nit": t["nit"],
                    "converged": "yes" if t["converged"] else "no",
                    "k_penalty": f"{t['k_penalty']:.2f}",
                    "retries": t["retries"],
                    "sweep": t["sweep"],
                }
            )

    lines = [
        f"# P3 NBA tightened twist PES ({tag})",
        "",
        f"- Protocol: {pack['protocol']}",
        f"- Pre-relaxed seed: theta={pre['theta_deg']:.1f} deg, "
        f"converged={pre['success']}, dmin={dseed:.3f} A",
        f"- E_min at theta = **{tmin['theta_deg']:.1f} deg** "
        f"(dE={tmin['dE_kcal']:+.3f}, dEe={tmin['dEe_kcal']:+.3f}, "
        f"dEN={tmin['dEN_kcal']:+.3f} kcal/mol vs theta={ref['theta_deg']:.1f})",
        f"- dE span = {de_span:.2f} kcal/mol (limit {DE_SPAN_LIMIT_KCAL})",
        f"- max |theta - target| = {max_theta_err:.2f} deg (limit {THETA_ERROR_LIMIT_DEG})",
        f"- topology intact at every point: {gate['topology_ok']}",
        f"- converged points: {gate['n_converged']}/{gate['n_points']}; "
        f"max up/down hysteresis = {gate['max_hysteresis_kcal']} kcal/mol "
        f"(limit {HYSTERESIS_LIMIT_KCAL})",
        f"- **Quality gate: {'PASSED' if gate['passed'] else 'FAILED'}**",
        f"- Checks: {checks}",
        f"- Auto agree flag: **{agree}** (None = gate failed; formal VERDICT only in deliverables/)",
        "",
    ]
    (out / "summary_tight.md").write_text("\n".join(lines), encoding="utf-8")
    write_json(
        out / "meta.json",
        {"proposition": "P3", "updated_utc": pack["timestamp_utc"], "latest_tag": tag},
    )
    print("\n".join(lines), flush=True)
    return pack


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geom-basis", default="3-21g")
    ap.add_argument("--geom-method", default="RHF", choices=["B3LYP", "RHF"])
    ap.add_argument("--energy-basis", default="6-31g*")
    ap.add_argument("--energy-method", default="B3LYP", choices=["B3LYP", "RHF"])
    ap.add_argument("--out", default=str(ROOT / "results" / "P3"))
    ap.add_argument("--angles", default="0,15,30,45,60,75,90")
    ap.add_argument(
        "--cascade",
        default="sto-3g:250,3-21g:80",
        help="comma list of basis:maxiter free pre-relaxation stages",
    )
    ap.add_argument("--scan-maxiter", type=int, default=30)
    ap.add_argument("--k-penalty", type=float, default=2.0)
    ap.add_argument("--seed-only", action="store_true", help="stop after pre-relaxation")
    ap.add_argument("--seed-xyz", default=None, help="reuse a pre-relaxed seed geometry")
    ap.add_argument("--direction", default="both", choices=["up", "down", "both"])
    ap.add_argument(
        "--refine-seed",
        action="store_true",
        help="apply the cascade stages starting from --seed-xyz",
    )
    args = ap.parse_args()

    out = Path(args.out)
    ensure_dir(out)
    angles = [float(x) for x in args.angles.split(",") if x.strip()]
    cascade = []
    for item in args.cascade.split(","):
        if not item.strip():
            continue
        bas, _, mit = item.partition(":")
        cascade.append((bas.strip(), int(mit)))
    run_tight(
        args.geom_basis,
        args.geom_method,
        args.energy_basis,
        args.energy_method,
        out,
        angles,
        cascade,
        args.scan_maxiter,
        args.k_penalty,
        Path(args.seed_xyz) if args.seed_xyz else None,
        args.seed_only,
        args.direction,
        args.refine_seed,
    )


if __name__ == "__main__":
    main()
