"""Benzotricyclobutadiene (C12H6) and benzocyclobutadiene (C8H6) builders + Δr metrics."""
from __future__ import annotations

import numpy as np

# Central ring: 0..5; outer doubles: (6,7) on endo 0-1, (8,9) on 2-3, (10,11) on 4-5
ENDO_BT = [(0, 1), (2, 3), (4, 5)]
EXO_BT = [(1, 2), (3, 4), (5, 0)]
ATOMS_A_BT = list(range(6))
ATOMS_B_BT = list(range(6, 12))

# Benzocyclobutadiene: central 0..5 benzene-like; outer double 6-7 on endo 0-1
ENDO_BCB = [(0, 1)]
EXO_BCB = [(1, 2), (5, 0)]  # exo adjacent to the fused bond
ATOMS_A_BCB = list(range(6))
ATOMS_B_BCB = [6, 7]


def _kekule_hexagon(r_long: float, r_short: float) -> np.ndarray:
    """Walk a planar hexagon with alternating bond lengths; center at origin."""
    pts = np.zeros((6, 3))
    x, y, theta = 0.0, 0.0, 0.0
    for i in range(6):
        pts[i, 0] = x
        pts[i, 1] = y
        L = r_long if i % 2 == 0 else r_short
        x += L * np.cos(theta)
        y += L * np.sin(theta)
        theta += np.pi / 3.0
    pts -= pts.mean(axis=0)
    return pts


def _attach_outer_double(
    c_i: np.ndarray,
    c_j: np.ndarray,
    ring_center: np.ndarray,
    r_side: float,
    r_outer: float,
) -> tuple[np.ndarray, np.ndarray]:
    mid = 0.5 * (c_i + c_j)
    bond = c_j - c_i
    bu = bond / (np.linalg.norm(bond) + 1e-16)
    radial = mid - ring_center
    n = np.array([-bu[1], bu[0], 0.0])
    if np.dot(n[:2], radial[:2]) < 0:
        n = -n
    n = n / (np.linalg.norm(n) + 1e-16)
    # Place outer C=C parallel to endo bond at distance ~r_side from mid
    om = mid + n * r_side
    ca = om - 0.5 * r_outer * bu
    cb = om + 0.5 * r_outer * bu
    return ca, cb


def build_benzotricyclobutadiene(
    r_endo: float = 1.50,
    r_exo: float = 1.34,
    r_side: float = 1.46,
    r_outer: float = 1.35,
    r_ch: float = 1.085,
) -> tuple[list[str], np.ndarray]:
    """C12H6 — benzotricyclobutadiene (Yu molecule 10-12)."""
    ring = _kekule_hexagon(r_endo, r_exo)  # even bonds = endo
    center = ring.mean(axis=0)
    outer = []
    for a, b in ENDO_BT:
        ca, cb = _attach_outer_double(ring[a], ring[b], center, r_side, r_outer)
        outer.extend([ca, cb])
    coords = np.vstack([ring, np.array(outer)])
    # H on each outer carbon, pointing further out
    hcoords = []
    for i in range(6, 12):
        v = coords[i] - center
        v = v / (np.linalg.norm(v) + 1e-16)
        hcoords.append(coords[i] + r_ch * v)
    coords = np.vstack([coords, np.array(hcoords)])
    symbols = ["C"] * 12 + ["H"] * 6
    return symbols, coords


def build_benzocyclobutadiene(
    r_endo: float = 1.50,
    r_exo: float = 1.38,
    r_side: float = 1.46,
    r_outer: float = 1.35,
    r_ch: float = 1.085,
) -> tuple[list[str], np.ndarray]:
    """C8H6 — mono-fused extension (sensitivity)."""
    # Nearly equal benzene with one long endo
    ring = _kekule_hexagon(r_endo, r_exo)
    # Make non-fused bonds more benzene-like by mild average
    center = ring.mean(axis=0)
    ca, cb = _attach_outer_double(ring[0], ring[1], center, r_side, r_outer)
    coords = np.vstack([ring, ca[None, :], cb[None, :]])
    h_atoms = []
    # H on C2,C3,C4,C5 (non-bridgehead) and on outer 6,7
    for i in (2, 3, 4, 5):
        v = coords[i] - center
        v = v / (np.linalg.norm(v) + 1e-16)
        h_atoms.append(coords[i] + r_ch * v)
    for i in (6, 7):
        v = coords[i] - center
        v = v / (np.linalg.norm(v) + 1e-16)
        h_atoms.append(coords[i] + r_ch * v)
    coords = np.vstack([coords, np.array(h_atoms)])
    symbols = ["C"] * 8 + ["H"] * 6
    return symbols, coords


def bond_length(coords: np.ndarray, i: int, j: int) -> float:
    return float(np.linalg.norm(coords[i] - coords[j]))


def delta_r(coords: np.ndarray, endo: list[tuple[int, int]], exo: list[tuple[int, int]]) -> dict:
    r_endo = [bond_length(coords, a, b) for a, b in endo]
    r_exo = [bond_length(coords, a, b) for a, b in exo]
    mean_endo = float(np.mean(r_endo))
    mean_exo = float(np.mean(r_exo))
    return {
        "r_endo": r_endo,
        "r_exo": r_exo,
        "r_endo_mean": mean_endo,
        "r_exo_mean": mean_exo,
        "delta_r": mean_endo - mean_exo,
    }


def dmin(coords: np.ndarray) -> float:
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    return float(d.min())


def check_topology_bt(coords: np.ndarray) -> list[str]:
    """Assert C12H6 benzotricyclobutadiene connectivity windows."""
    defects: list[str] = []
    if coords.shape[0] != 18:
        defects.append(f"natm={coords.shape[0]} != 18")
        return defects
    if dmin(coords) < 0.85:
        defects.append(f"dmin={dmin(coords):.3f} < 0.85")
    # endo / exo / outer / radial
    for a, b in ENDO_BT:
        r = bond_length(coords, a, b)
        if not (1.25 <= r <= 1.70):
            defects.append(f"endo {a}-{b}={r:.3f} out of window")
    for a, b in EXO_BT:
        r = bond_length(coords, a, b)
        if not (1.25 <= r <= 1.65):
            defects.append(f"exo {a}-{b}={r:.3f} out of window")
    for a, b in [(6, 7), (8, 9), (10, 11)]:
        r = bond_length(coords, a, b)
        if not (1.20 <= r <= 1.55):
            defects.append(f"outer {a}-{b}={r:.3f} out of window")
    for i, (a, b) in enumerate(ENDO_BT):
        o0, o1 = 6 + 2 * i, 7 + 2 * i
        for c, o in ((a, o0), (b, o1)):
            r = bond_length(coords, c, o)
            if not (1.30 <= r <= 1.70):
                defects.append(f"radial {c}-{o}={r:.3f} out of window")
    # four-membered rings: max diagonal not too short
    for i, (a, b) in enumerate(ENDO_BT):
        o0, o1 = 6 + 2 * i, 7 + 2 * i
        ring4 = [a, b, o1, o0]
        for j in range(4):
            r = bond_length(coords, ring4[j], ring4[(j + 1) % 4])
            if r > 1.85:
                defects.append(f"4-ring bond too long {ring4[j]}-{ring4[(j+1)%4]}={r:.3f}")
    return defects


def check_topology_bcb(coords: np.ndarray) -> list[str]:
    defects: list[str] = []
    if coords.shape[0] != 14:
        defects.append(f"natm={coords.shape[0]} != 14")
        return defects
    if dmin(coords) < 0.85:
        defects.append(f"dmin={dmin(coords):.3f}")
    r_endo = bond_length(coords, 0, 1)
    if not (1.25 <= r_endo <= 1.70):
        defects.append(f"endo={r_endo:.3f}")
    r_out = bond_length(coords, 6, 7)
    if not (1.20 <= r_out <= 1.55):
        defects.append(f"outer={r_out:.3f}")
    return defects
