"""Five-membered heterocycle builders (furan, pyrrole) for P8 LDE."""
from __future__ import annotations

import numpy as np

R_CH = 1.085
R_NH = 1.01


def _walk_ring(
    bond_lengths: list[float],
    start_angle: float = np.pi / 2.0,
) -> np.ndarray:
    """Place n atoms on a closed ring in the xy-plane."""
    n = len(bond_lengths)
    coords = np.zeros((n, 3))
    x = y = 0.0
    angle = start_angle
    step = 2.0 * np.pi / n
    for i in range(n):
        coords[i, 0] = x
        coords[i, 1] = y
        bl = bond_lengths[i]
        x += bl * np.cos(angle)
        y += bl * np.sin(angle)
        angle -= step
    coords -= coords.mean(axis=0)
    coords[:, 2] = 0.0
    return coords


def _add_h_on_c(coords: np.ndarray, ci: int, ring_center: np.ndarray, r_ch: float = R_CH) -> np.ndarray:
    v = coords[ci] - ring_center
    v = v / (np.linalg.norm(v) + 1e-16)
    return coords[ci] + r_ch * v


def build_furan(
    r_co: float = 1.360,
    r_cc_d: float = 1.340,
    r_cc_s: float = 1.430,
    r_ch: float = R_CH,
) -> tuple[list[str], np.ndarray, list[tuple[int, int]], list[tuple[int, int]]]:
    """Planar furan. Atom 0=O; ring 0-1-2-3-4-0. Doubles (1,2),(3,4)."""
    ring = _walk_ring([r_co, r_cc_d, r_cc_s, r_cc_d, r_co])
    center = ring.mean(axis=0)
    hcoords = [_add_h_on_c(ring, i, center, r_ch) for i in range(1, 5)]
    coords = np.vstack([ring, np.array(hcoords)])
    symbols = ["O", "C", "C", "C", "C"] + ["H"] * 4
    doubles = [(1, 2), (3, 4)]
    singles = [(0, 1), (2, 3), (4, 0)]
    return symbols, coords, doubles, singles


def build_pyrrole(
    r_nh: float = 1.010,
    r_cc_d: float = 1.340,
    r_cc_s: float = 1.430,
    r_cn: float = 1.370,
    r_ch: float = R_CH,
) -> tuple[list[str], np.ndarray, list[tuple[int, int]], list[tuple[int, int]]]:
    """Planar pyrrole. Atom 0=N; ring 0-1-2-3-4-0. Doubles (1,2),(3,4)."""
    ring = _walk_ring([r_cn, r_cc_d, r_cc_s, r_cc_d, r_cn])
    center = ring.mean(axis=0)
    # N-H outward from ring
    nh = ring[0] + r_nh * (ring[0] - center) / (np.linalg.norm(ring[0] - center) + 1e-16)
    h_c = [_add_h_on_c(ring, i, center, r_ch) for i in range(1, 5)]
    coords = np.vstack([ring, nh[None, :], np.array(h_c)])
    symbols = ["N", "C", "C", "C", "C", "H"] + ["H"] * 4
    doubles = [(1, 2), (3, 4)]
    singles = [(0, 1), (2, 3), (4, 0)]
    return symbols, coords, doubles, singles


def bond_length(coords: np.ndarray, i: int, j: int) -> float:
    return float(np.linalg.norm(coords[i] - coords[j]))


def dmin(coords: np.ndarray) -> float:
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    return float(d.min())


def check_furan_topology(coords: np.ndarray) -> list[str]:
    defects: list[str] = []
    if coords.shape[0] != 9:
        defects.append(f"natm={coords.shape[0]}!=9")
        return defects
    if dmin(coords) < 0.85:
        defects.append(f"dmin={dmin(coords):.3f}")
    for a, b in [(1, 2), (3, 4)]:
        r = bond_length(coords, a, b)
        if not (1.20 <= r <= 1.55):
            defects.append(f"double {a}-{b}={r:.3f}")
    for a, b in [(0, 1), (0, 4), (2, 3)]:
        r = bond_length(coords, a, b)
        if not (1.20 <= r <= 1.65):
            defects.append(f"single {a}-{b}={r:.3f}")
    if float(np.max(np.abs(coords[:, 2]))) > 0.05:
        defects.append("non-planar")
    return defects


def check_pyrrole_topology(coords: np.ndarray) -> list[str]:
    defects: list[str] = []
    if coords.shape[0] != 10:
        defects.append(f"natm={coords.shape[0]}!=10")
        return defects
    if dmin(coords) < 0.85:
        defects.append(f"dmin={dmin(coords):.3f}")
    for a, b in [(1, 2), (3, 4)]:
        r = bond_length(coords, a, b)
        if not (1.20 <= r <= 1.55):
            defects.append(f"double {a}-{b}={r:.3f}")
    r_nh = bond_length(coords, 0, 5)
    if not (0.95 <= r_nh <= 1.15):
        defects.append(f"N-H={r_nh:.3f}")
    return defects

def build_oxazole(
    r_co: float = 1.360,
    r_cn: float = 1.300,
    r_cc_d: float = 1.350,
    r_nc: float = 1.390,
    r_cc_s: float = 1.360,
    r_ch: float = R_CH,
) -> tuple[list[str], np.ndarray, list[tuple[int, int]], list[tuple[int, int]]]:
    """Planar oxazole. 0=O, 1=C2, 2=N3, 3=C4, 4=C5. Doubles (1,2)=C=N, (3,4)=C=C."""
    # bond order around ring: O-C2, C2=N3, N3-C4, C4=C5, C5-O
    ring = _walk_ring([r_co, r_cn, r_nc, r_cc_d, r_cc_s])
    center = ring.mean(axis=0)
    # H on C2, C4, C5 (N has no H)
    hcoords = [
        _add_h_on_c(ring, 1, center, r_ch),
        _add_h_on_c(ring, 3, center, r_ch),
        _add_h_on_c(ring, 4, center, r_ch),
    ]
    coords = np.vstack([ring, np.array(hcoords)])
    symbols = ["O", "C", "N", "C", "C"] + ["H"] * 3
    doubles = [(1, 2), (3, 4)]
    singles = [(0, 1), (2, 3), (4, 0)]
    return symbols, coords, doubles, singles


def check_oxazole_topology(coords: np.ndarray) -> list[str]:
    defects: list[str] = []
    if coords.shape[0] != 8:
        defects.append(f"natm={coords.shape[0]}!=8")
        return defects
    if dmin(coords) < 0.85:
        defects.append(f"dmin={dmin(coords):.3f}")
    for a, b in [(1, 2), (3, 4)]:
        r = bond_length(coords, a, b)
        if not (1.20 <= r <= 1.55):
            defects.append(f"double {a}-{b}={r:.3f}")
    for a, b in [(0, 1), (0, 4), (2, 3)]:
        r = bond_length(coords, a, b)
        if not (1.20 <= r <= 1.65):
            defects.append(f"single {a}-{b}={r:.3f}")
    if float(np.max(np.abs(coords[:, 2]))) > 0.05:
        defects.append("non-planar")
    return defects


# Fragment layout for hetero GL (Ch6 Fig 6-8): [double1, double2, heteroatom]
HETERO_PI_ATOM = {"furan": [[0]], "pyrrole": [[0]], "oxazole": [[0]]}
# GE-1: two doubles; GE-2/3: hetero + each double
GE_PAIRS_FURAN = [(0, 1), (0, 2), (1, 2)]
