"""Planar Kekulé [N]annulene builders for P9 VDE scan."""
from __future__ import annotations

import numpy as np

R_CH = 1.085


def build_annulene(
    n: int,
    r_d: float = 1.350,
    r_s: float = 1.450,
    r_ch: float = R_CH,
) -> tuple[list[str], np.ndarray, list[tuple[int, int]], list[tuple[int, int]]]:
    """Planar even-[N]annulene with alternating double/single bonds.

    Doubles: (0,1),(2,3),...,(n-2,n-1). Carbons 0..n-1, hydrogens n..2n-1.
    """
    if n < 6 or n % 2 != 0:
        raise ValueError(f"N must be even >=6, got {n}")
    coords_c = np.zeros((n, 3))
    x = y = 0.0
    angle = 0.0
    turn = 2.0 * np.pi / n
    for i in range(n):
        coords_c[i] = (x, y, 0.0)
        L = r_d if i % 2 == 0 else r_s
        x += L * np.cos(angle)
        y += L * np.sin(angle)
        angle -= turn
    coords_c -= coords_c.mean(axis=0)
    # H outward in plane
    center = coords_c.mean(axis=0)
    coords_h = np.zeros((n, 3))
    for i in range(n):
        v = coords_c[i] - center
        v = v / (np.linalg.norm(v) + 1e-16)
        coords_h[i] = coords_c[i] + r_ch * v
    coords = np.vstack([coords_c, coords_h])
    symbols = ["C"] * n + ["H"] * n
    doubles = [(i, i + 1) for i in range(0, n, 2)]
    singles = [(i + 1, (i + 2) % n) for i in range(0, n, 2)]
    return symbols, coords, doubles, singles


def adjacent_ge_pairs(n_doubles: int) -> list[tuple[int, int]]:
    """Cyclic adjacent double-bond fragment pairs for GE-m."""
    return [(i, (i + 1) % n_doubles) for i in range(n_doubles)]


def bond_length(coords: np.ndarray, i: int, j: int) -> float:
    return float(np.linalg.norm(coords[i] - coords[j]))


def dmin(coords: np.ndarray) -> float:
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    return float(d.min())


def check_annulene_topology(coords: np.ndarray, n: int, doubles) -> list[str]:
    defects: list[str] = []
    if coords.shape[0] != 2 * n:
        defects.append(f"natm={coords.shape[0]}!={2*n}")
        return defects
    if dmin(coords) < 0.85:
        defects.append(f"dmin={dmin(coords):.3f}")
    for a, b in doubles:
        r = bond_length(coords, a, b)
        if not (1.20 <= r <= 1.55):
            defects.append(f"double {a}-{b}={r:.3f}")
    for i in range(n):
        j = (i + 1) % n
        r = bond_length(coords, i, j)
        if not (1.20 <= r <= 1.70):
            defects.append(f"ring {i}-{j}={r:.3f}")
    if float(np.max(np.abs(coords[:n, 2]))) > 0.05:
        defects.append("non-planar")
    return defects


def is_4n_plus_2(n: int) -> bool:
    return (n - 2) % 4 == 0


def class_label(n: int) -> str:
    return "4n+2" if is_4n_plus_2(n) else "4n"
