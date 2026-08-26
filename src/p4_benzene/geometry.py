"""Benzene geometry helpers for P4 BLA / Ee–EN scans."""
from __future__ import annotations

import numpy as np

# Ring carbons 0–5, hydrogens 6–11 (H_i on C_i)
N_C = 6
N_H = 6
R_CH = 1.085


def build_benzene_d3h(r_a: float, r_b: float, r_ch: float = R_CH) -> tuple[list[str], np.ndarray]:
    """Planar Kekulé benzene with alternating CC lengths r_a, r_b (Å).

    Construction: successive bonds a,b,a,b,a,b with exterior turn 60° (interior 120°).
    Closes for any a,b (book prerequisite dra = −drb when scanning toward equalization).
    """
    coords = np.zeros((N_C + N_H, 3), dtype=float)
    angle = 0.0
    x = y = 0.0
    for i in range(N_C):
        coords[i] = (x, y, 0.0)
        length = r_a if (i % 2 == 0) else r_b
        x += length * np.cos(angle)
        y += length * np.sin(angle)
        angle += np.pi / 3.0
    # translate to centroid
    c = coords[:N_C].mean(axis=0)
    coords[:N_C] -= c
    for i in range(N_C):
        v = coords[i].copy()
        n = np.linalg.norm(v)
        if n < 1e-8:
            v = np.array([1.0, 0.0, 0.0])
            n = 1.0
        coords[N_C + i] = coords[i] + r_ch * v / n
    symbols = ["C"] * N_C + ["H"] * N_H
    return symbols, coords


def build_benzene_equal(r_cc: float, r_ch: float = R_CH) -> tuple[list[str], np.ndarray]:
    return build_benzene_d3h(r_cc, r_cc, r_ch)


def bla_delta(coords: np.ndarray) -> float:
    """δ = mean(odd bonds) − mean(even bonds); ≈ r_long − r_short for D3h."""
    bonds = []
    for i in range(N_C):
        j = (i + 1) % N_C
        bonds.append(float(np.linalg.norm(coords[i] - coords[j])))
    even = np.mean(bonds[0::2])
    odd = np.mean(bonds[1::2])
    return float(odd - even) if odd >= even else float(even - odd)


def cc_bonds(coords: np.ndarray) -> list[float]:
    return [
        float(np.linalg.norm(coords[i] - coords[(i + 1) % N_C])) for i in range(N_C)
    ]


def check_topology(symbols: list[str], coords: np.ndarray) -> list[str]:
    """Return defects; empty ⇒ C6H6 six-membered ring with expected bonding."""
    problems: list[str] = []
    if len(symbols) != 12 or symbols.count("C") != 6 or symbols.count("H") != 6:
        problems.append(f"formula: expected C6H6 got {''.join(symbols)}")
        return problems
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    dmin = float(d.min())
    if dmin < 0.85:
        problems.append(f"dmin={dmin:.3f} Å < 0.85")

    # each C bonded to exactly 2 C (1.20–1.70) and 1 H (0.90–1.30)
    for i in range(N_C):
        c_neighbors = [j for j in range(N_C) if 1.15 < d[i, j] < 1.75]
        h_neighbors = [j for j in range(N_C, N_C + N_H) if 0.90 < d[i, j] < 1.30]
        if len(c_neighbors) != 2:
            problems.append(f"C{i} has {len(c_neighbors)} C neighbors (want 2)")
        if len(h_neighbors) != 1:
            problems.append(f"C{i} has {len(h_neighbors)} H neighbors (want 1)")

    # ring size: follow C–C around once → 6
    try:
        start = 0
        path = [start]
        prev, cur = -1, start
        for _ in range(6):
            nbrs = [j for j in range(N_C) if 1.15 < d[cur, j] < 1.75 and j != prev]
            if not nbrs:
                problems.append("ring walk broken")
                break
            nxt = nbrs[0] if nbrs[0] != path[0] or len(path) == 1 else (
                nbrs[1] if len(nbrs) > 1 else nbrs[0]
            )
            if len(path) > 1 and nxt == path[0]:
                break
            path.append(nxt)
            prev, cur = cur, nxt
        if len(set(path)) != 6:
            problems.append(f"ring carbons visited {len(set(path))} ≠ 6")
    except Exception as exc:  # noqa: BLE001
        problems.append(f"ring walk error: {exc}")

    # planarity
    z = coords[:, 2]
    if float(np.max(np.abs(z - z.mean()))) > 0.05:
        problems.append("non-planar (|z| > 0.05 Å)")

    return problems


def to_xyz(symbols: list[str], coords: np.ndarray) -> str:
    lines = [f"{len(symbols)}", "benzene"]
    for s, c in zip(symbols, coords):
        lines.append(f"{s} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}")
    return "\n".join(lines) + "\n"
