"""Polyene / benzene builders for P5 GL–GE-m pilots."""
from __future__ import annotations

import numpy as np

R_CH = 1.085


def _place_chain(cc_bonds: list[float], angles_deg: list[float]) -> np.ndarray:
    """Place n+1 carbons given n bond lengths and n-1 interior turn supplements.

    angles_deg[i] = exterior direction change after atom i+1 (degrees).
    For all-trans polyene use ±60° zig-zag in plane... actually for 120° C-C-C
    use exterior turn of 60°.
    """
    n_bonds = len(cc_bonds)
    coords = np.zeros((n_bonds + 1, 3))
    angle = 0.0
    x = y = 0.0
    for i, length in enumerate(cc_bonds):
        coords[i] = (x, y, 0.0)
        x += length * np.cos(angle)
        y += length * np.sin(angle)
        if i < len(angles_deg):
            angle += np.radians(angles_deg[i])
    coords[n_bonds] = (x, y, 0.0)
    coords -= coords.mean(axis=0)
    return coords


def build_butadiene(
    r_d: float = 1.340, r_s: float = 1.450, r_ch: float = R_CH
) -> tuple[list[str], np.ndarray, list[tuple[int, int]], list[tuple[int, int]]]:
    """trans-1,3-butadiene. Doubles (0,1)(2,3); single (1,2)."""
    cc = _place_chain([r_d, r_s, r_d], [60.0, -60.0])
    # H: terminal CH2 (2H each), internal CH (1H each)
    coords = np.zeros((10, 3))
    coords[:4] = cc
    # rough H placement outward / in-plane
    def add_h(ci, direction, idx):
        n = direction / (np.linalg.norm(direction) + 1e-16)
        coords[idx] = coords[ci] + r_ch * n

    # C0: two H
    t0 = coords[0] - coords[1]
    perp = np.array([-t0[1], t0[0], 0.0])
    add_h(0, t0 + 0.6 * perp, 4)
    add_h(0, t0 - 0.6 * perp, 5)
    # C1: one H
    add_h(1, np.array([-(coords[2] - coords[0])[1], (coords[2] - coords[0])[0], 0.0]), 6)
    # C2: one H
    add_h(2, np.array([-(coords[3] - coords[1])[1], (coords[3] - coords[1])[0], 0.0]), 7)
    # C3: two H
    t3 = coords[3] - coords[2]
    perp = np.array([-t3[1], t3[0], 0.0])
    add_h(3, t3 + 0.6 * perp, 8)
    add_h(3, t3 - 0.6 * perp, 9)
    symbols = ["C"] * 4 + ["H"] * 6
    doubles = [(0, 1), (2, 3)]
    singles = [(1, 2)]
    return symbols, coords, doubles, singles


def build_hexatriene(
    r_d: float = 1.340, r_s: float = 1.450, r_ch: float = R_CH
) -> tuple[list[str], np.ndarray, list[tuple[int, int]], list[tuple[int, int]]]:
    """All-trans hexatriene with approximate Ci symmetry (fixes GE1/GE2 asymmetry).

    Doubles (0,1)(2,3)(4,5); singles (1,2)(3,4). H indices: 6–7 on C0, 8–11 on
    C1–C4, 12–13 on C5.
    """
    cc = _place_chain(
        [r_d, r_s, r_d, r_s, r_d],
        [60.0, -60.0, 60.0, -60.0],
    )
    n_c = 6
    coords = np.zeros((n_c + 8, 3))
    coords[:n_c] = cc

    def add_h(ci: int, direction: np.ndarray, idx: int) -> None:
        n = direction / (np.linalg.norm(direction) + 1e-16)
        coords[idx] = coords[ci] + r_ch * n

    # Terminal CH2: bisect exterior along -bond and ±perp
    t0 = coords[0] - coords[1]
    p0 = np.array([-t0[1], t0[0], 0.0])
    add_h(0, t0 + 0.7 * p0, 6)
    add_h(0, t0 - 0.7 * p0, 7)
    t5 = coords[5] - coords[4]
    p5 = np.array([-t5[1], t5[0], 0.0])
    add_h(5, t5 + 0.7 * p5, 12)
    add_h(5, t5 - 0.7 * p5, 13)

    # Internal vinyl H: alternate sides for all-trans, then enforce Ci on H
    # Place C1,C2 then copy by inversion for C4,C3
    for ci, hi, sign in ((1, 8, +1.0), (2, 9, -1.0)):
        prev, nxt = coords[ci - 1], coords[ci + 1]
        tang = nxt - prev
        perp = np.array([-tang[1], tang[0], 0.0])
        add_h(ci, sign * perp, hi)
    # Ci images: C3↔C2, C4↔C1 (carbon frame already ~Ci about origin)
    coords[10] = -coords[9]  # H on C3 ← −H on C2
    coords[11] = -coords[8]  # H on C4 ← −H on C1
    # Terminal H already ~Ci if chain is Ci; rebuild C5 from C0 images if needed
    coords[12] = -coords[6]
    coords[13] = -coords[7]

    symbols = ["C"] * 6 + ["H"] * 8
    doubles = [(0, 1), (2, 3), (4, 5)]
    singles = [(1, 2), (3, 4)]
    return symbols, coords, doubles, singles


def build_benzene_kekule(
    r_d: float = 1.340, r_s: float = 1.450, r_ch: float = R_CH
) -> tuple[list[str], np.ndarray, list[tuple[int, int]], list[tuple[int, int]]]:
    """Planar Kekulé benzene: alternating r_d / r_s (exterior turns 60°).

    Doubles (0,1)(2,3)(4,5); singles (1,2)(3,4)(5,0).
    """
    from src.p4_benzene.geometry import build_benzene_d3h

    symbols, coords = build_benzene_d3h(r_d, r_s, r_ch)
    doubles = [(0, 1), (2, 3), (4, 5)]
    singles = [(1, 2), (3, 4), (5, 0)]
    return symbols, coords, doubles, singles


def set_bond_length(coords: np.ndarray, i: int, j: int, length: float, move: str = "j") -> np.ndarray:
    """Scale bond i–j to `length`, moving atom j (and optionally a fragment)."""
    c = coords.copy()
    v = c[j] - c[i]
    n = np.linalg.norm(v)
    if n < 1e-12:
        return c
    scale = length / n
    if move == "j":
        c[j] = c[i] + v * scale
    else:
        c[i] = c[j] - v * scale
    c[:, 2] = 0.0
    return c


def shift_fragment(coords: np.ndarray, atoms: list[int], shift: np.ndarray) -> np.ndarray:
    c = coords.copy()
    for a in atoms:
        c[a] = c[a] + shift
    c[:, 2] = 0.0
    return c


def set_single_bond_keep_sides(
    coords: np.ndarray,
    i: int,
    j: int,
    length: float,
    right_atoms: list[int],
) -> np.ndarray:
    """Set |rj−ri|=length by translating the right-hand fragment."""
    c = coords.copy()
    v = c[j] - c[i]
    n = np.linalg.norm(v) + 1e-16
    target = v / n * length
    shift = target - v
    for a in right_atoms:
        c[a] = c[a] + shift
    c[:, 2] = 0.0
    return c


def check_polyene_topology(
    symbols: list[str],
    coords: np.ndarray,
    n_c: int,
    doubles: list[tuple[int, int]],
) -> list[str]:
    problems: list[str] = []
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    if float(d.min()) < 0.85:
        problems.append(f"dmin={float(d.min()):.3f}")
    for a, b in doubles:
        rij = float(np.linalg.norm(coords[a] - coords[b]))
        if not (1.20 < rij < 1.55):
            problems.append(f"double {a}-{b} r={rij:.3f}")
    # chain connectivity
    for i in range(n_c - 1):
        rij = float(np.linalg.norm(coords[i] - coords[i + 1]))
        if not (1.20 < rij < 1.70):
            problems.append(f"chain {i}-{i+1} r={rij:.3f}")
    if float(np.max(np.abs(coords[:, 2]))) > 0.05:
        problems.append("non-planar")
    return problems


def to_xyz(symbols: list[str], coords: np.ndarray, comment: str = "") -> str:
    lines = [f"{len(symbols)}", comment]
    for s, c in zip(symbols, coords):
        lines.append(f"{s} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}")
    return "\n".join(lines) + "\n"
