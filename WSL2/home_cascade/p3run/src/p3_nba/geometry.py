"""NBA (Ph–N=CH–Ph) geometry and twist dihedral helpers."""
from __future__ import annotations

import numpy as np

IDX_N = 0
IDX_C_IMINE = 1
IDX_CIPSO_N = 3
IDX_CORTHO_N = 4
NPH_FRAGMENT = list(range(4, 14))


def _phenyl(origin: np.ndarray, x_axis: np.ndarray, y_axis: np.ndarray) -> np.ndarray:
    """C6H5 with C_ipso at origin; ring in plane of x_axis/y_axis."""
    x = x_axis / (np.linalg.norm(x_axis) + 1e-16)
    y = y_axis - np.dot(y_axis, x) * x
    y = y / (np.linalg.norm(y) + 1e-16)
    r_cc, r_ch = 1.390, 1.080
    coords = [origin.copy()]
    for a in (60, 120, 180, 240, 300):
        th = np.radians(a)
        coords.append(origin + r_cc * (np.cos(th) * x + np.sin(th) * y))
    for a in (60, 120, 180, 240, 300):
        th = np.radians(a)
        coords.append(origin + (r_cc + r_ch) * (np.cos(th) * x + np.sin(th) * y))
    return np.asarray(coords, dtype=float)


def build_nba(twist_deg: float = 40.0) -> tuple[list[str], np.ndarray]:
    """Ph–N=CH–Ph with bent ∠C–N–C (~120°) so twist dihedral is well-defined."""
    n = np.zeros(3)
    c_im = np.array([1.280, 0.0, 0.0])
    h_im = np.array([1.820, 0.950, 0.0])
    c_ipso_n = np.array([-0.700, 1.212, 0.0])
    x_nph = n - c_ipso_n
    y_nph = np.array([0.0, 0.0, 1.0])
    xyz_nph = _phenyl(c_ipso_n, x_nph, y_nph)
    c_ipso_c = np.array([2.050, -1.200, 0.0])
    x_cph = c_im - c_ipso_c
    y_cph = np.array([0.0, 0.0, 1.0])
    xyz_cph = _phenyl(c_ipso_c, x_cph, y_cph)
    symbols = ["N", "C", "H"] + ["C"] * 6 + ["H"] * 5 + ["C"] * 6 + ["H"] * 5
    coords = np.vstack([n, c_im, h_im, xyz_nph, xyz_cph])
    return symbols, set_twist_deg(coords, twist_deg)


def format_xyz(symbols: list[str], coords: np.ndarray) -> str:
    return "\n".join(
        f"{s} {coords[i, 0]:.8f} {coords[i, 1]:.8f} {coords[i, 2]:.8f}"
        for i, s in enumerate(symbols)
    )


def dihedral_deg(coords: np.ndarray, i: int, j: int, k: int, l: int) -> float:
    b0 = coords[i] - coords[j]
    b1 = coords[k] - coords[j]
    b2 = coords[l] - coords[k]
    b1u = b1 / (np.linalg.norm(b1) + 1e-16)
    v = b0 - np.dot(b0, b1u) * b1u
    w = b2 - np.dot(b2, b1u) * b1u
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1u, v), w), np.dot(v, w))))


def dihedral_grad_deg(coords: np.ndarray, i: int, j: int, k: int, l: int) -> np.ndarray:
    g = np.zeros_like(coords)
    eps = 1.0e-4
    for a in (i, j, k, l):
        for t in range(3):
            c = coords.copy()
            c[a, t] += eps
            f1 = dihedral_deg(c, i, j, k, l)
            c[a, t] -= 2 * eps
            f0 = dihedral_deg(c, i, j, k, l)
            d = f1 - f0
            if d > 180:
                d -= 360
            if d < -180:
                d += 360
            g[a, t] = d / (2 * eps)
    return g


def rotate_fragment(
    coords: np.ndarray, axis_a: int, axis_b: int, frag: list[int], angle_deg: float
) -> np.ndarray:
    out = coords.copy()
    origin = out[axis_a]
    axis = out[axis_b] - origin
    axis = axis / (np.linalg.norm(axis) + 1e-16)
    th = np.radians(angle_deg)
    c, s = np.cos(th), np.sin(th)
    K = np.array(
        [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
    )
    R = np.eye(3) + s * K + (1 - c) * (K @ K)
    for idx in frag:
        out[idx] = origin + R @ (out[idx] - origin)
    return out


def set_twist_deg(coords: np.ndarray, target_deg: float) -> np.ndarray:
    cur = dihedral_deg(coords, IDX_CORTHO_N, IDX_CIPSO_N, IDX_N, IDX_C_IMINE)
    delta = target_deg - cur
    while delta > 180:
        delta -= 360
    while delta < -180:
        delta += 360
    return rotate_fragment(coords, IDX_CIPSO_N, IDX_N, NPH_FRAGMENT, -delta)
