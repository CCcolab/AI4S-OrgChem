"""NBA (Ph–N=CH–Ph) geometry and twist dihedral helpers."""
from __future__ import annotations

import numpy as np

IDX_N = 0
IDX_C_IMINE = 1
IDX_CIPSO_N = 3
IDX_CORTHO_N = 4
NPH_FRAGMENT = list(range(4, 14))


def _phenyl(origin: np.ndarray, x_axis: np.ndarray, y_axis: np.ndarray) -> np.ndarray:
    """C6H5 as a regular hexagon, C_ipso at `origin`.

    `x_axis` points from C_ipso towards its substituent, so the ring is built
    on the opposite side. A regular hexagon of side r_cc has circumradius
    r_cc, so its centre sits one bond length behind C_ipso. Returned order is
    ipso, ortho, meta, para, meta', ortho' followed by the five H.
    """
    x = x_axis / (np.linalg.norm(x_axis) + 1e-16)
    y = y_axis - np.dot(y_axis, x) * x
    y = y / (np.linalg.norm(y) + 1e-16)
    r_cc, r_ch = 1.390, 1.080
    centre = origin - r_cc * x
    ring = [
        centre + r_cc * (np.cos(np.radians(a)) * x + np.sin(np.radians(a)) * y)
        for a in (0, 60, 120, 180, 240, 300)
    ]
    hydrogens = [
        centre + (r_cc + r_ch) * (np.cos(np.radians(a)) * x + np.sin(np.radians(a)) * y)
        for a in (60, 120, 180, 240, 300)
    ]
    return np.asarray(ring + hydrogens, dtype=float)


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


# Index layout: 0 N, 1 C_imine, 2 H_imine,
# 3-8 N-phenyl C (ipso, ortho, meta, para, meta', ortho'), 9-13 its H,
# 14-19 C-phenyl C (same order), 20-24 its H.
RING_N = (3, 4, 5, 6, 7, 8)
RING_C = (14, 15, 16, 17, 18, 19)
EXPECTED_BONDS = {
    (0, 1),  # N=C imine
    (0, 3),  # N-C_ipso(N-phenyl)
    (1, 2),  # C_imine-H
    (1, 14),  # C_imine-C_ipso(C-phenyl)
}


def check_topology(symbols: list[str], coords: np.ndarray, tol: float = 0.25) -> list[str]:
    """Return a list of structural defects; empty means the graph is NBA.

    Guards against a builder or an optimizer silently producing a different
    molecule (e.g. a ring of the wrong size) whose energies would look
    perfectly plausible.
    """
    problems: list[str] = []
    n = len(symbols)
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)

    def bonded(i: int, j: int) -> bool:
        cut = 1.25 if "H" in (symbols[i], symbols[j]) else 1.75
        return bool(d[i, j] < cut)

    if float(d.min()) < 0.85:
        i, j = np.unravel_index(int(np.argmin(d)), d.shape)
        problems.append(f"atoms {i},{j} only {d[i, j]:.3f} A apart")

    for ring, name in ((RING_N, "N-phenyl"), (RING_C, "C-phenyl")):
        for pos in range(6):
            a, b = ring[pos], ring[(pos + 1) % 6]
            if not bonded(a, b):
                problems.append(f"{name}: ring bond {a}-{b} missing (d={d[a, b]:.3f})")
        ipso, ortho, meta, para = ring[0], ring[1], ring[2], ring[3]
        for label, idx, ref in (
            ("ortho", ortho, 1.39),
            ("meta", meta, 2.41),
            ("para", para, 2.78),
        ):
            if abs(d[ipso, idx] - ref) > tol + (0.25 if label == "para" else 0.0):
                problems.append(
                    f"{name}: d(ipso,{label})={d[ipso, idx]:.3f} A, expected ~{ref}"
                )
        # A benzene carbon has three neighbours; the ipso one carries no H.
        for idx in ring:
            deg = sum(1 for k in range(n) if k != idx and bonded(idx, k))
            if deg != 3:
                problems.append(f"{name}: C{idx} has {deg} neighbours, expected 3")

    for i, j in sorted(EXPECTED_BONDS):
        if not bonded(i, j):
            problems.append(f"expected bond {i}-{j} missing (d={d[i, j]:.3f})")

    for idx, sym in enumerate(symbols):
        if sym != "H":
            continue
        deg = sum(1 for k in range(n) if k != idx and bonded(idx, k))
        if deg != 1:
            problems.append(f"H{idx} has {deg} neighbours, expected 1")
    return problems


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
