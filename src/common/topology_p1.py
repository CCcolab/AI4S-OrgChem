"""Shared helpers for P1-style small-molecule topology / formula checks."""
from __future__ import annotations

from collections import Counter

import numpy as np

# Expected stoichiometry for P1 species
P1_FORMULAS: dict[str, Counter] = {
    "h2": Counter({"H": 2}),
    "butadiene": Counter({"C": 4, "H": 6}),
    "1_butene": Counter({"C": 4, "H": 8}),
    "trans_2_butene": Counter({"C": 4, "H": 8}),
    "n_butane": Counter({"C": 4, "H": 10}),
}


def parse_xyz_block(xyz: str) -> tuple[list[str], np.ndarray]:
    lines = [ln for ln in xyz.strip().splitlines() if ln.strip()]
    # bare "Sym x y z" block (no natoms header)
    if len(lines[0].split()) >= 4 and lines[0].split()[0].isalpha():
        body = lines
    else:
        n = int(lines[0].split()[0])
        body = lines[2 : 2 + n]
    syms, rows = [], []
    for ln in body:
        p = ln.split()
        syms.append(p[0])
        rows.append([float(x) for x in p[1:4]])
    return syms, np.asarray(rows, dtype=float)


def dmat(coords: np.ndarray) -> np.ndarray:
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    return d


def check_p1_molecule(name: str, xyz: str) -> list[str]:
    """Return defect list; empty means topology/formula look like the named species."""
    issues: list[str] = []
    if name not in P1_FORMULAS:
        return [f"unknown species {name}"]
    syms, coords = parse_xyz_block(xyz)
    form = Counter(syms)
    if form != P1_FORMULAS[name]:
        issues.append(f"formula {dict(form)} != {dict(P1_FORMULAS[name])}")
    d = dmat(coords)
    # clash excluding normal H–H / C–H
    for i in range(len(syms)):
        for j in range(i + 1, len(syms)):
            rij = float(d[i, j])
            pair = {syms[i], syms[j]}
            if pair == {"H"} and rij > 0.60:
                continue
            if "H" in pair and rij > 0.90:
                continue
            if rij < 0.85:
                issues.append(f"clash {i}-{j} d={rij:.3f}")
    if name == "h2":
        return issues
    cs = [i for i, s in enumerate(syms) if s == "C"]
    cc = [
        (a, b, float(d[a, b]))
        for ia, a in enumerate(cs)
        for b in cs[ia + 1 :]
        if d[a, b] < 1.70
    ]
    if len(cc) != 3:
        issues.append(f"C–C count {len(cc)} != 3")
    else:
        lens = sorted(r for _, _, r in cc)
        if name == "butadiene" and not (
            lens[0] < 1.42 and lens[1] < 1.42 and lens[2] > 1.44
        ):
            issues.append(f"butadiene not D-S-D: {lens}")
        if name == "1_butene" and not (
            lens[0] < 1.40 and lens[1] > 1.45 and lens[2] > 1.45
        ):
            issues.append(f"1-butene pattern odd: {lens}")
        if name == "n_butane" and any(r < 1.45 for r in lens):
            issues.append(f"butane short CC: {lens}")
        if name == "trans_2_butene":
            a, b, _ = min(cc, key=lambda t: t[2])
            ends = [c for c in cs if c not in (a, b)]
            if len(ends) == 2:
                th = _dihedral(coords[ends[0]], coords[a], coords[b], coords[ends[1]])
                if abs(abs(th) - 180.0) > 30.0:
                    issues.append(f"2-butene dihedral {th:.1f} not trans")
    return issues


def _dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1u = b1 / (np.linalg.norm(b1) + 1e-16)
    v = b0 - np.dot(b0, b1u) * b1u
    w = b2 - np.dot(b2, b1u) * b1u
    x = np.dot(v, w)
    y = np.dot(np.cross(b1u, v), w)
    return float(np.degrees(np.arctan2(y, x)))
