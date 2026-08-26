"""
Localized GL / GE-m SCF — 2007 public matrix-deletion definition (Ch6).

GL: zero all inter-double-bond π–π Fock and overlap blocks.
GE-m: merge one allowed conjugated pair of double bonds into a super-fragment;
      zero π–π blocks only between distinct remaining fragments.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from pyscf import dft, gto, scf


def identify_pi_ao_indices(mol: gto.Mole, thresh: float = 1e-9) -> list[int]:
    """π-AOs as in Ch6: |hcore(0, i)| < thresh and decoupled from all σ-AOs."""
    h = scf.hf.get_hcore(mol)
    n = mol.nao
    candidates = [i for i in range(n) if abs(h[0, i]) < thresh]
    sigma = [i for i in range(n) if i not in candidates]
    pi = []
    for kp in candidates:
        if all(abs(h[kp, ks]) < thresh for ks in sigma):
            pi.append(kp)
    return pi


def assign_pi_to_double_bonds(
    mol: gto.Mole,
    pi_aos: Sequence[int],
    double_bond_atoms: Sequence[tuple[int, int]],
) -> list[list[int]]:
    """Assign each π-AO to the double bond that owns its atom center."""
    frags: list[list[int]] = [[] for _ in double_bond_atoms]
    atom_to_frag: dict[int, int] = {}
    for fi, (a, b) in enumerate(double_bond_atoms):
        atom_to_frag[a] = fi
        atom_to_frag[b] = fi
    aoslices = mol.aoslice_by_atom()
    for mu in pi_aos:
        placed = False
        for atom, (_p0, _p1, ao0, ao1) in enumerate(aoslices):
            if ao0 <= mu < ao1 and atom in atom_to_frag:
                frags[atom_to_frag[atom]].append(mu)
                placed = True
                break
        if not placed:
            lab = mol.ao_labels(fmt=True)[mu]
            atom = int(lab.split()[0])
            if atom in atom_to_frag:
                frags[atom_to_frag[atom]].append(mu)
    return frags


def merge_fragments(
    frags: Sequence[Sequence[int]], allow_pair: tuple[int, int] | None
) -> list[list[int]]:
    """For GE-m: merge fragment indices allow_pair into one super-fragment."""
    if allow_pair is None:
        return [list(f) for f in frags]
    i, j = allow_pair
    if i > j:
        i, j = j, i
    out: list[list[int]] = []
    merged = list(frags[i]) + list(frags[j])
    for k, f in enumerate(frags):
        if k == i:
            out.append(merged)
        elif k == j:
            continue
        else:
            out.append(list(f))
    return out


def zero_inter_fragment_pi(
    mat: np.ndarray, fragments: Sequence[Sequence[int]]
) -> np.ndarray:
    out = np.array(mat, dtype=float, copy=True)
    for i, fi in enumerate(fragments):
        for fj in fragments[i + 1 :]:
            for a in fi:
                for b in fj:
                    out[a, b] = 0.0
                    out[b, a] = 0.0
    return out


def make_localized_mf(
    mol: gto.Mole,
    method: str,
    double_bond_atoms: Sequence[tuple[int, int]],
    *,
    allow_pair: tuple[int, int] | None = None,
    zero_overlap: bool = True,
    zero_exchange: bool = False,
):
    """Mean-field with 2007 GL/GE-m; optional 2011-lite exchange-block deletion.

    zero_exchange: zero inter-fragment π–π blocks of the exchange matrix K
    (public 2011 idea: delete inter-double-bond exchange). Independent of
    copying any original code.
    """
    if method.upper() == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
    else:
        mf = scf.RHF(mol)

    pi = identify_pi_ao_indices(mol)
    base = assign_pi_to_double_bonds(mol, pi, double_bond_atoms)
    if any(len(f) == 0 for f in base):
        raise RuntimeError(
            f"π-AO assignment failed: pi={pi} frags={base} "
            f"labels={[mol.ao_labels(fmt=True)[i] for i in pi]}"
        )
    frags = merge_fragments(base, allow_pair)

    _get_ovlp = mf.get_ovlp
    _get_fock = mf.get_fock
    _get_jk = mf.get_jk

    def get_ovlp(mol_=None):
        s = np.asarray(_get_ovlp(mol_), dtype=float)
        if zero_overlap:
            return zero_inter_fragment_pi(s, frags)
        return s

    def get_fock(h1e=None, s1e=None, vhf=None, dm=None, cycle=-1, mf_prev=None, **kw):
        fock = _get_fock(h1e, s1e, vhf, dm, cycle, mf_prev, **kw)
        return zero_inter_fragment_pi(np.asarray(fock, dtype=float), frags)

    def get_jk(mol_=None, dm=None, hermi=1, with_j=True, with_k=True, omega=None, **kw):
        vj, vk = _get_jk(mol_, dm, hermi=hermi, with_j=with_j, with_k=with_k, omega=omega, **kw)
        if zero_exchange and with_k and vk is not None:
            vk = zero_inter_fragment_pi(np.asarray(vk, dtype=float), frags)
        return vj, vk

    mf.get_ovlp = get_ovlp
    mf.get_fock = get_fock
    if zero_exchange:
        mf.get_jk = get_jk
    mf.level_shift = 0.5 if zero_overlap else 0.2
    mf.max_cycle = 200
    mf.conv_tol = 1e-8
    mf._gl_pi_aos = pi
    mf._gl_fragments_base = base
    mf._gl_fragments = frags
    mf._gl_allow_pair = allow_pair
    mf._gl_zero_overlap = zero_overlap
    mf._gl_zero_exchange = zero_exchange
    return mf


# Back-compat alias used by P1
def make_gl_mf(
    mol: gto.Mole,
    method: str,
    double_bond_atoms: Sequence[tuple[int, int]],
    *,
    zero_overlap: bool = True,
):
    return make_localized_mf(
        mol, method, double_bond_atoms, allow_pair=None, zero_overlap=zero_overlap
    )
