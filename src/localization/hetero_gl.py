"""
Heterocycle GL / GE-m — extend 2007 with extra π-only fragments (O/N lone pair).

Fragment order: one per double bond, then each extra_pi_atom_group.
Furan GE-m (public Ch6 Fig 6-8): merge (C2=C3,C4=C5), (O,C2=C3), (O,C4=C5).
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from pyscf import dft, gto, scf

from src.localization.gl_2007 import (
    assign_pi_to_double_bonds,
    identify_pi_ao_indices,
    merge_fragments,
    zero_inter_fragment_pi,
)
from src.localization.plg import pi_aos_on_atoms


def build_fragments(
    mol: gto.Mole,
    double_bond_atoms: Sequence[tuple[int, int]],
    extra_pi_atom_groups: Sequence[Sequence[int]] | None = None,
) -> list[list[int]]:
    pi = identify_pi_ao_indices(mol)
    frags = assign_pi_to_double_bonds(mol, pi, double_bond_atoms)
    if extra_pi_atom_groups:
        for atoms in extra_pi_atom_groups:
            frag = pi_aos_on_atoms(mol, pi, atoms)
            if not frag:
                raise RuntimeError(
                    f"empty hetero π fragment for atoms={atoms} pi={pi} "
                    f"labels={[mol.ao_labels(fmt=True)[i] for i in pi]}"
                )
            frags.append(frag)
    if any(len(f) == 0 for f in frags):
        raise RuntimeError(
            f"π assignment failed: pi={pi} frags={frags} "
            f"doubles={double_bond_atoms} extra={extra_pi_atom_groups}"
        )
    return frags


def make_hetero_localized_mf(
    mol: gto.Mole,
    method: str,
    double_bond_atoms: Sequence[tuple[int, int]],
    *,
    extra_pi_atom_groups: Sequence[Sequence[int]] | None = None,
    allow_pair: tuple[int, int] | None = None,
    zero_overlap: bool = True,
    zero_exchange: bool = False,
):
    if method.upper() == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
    else:
        mf = scf.RHF(mol)

    base = build_fragments(mol, double_bond_atoms, extra_pi_atom_groups)
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
    mf.level_shift = 0.5
    mf.max_cycle = 200
    mf.conv_tol = 1e-8
    mf._hetero_frags_base = base
    mf._hetero_frags = frags
    mf._hetero_allow_pair = allow_pair
    return mf


def make_hetero_gl_mf(mol, method, doubles, extra_pi=None, **kw):
    return make_hetero_localized_mf(
        mol,
        method,
        doubles,
        extra_pi_atom_groups=extra_pi,
        allow_pair=None,
        **kw,
    )


def make_hetero_ge_mf(mol, method, doubles, extra_pi, allow_pair, **kw):
    return make_hetero_localized_mf(
        mol,
        method,
        doubles,
        extra_pi_atom_groups=extra_pi,
        allow_pair=allow_pair,
        **kw,
    ).newton()


def run_hetero_scf(mol, method, doubles, extra_pi, *, allow_pair=None, dm0=None):
    """GL/GE with DIIS then Newton fallback."""
    if allow_pair is None:
        mf = make_hetero_gl_mf(mol, method, doubles, extra_pi)
    else:
        mf = make_hetero_ge_mf(mol, method, doubles, extra_pi, allow_pair)
    e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    if mf.converged:
        return e, mf
    if allow_pair is None:
        mf2 = make_hetero_gl_mf(mol, method, doubles, extra_pi)
    else:
        mf2 = make_hetero_ge_mf(mol, method, doubles, extra_pi, allow_pair)
    mf2.level_shift = 0.3
    mf2.damp = 0.0
    mf2 = mf2.newton()
    seed = dm0 if dm0 is not None else mf.make_rdm1()
    e2 = float(mf2.kernel(dm0=seed))
    return e2, mf2
