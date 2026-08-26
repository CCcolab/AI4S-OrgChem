"""
Particular Localized Geometry (PLG) SCF — public Ch10 definition.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from pyscf import dft, gto, scf

from src.localization.gl_2007 import identify_pi_ao_indices, zero_inter_fragment_pi


def pi_aos_on_atoms(mol: gto.Mole, pi_aos: Sequence[int], atoms: Sequence[int]) -> list[int]:
    atom_set = set(int(a) for a in atoms)
    aoslices = mol.aoslice_by_atom()
    out: list[int] = []
    for mu in pi_aos:
        for atom, (_p0, _p1, ao0, ao1) in enumerate(aoslices):
            if ao0 <= mu < ao1 and atom in atom_set:
                out.append(mu)
                break
    return out


def make_plg_mf(
    mol: gto.Mole,
    method: str,
    atoms_A: Sequence[int],
    atoms_B: Sequence[int],
    *,
    zero_overlap: bool = True,
    zero_exchange: bool = True,
):
    if method.upper() == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
    else:
        mf = scf.RHF(mol)

    pi = identify_pi_ao_indices(mol)
    frag_a = pi_aos_on_atoms(mol, pi, atoms_A)
    frag_b = pi_aos_on_atoms(mol, pi, atoms_B)
    if not frag_a or not frag_b:
        raise RuntimeError(
            f"PLG π assignment failed: pi={pi} A={frag_a} B={frag_b} "
            f"labels={[mol.ao_labels(fmt=True)[i] for i in pi]}"
        )
    frags = [frag_a, frag_b]

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
    mf.damp = 0.3
    mf.diis_space = 12
    mf.max_cycle = 80
    mf.conv_tol = 1e-7
    mf._plg_pi = pi
    mf._plg_frags = frags
    mf._plg_atoms_A = list(atoms_A)
    mf._plg_atoms_B = list(atoms_B)
    return mf


def run_plg_scf(mol: gto.Mole, method: str, atoms_A, atoms_B, dm0=None):
    """DIIS warm-up then Newton (required for reliable PLG convergence at 6-31G*)."""
    mf = make_plg_mf(mol, method, atoms_A, atoms_B)
    e = float(mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel())
    if mf.converged:
        return e, mf
    mf2 = make_plg_mf(mol, method, atoms_A, atoms_B)
    mf2.level_shift = 0.2
    mf2.damp = 0.0
    mf2 = mf2.newton()
    seed = dm0 if dm0 is not None else mf.make_rdm1()
    e2 = float(mf2.kernel(dm0=seed))
    return e2, mf2
