"""Minimal PySCF molecule / SCF helpers (no third-party extras)."""
from __future__ import annotations

from typing import Any

from pyscf import dft, gto, scf


def build_mol(atom: str, basis: str = "6-31g*", charge: int = 0, spin: int = 0, unit: str = "Angstrom") -> gto.Mole:
    mol = gto.Mole()
    mol.atom = atom
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.unit = unit
    mol.verbose = 0
    mol.build()
    return mol


def energy_rks(mol: gto.Mole, xc: str = "B3LYP") -> tuple[float, Any]:
    mf = dft.RKS(mol)
    mf.xc = xc
    e = mf.kernel()
    return float(e), mf


def energy_rhf(mol: gto.Mole) -> tuple[float, Any]:
    mf = scf.RHF(mol)
    e = mf.kernel()
    return float(e), mf


def energy_mp2(mol: gto.Mole) -> tuple[float, float]:
    """Return (E_HF, E_MP2_total)."""
    from pyscf import mp

    mf = scf.RHF(mol).run()
    pt = mp.MP2(mf).run()
    return float(mf.e_tot), float(pt.e_tot)
