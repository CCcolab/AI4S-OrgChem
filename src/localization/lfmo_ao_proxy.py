"""
LFMO-lite: public Fig. 5-15 energy effects via fragment-AO densities.

Not Kost LFMO. States are constructed from the ground-state Fock by
absolutely localized (block) subspace diagonalization.

  EV   = E(FUD) − E(DSI)     π–π (delocalized π vs fragment π densities)
  Enσσ = 2∑ D_μν F_μν        nonbonded σ–σ Fock coupling (link atoms excluded)
  Eπσ  = E(G) − E(FUD)       π–σ; at planar NBA this is ~0 (Table 5-16)
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from pyscf import dft, gto, scf

FRAG_A = list(range(3, 14))
FRAG_BC = [0, 1, 2] + list(range(14, 25))
RING_N = (3, 4, 5, 6, 7, 8)
RING_C = (14, 15, 16, 17, 18, 19)
# C13H11N = 96e → 48 occupied MOs; 14 π e → 7 occ π; 82 σ e → 41 occ σ
NOCC_PI = 7
NOCC_SIGMA = 41
NOCC_PI_A = 3  # N-phenyl 6π
NOCC_PI_BC = 4  # imine 2π + C-phenyl 6π
# σ ALMO split (closed-shell): A takes linking C–N pair (18) + BC (23) = 41
NOCC_SIG_A = 18
NOCC_SIG_BC = 23


def _plane_normal(coords: np.ndarray, i: int, j: int, k: int) -> np.ndarray:
    n = np.cross(coords[j] - coords[i], coords[k] - coords[i])
    return n / (np.linalg.norm(n) + 1e-16)


def _p_kind(label: str) -> str | None:
    lab = label.lower().replace("_", "")
    # require px/py/pz token; do not match trailing x in '2px' via endswith only after 'p'
    for kind in ("px", "py", "pz"):
        if kind in lab:
            return kind[-1]
    return None


def identify_pi_sigma_aos(
    mol: gto.Mole, coords: np.ndarray
) -> tuple[list[int], list[int], dict[int, str]]:
    n_n = _plane_normal(coords, RING_N[0], RING_N[1], RING_N[2])
    n_c = _plane_normal(coords, RING_C[0], RING_C[1], RING_C[2])
    n_im = _plane_normal(coords, 3, 0, 1)
    atom_normal = {a: n_n for a in RING_N}
    atom_normal.update({a: n_c for a in RING_C})
    atom_normal[0] = n_im
    atom_normal[1] = n_im

    labels = mol.ao_labels(fmt=True)
    aoslices = mol.aoslice_by_atom()
    axes = {
        "x": np.array([1.0, 0.0, 0.0]),
        "y": np.array([0.0, 1.0, 0.0]),
        "z": np.array([0.0, 0.0, 1.0]),
    }
    pi: list[int] = []
    detail: dict[int, str] = {}
    for atom, normal in atom_normal.items():
        ao0, ao1 = aoslices[atom][2], aoslices[atom][3]
        best_i, best_abs, best_lab = None, -1.0, ""
        for mu in range(ao0, ao1):
            kind = _p_kind(labels[mu])
            if kind is None:
                continue
            score = abs(float(np.dot(axes[kind], normal)))
            if score > best_abs:
                best_abs, best_i, best_lab = score, mu, labels[mu]
        if best_i is not None and best_abs >= 0.40:
            pi.append(int(best_i))
            detail[int(best_i)] = f"{best_lab}|align={best_abs:.3f}"
    pi = sorted(set(pi))
    sigma = [i for i in range(mol.nao) if i not in set(pi)]
    return pi, sigma, detail


def _atom_aos(mol: gto.Mole, atoms: Sequence[int]) -> list[int]:
    aoslices = mol.aoslice_by_atom()
    out: list[int] = []
    for a in atoms:
        out.extend(range(aoslices[a][2], aoslices[a][3]))
    return out


def _eigh_sub(F: np.ndarray, S: np.ndarray, nocc: int) -> tuple[np.ndarray, np.ndarray]:
    """Occupied MO coefficients in the given AO subspace (canonical orthog.)."""
    se, su = np.linalg.eigh(S)
    keep = se > 1e-8
    x = su[:, keep] / np.sqrt(se[keep])
    fx = x.T @ F @ x
    e, u = np.linalg.eigh(fx)
    c = x @ u
    if nocc > c.shape[1]:
        raise RuntimeError(f"nocc={nocc} > nbf={c.shape[1]}")
    return c[:, :nocc], e


def _embed_mos(nao: int, aos: Sequence[int], c_sub: np.ndarray) -> np.ndarray:
    c = np.zeros((nao, c_sub.shape[1]))
    c[list(aos), :] = c_sub
    return c


def _d_ao_blocks(nao: int, parts: Sequence[tuple[Sequence[int], np.ndarray]]) -> np.ndarray:
    """Block-diagonal AO density from disjoint subspace occupied MOs."""
    d = np.zeros((nao, nao))
    for aos, c_sub in parts:
        aos = list(aos)
        d[np.ix_(aos, aos)] += 2.0 * (c_sub @ c_sub.T)
    return d


def _d_hl(c_cols: np.ndarray, s_full: np.ndarray) -> np.ndarray:
    """Heitler–London / dual density from (possibly overlapping) occupied MOs."""
    smo = c_cols.T @ s_full @ c_cols
    smo = 0.5 * (smo + smo.T)
    return 2.0 * (c_cols @ np.linalg.inv(smo) @ c_cols.T)


def ground_mf(mol: gto.Mole, method: str = "RHF"):
    if method.upper() == "B3LYP":
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"
    else:
        mf = scf.RHF(mol)
    mf.conv_tol = 1e-8
    mf.max_cycle = 80
    mf.kernel()
    return mf


def construct_states(mf) -> tuple[dict[str, np.ndarray], dict]:
    """Build G/FUD/DSI/PDSI/FUL densities from ground-state Fock."""
    mol = mf.mol
    coords = np.asarray(mol.atom_coords(unit="Angstrom"))
    pi, sigma, detail = identify_pi_sigma_aos(mol, coords)
    if len(pi) != 14:
        raise RuntimeError(f"expected 14 π AOs, got {len(pi)}: {detail}")

    ao_a = set(_atom_aos(mol, FRAG_A))
    ao_bc = set(_atom_aos(mol, FRAG_BC))
    pi_a = [i for i in pi if i in ao_a]
    pi_bc = [i for i in pi if i in ao_bc]
    sig_a = [i for i in sigma if i in ao_a]
    sig_bc = [i for i in sigma if i in ao_bc]
    if len(pi_a) != 6 or len(pi_bc) != 8:
        raise RuntimeError(f"π split A/BC = {len(pi_a)}/{len(pi_bc)}, want 6/8")

    s = np.asarray(mf.get_ovlp(), dtype=float)
    dm_g = np.asarray(mf.make_rdm1(), dtype=float)
    fock = np.asarray(mf.get_fock(dm=dm_g), dtype=float)
    nao = mol.nao

    def sub(aos: Sequence[int], nocc: int) -> tuple[list[int], np.ndarray]:
        aos = list(aos)
        csub, _ = _eigh_sub(fock[np.ix_(aos, aos)], s[np.ix_(aos, aos)], nocc)
        return aos, csub

    pi_all, c_pi = sub(pi, NOCC_PI)
    sg_all, c_sg = sub(sigma, NOCC_SIGMA)
    pia, c_pia = sub(pi_a, NOCC_PI_A)
    pibc, c_pibc = sub(pi_bc, NOCC_PI_BC)
    sa, c_sa = sub(sig_a, NOCC_SIG_A)
    sbc, c_sbc = sub(sig_bc, NOCC_SIG_BC)

    d_g = dm_g
    d_fud = _d_ao_blocks(nao, [(pi_all, c_pi), (sg_all, c_sg)])
    d_dsi = _d_ao_blocks(nao, [(pia, c_pia), (pibc, c_pibc), (sg_all, c_sg)])
    d_ful = _d_ao_blocks(nao, [(pia, c_pia), (pibc, c_pibc), (sa, c_sa), (sbc, c_sbc)])
    c_sig_hl = np.hstack([_embed_mos(nao, sa, c_sa), _embed_mos(nao, sbc, c_sbc)])
    d_pdsi = _d_ao_blocks(nao, [(pia, c_pia), (pibc, c_pibc)]) + _d_hl(c_sig_hl, s)

    dens = {"G": d_g, "FUD": d_fud, "DSI": d_dsi, "PDSI": d_pdsi, "FUL": d_ful}
    meta = {
        "n_pi": len(pi),
        "n_pi_a": len(pi_a),
        "n_pi_bc": len(pi_bc),
        "pi_detail": detail,
        "nocc_sigma_A": NOCC_SIG_A,
        "nocc_sigma_BC": NOCC_SIG_BC,
    }
    return dens, meta


def fock_nonbonded_ss(mf, dm: np.ndarray, coords: np.ndarray) -> float:
    """Occupied nonbonded σ–σ Fock coupling between A and (B+C), hartree.

    Link atoms (N, C_ipso of N-phenyl) excluded so the C–N bond is not
    counted as 'nonbonded'. Public-definition proxy of Enσσ, not Kost LFMO.
    """
    mol = mf.mol
    pi, sigma, _ = identify_pi_sigma_aos(mol, coords)
    fock = np.asarray(mf.get_fock(dm=dm), dtype=float)
    ao_a = set(_atom_aos(mol, FRAG_A))
    ao_bc = set(_atom_aos(mol, FRAG_BC))
    link = set(_atom_aos(mol, [0, 3]))
    sig_a = [i for i in sigma if i in ao_a and i not in link]
    sig_bc = [i for i in sigma if i in ao_bc and i not in link]
    return 2.0 * float(
        np.sum(dm[np.ix_(sig_a, sig_bc)] * fock[np.ix_(sig_a, sig_bc)])
    )


def energy_of_density(mf, dm: np.ndarray) -> float:
    return float(mf.energy_tot(dm=dm))


def energy_effects(energies: dict[str, float]) -> dict[str, float]:
    return {
        "EV": energies["FUD"] - energies["DSI"],
        "En_ss": energies["PDSI"] - energies["FUL"],
        "E_ps": energies["G"] - energies["FUD"],
        "E_ss_full": energies["DSI"] - energies["FUL"],
    }
