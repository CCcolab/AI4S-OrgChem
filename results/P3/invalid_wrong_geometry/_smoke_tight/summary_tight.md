# P3 NBA tightened twist PES (tight_RHF_sto-3g_on_RHF_sto-3g)

- Protocol: A: cascaded free relax RHF/sto-3g(maxiter=2); B: theta-constrained relax from that seed at RHF/sto-3g (maxiter=2, k0=2.0 Ha/rad^2, stiffening retries until |dtheta|<=3.0 deg); C: RHF/sto-3g single points
- Pre-relaxed seed: theta=56.5 deg, converged=False, dmin=0.976 A
- E_min at theta = **44.5 deg** (dE=-50.920, dEe=-1119.745, dEN=+1068.825 kcal/mol vs theta=0.8)
- dE span = 50.92 kcal/mol (limit 25.0)
- max |theta - target| = 0.76 deg (limit 3.0)
- **Quality gate: FAILED**
- Checks: {'E_min_in_30_60': True, 'at_min_EN_up_Ee_down': True, 'near_planar_not_global_min': True}
- Auto agree flag: **None** (None = gate failed; formal VERDICT only in deliverables/)
