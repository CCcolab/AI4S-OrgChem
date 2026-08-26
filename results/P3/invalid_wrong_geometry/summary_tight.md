# P3 NBA tightened twist PES (tight_B3LYP_6-31gs_on_RHF_3-21g)

- Protocol: A: cascaded free relax RHF/sto-3g(maxiter=250) -> RHF/3-21g(maxiter=80); B: theta-constrained relax from that seed at RHF/3-21g (maxiter=20, k0=2.0 Ha/rad^2, stiffening retries until |dtheta|<=3.0 deg); C: B3LYP/6-31g* single points
- Pre-relaxed seed: theta=-92.2 deg, converged=True, dmin=1.066 A
- E_min at theta = **90.0 deg** (dE=-5.356, dEe=+10271.546, dEN=-10276.902 kcal/mol vs theta=-0.3)
- dE span = 14.44 kcal/mol (limit 25.0)
- max |theta - target| = 0.31 deg (limit 3.0)
- **Quality gate: PASSED**
- Checks: {'E_min_in_30_60': False, 'at_min_EN_up_Ee_down': False, 'near_planar_not_global_min': True}
- Auto agree flag: **False** (None = gate failed; formal VERDICT only in deliverables/)
