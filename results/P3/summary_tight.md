# P3 NBA tightened twist PES (tight_B3LYP_6-31gs_on_RHF_3-21g_both)

- Protocol: A: cascaded free relax RHF/sto-3g(maxiter=250) -> RHF/3-21g(maxiter=80); B: sequential theta sweep (both) at RHF/3-21g (maxiter=40, k0=2.0 Ha/rad^2, stiffening retries until |dtheta|<=3.0 deg, each point seeded from the previous one); C: B3LYP/6-31g* single points; per-angle energy = lower of the two sweeps
- Pre-relaxed seed: theta=149.2 deg, converged=True, dmin=1.070 A
- E_min at theta = **44.9 deg** (dE=-0.965, dEe=-1407.214, dEN=+1406.249 kcal/mol vs theta=-0.1)
- dE span = 1.94 kcal/mol (limit 25.0)
- max |theta - target| = 0.19 deg (limit 3.0)
- topology intact at every point: True
- converged points: 5/7; max up/down hysteresis = 0.275442099676744 kcal/mol (limit 1.0)
- **Quality gate: PASSED**
- Checks: {'E_min_in_30_60': True, 'at_min_EN_up_Ee_down': True, 'near_planar_not_global_min': True}
- Auto agree flag: **True** (None = gate failed; formal VERDICT only in deliverables/)
