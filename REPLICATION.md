# Replication record

**2026-09-24.** Full re-run of the model-fitting pipeline on an independent
machine (compute.kcir.se, MATLAB R2026a, 26.1.0.3203278, 64-bit glnxa64), from a
fresh copy of the code with inputs uploaded from the local copy of the original
HALFpipe timeseries export. The MATLAB files in this repository are
byte-identical to those used for that run:

```
    34d5ad4b45b85fda7b509367b08cfc77314ae55b83fdc910c8af15c8e29b5fb1  matlab/run_hopf_all_psipet.sh
    95c4f023c17f3ef9fb00b777ab3781ecaa7e1a5c2b316cb3a2d882f62a1ef4c4  matlab/run_hopf_grid_group_psipet.m
    8bbc6cab89091c4d4fb010077508464036def4ae1f65a7316dc171bb0db4668a  matlab/run_hopf_subject_refinement_psipet.m
    901e7698942deeb032108130155742c3deabc3adbcd533cf78ed459eaf8cabe6  matlab/fit_group_level.m
    c9d1aa729d48c2bcceaad5059075ed9c634f9b92b9c997bd09473350cc4c7106  matlab/compute_group_empirical.m
```

(`sha256sum matlab/*.m matlab/*.sh` to check.)

- The grid search (g = 0.20 to 0.80 in steps of 0.005, a = -0.02) reselected
  g = 0.435 / 0.420 / 0.405 / 0.400 (psilocybin ses-01/02, placebo ses-01/02)
  with total errors 0.03174 / 0.026982 / 0.030844 / 0.03173, identical to the
  published fits.
- Subject refinement was re-run at the reselected g, and the analysis notebook
  was executed against the replicated results with outputs redirected.
- Every CSV the notebook produces (group tables, eTables, permutation results)
  was byte-identical to the published outputs. Figures were not compared at
  byte level (they embed timestamps); every number feeding the manuscript is in
  the CSVs.

The model is deterministic (linear noise approximation via the continuous-time
Lyapunov equation; no stochastic simulation), so exact reproduction is the
expected outcome for correct code and inputs. Byte-identity was achievable, not
trivial: it required the same inputs, the same code, and the grid actually
being searched.
