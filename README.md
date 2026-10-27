# psipet-hierarchy

Analysis code for *Directed Brain Hierarchy and Symptom Improvement After
Psilocybin in Major Depression: A Secondary Analysis of a Randomized Clinical
Trial*.

Resting-state fMRI and diffusion MRI were acquired before and a median of two
days after a single 25 mg dose of psilocybin or 100 mg niacin (active placebo)
in adults with major depressive disorder (ClinicalTrials.gov NCT04630964).
Whole-brain effective connectivity was estimated with a linearized Hopf model
fitted to functional connectivity and lagged covariance, and its hierarchical
organization quantified by trophic analysis.

## Layout

```
config.py              default paths and analysis switches
config_local.py        your machine's paths (not tracked; overrides config.py)
notebooks/main.ipynb   the analysis: all figures, tables and statistics
src/scan_intervals.py  dosing-to-MRI intervals reported in Methods
matlab/                model fitting (run first, produces the .mat inputs)
matlab/external/       third-party dependencies you must obtain yourself
tools/                 export the working notebook into this repository
results/figures/       generated figures
results/tables/        generated tables
```

## Running it

1. **Fit the model** (MATLAB). Copy `HopfModel_LNA.m` and `gaussfilt.m` into
   `matlab/external/` — see the README there — and lay the inputs out as

   ```
   $PSIPET_MODEL_ROOT/
   ├── data/
   │   ├── atlas_SchaeferTian116v3_sift_invnodevol_radius2_count_connectivity_ses-01.mat  (and ses-02)
   │   ├── PSIPET-ID_group-active_task-rest.mat  (and placebo)
   │   └── tseries/rest/ses-01/…_feature-corrMatrixICA24Motion_atlas-SchaeferTian116_timeseries.tsv  (and ses-02)
   └── results/            created by the scripts, with data/derived_data/
   ```

   then:

   ```bash
   export PSIPET_MODEL_ROOT=/path/to/run/folder
   cd matlab
   ./run_hopf_all_psipet.sh -group 'active'  -TASK 'rest' -SES '01' -Tau 1 -GRIDON 1
   ./run_hopf_all_psipet.sh -group 'active'  -TASK 'rest' -SES '02' -Tau 1 -GRIDON 1
   ./run_hopf_all_psipet.sh -group 'placebo' -TASK 'rest' -SES '01' -Tau 1 -GRIDON 1
   ./run_hopf_all_psipet.sh -group 'placebo' -TASK 'rest' -SES '02' -Tau 1 -GRIDON 1
   ```

   Each run grid-searches g (0.20–0.80 in steps of 0.005, bifurcation
   parameter fixed at −0.02, 121 fits) and then refits every subject at the
   selected g. Saved per run: the selected `gBest`, `Ceff_best`, the error
   surface, convergence trajectories, and the simulated and empirical FC and
   COVτ at the selected solution (`FCsim_best`, `COVtauSim_best`,
   `FCemp_group`, `COVtauemp_group`).

   `GRIDON 1` runs the grid search over global coupling, g = 0.20 to 0.80 in
   steps of 0.005, with the bifurcation parameter fixed at a = -0.02. Selected
   values were g = 0.435 and 0.420 (psilocybin, before and after) and 0.405 and
   0.400 (placebo).

2. **Run the analysis** (Python).

   ```bash
   pip install -r requirements.txt
   export PSIPET_ROOT=/path/to/your/data
   jupyter lab notebooks/main.ipynb
   ```

   Run the cells in order. Everything in the manuscript and supplement is
   produced by this notebook, writing to `results/`, with one exception:

3. **Scan intervals.** The dosing-to-MRI intervals in Methods come from

   ```bash
   python src/scan_intervals.py
   ```

   It writes aggregate statistics to `results/tables/scan_intervals_summary.csv`.
   The per-participant table contains treatment and scan dates and is written
   next to the source spreadsheet, never inside the repository.

## Replication

The pipeline in this repository was re-run end to end on an independent
machine and reproduced every published number byte for byte. See
`REPLICATION.md`.

## Data availability

The inputs are **not** distributed here: the fitted connectivity matrices
(`.mat`), the clinical outcome data, the motion and tSNR derivatives, and the
participant identifier files. They contain individual-level data from a clinical
trial. Requests for the derived model outputs should be addressed to the
corresponding author.

The principal functional connectivity gradient is fetched at runtime by
`neuromaps` (Margulies et al. 2016). The parcellation is Schaefer-100 combined
with 16 subcortical regions from the Tian S1 atlas, and is not redistributed.

## Notes on the code

- Sync the notebook with `python tools/export_notebook.py`, which reads
  `WORKING_NOTEBOOK` from `config_local.py`. It refuses to write if any
  absolute path remains in the source. Add `--keep-outputs` to publish an
  executed notebook; participant identifiers in the outputs are replaced by
  stable pseudonyms (P01, P02, ...) unless you pass `--no-pseudonymise`. Add
  `--with-results` to refresh `results/` from the directories the notebook
  wrote to.
- `notebooks/main.ipynb` is self-contained. Every function it calls is defined
  within it; the only imports are the packages in `requirements.txt`.
- Outputs are cleared from the committed notebook, both to keep diffs readable
  and because executed outputs contained participant identifiers.
- The MATLAB scripts here differ from the archived originals in three
  deliberate ways: the run folder comes from `PSIPET_MODEL_ROOT` instead of a
  hardcoded server path; the grid script performs the full search (the archived
  copy had the selected g values pinned for trajectory refits); and per-cell
  convergence trajectories plus the simulated FC/COVτ at the selected solution
  are saved. Model fitting itself is unchanged, so a re-run reproduces the
  published `gBest` values and fits.
- `CONFIG = 'grid_ga'` in `config.py` selects a superseded two-dimensional
  search over both g and a. It is retained for reference only. The manuscript
  reports `grid_g`.

## License

MIT, see `LICENSE`. This covers the code in this repository only, not the
third-party functions described in `matlab/external/README.md`.
