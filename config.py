"""Every filesystem path used by the analysis, in one place.

Two ways to point this at your data:

  1. Set PSIPET_ROOT and lay the data out as described in README.md, or
  2. Create `config_local.py` next to this file, assigning whichever of the
     constants below you need. It is imported last and wins. It is listed in
     .gitignore, so machine-specific paths never enter version control.

The inputs these paths point at are not distributed with this repository
(see README.md, "Data availability").
"""
import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Analysis configuration
# ---------------------------------------------------------------------------
SEED = 42
N_SUBCORTICAL = 16

# 24 = every participant with both sessions (no outlier removal)  <- manuscript
# 22 = MAD-median rule applied to delta trophic coherence
N_SUBJECTS = 24
SUBJECT_FILTER = {24: "paired_sessions", 22: "trophic_outlier"}[N_SUBJECTS]

# 'grid_g'  = coupling only, bifurcation parameter fixed at -0.02  <- manuscript
# 'grid_ga' = superseded two-dimensional search over g and a
CONFIG = "grid_g"
TAG = {"grid_g": "grid-g", "grid_ga": "grid-ga"}[CONFIG]

# ---------------------------------------------------------------------------
# Default layout, relative to PSIPET_ROOT
# ---------------------------------------------------------------------------
ROOT = Path(os.environ.get("PSIPET_ROOT", Path.home() / "psipet"))
_DATA, _MODEL, _ATLAS = ROOT / "data", ROOT / "model", ROOT / "atlas"

_RES = {"grid_g":  _MODEL / "results-gridon-only-g" / "results",
        "grid_ga": _MODEL / "results-gridon"}[CONFIG]

RES_PATH  = f"{_RES / 'results_subject_refinement'}/"           # per-subject refits
GRID_PATH = f"{_RES / 'results_group_only'}/"                    # grid search output
TRAJ_PATH = f"{_RES / 'results_best_traj' / 'results_group_only'}/"  # refit at selected g

_OUT      = Path(__file__).resolve().parent / "results"
SAVEFIG   = f"{_OUT / 'figures' / f'{TAG}_n{N_SUBJECTS}'}/"
SAVETABLE = f"{_OUT / 'tables'}/"

ATLAS_DLABEL   = str(_ATLAS / "schaefertian100" /
                     "Schaefer2018_100Parcels_7Networks_order_Tian_Subcortex_S1.dlabel.nii")
PET_TO_MR_XLSX = str(ROOT / "PSIPET PET to MR mapping.xlsx")
MADRS_XLSX     = str(_DATA / "BEH" / "MADRS_treatment.xlsx")
MADRSS_XLSX    = str(_DATA / "BEH" / "MADRSS.xlsx")       # self-rated, daily
MOTION_DIR     = f"{_DATA / 'MOTION'}/"
TS_DIR         = f"{ROOT / 'ts'}/"
JSON_DIR       = str(_DATA / "JSON-files")
DWI_DIR        = str(_DATA / "qsirecon-MRtrix3_fork-SS3T_act-None")
DWI_QC_CSV     = str(_DATA / "qsiprep_fd_qc_subject_level.csv")
SCAN_DATES_XLSX = str(ROOT / "scan_dates.xlsx")   # src/scan_intervals.py

# ---------------------------------------------------------------------------
# Machine-specific overrides (not tracked)
# ---------------------------------------------------------------------------
try:
    from config_local import *          # noqa: F401,F403
except ImportError:
    pass

for _d in (Path(SAVEFIG), Path(SAVETABLE)):
    try:
        _d.mkdir(parents=True, exist_ok=True)
    except OSError:
        pass   # output root not mounted; the notebook will report on first write
