#!/bin/bash

# Adjust this path to your actual MATLAB installation:
MATLABCMD="matlab"

: "${PSIPET_MODEL_ROOT:?Set PSIPET_MODEL_ROOT to the run folder (see README.md)}"
cd "$(dirname "$0")"

# -------- defaults (can be overridden by flags) --------
GROUP='active'
TASK='rest'
SES='02'
Tau=1
ATLASNAME='116'
FEATURE='corrMatrixICA24Motion'
GRIDON=1

# -------------------------------------------------------

# Parse CLI flags (minimal, POSIX-friendly)
usage() {
  echo "Usage: $0 [-group GROUP] [-TASK TASK] [-SES SES] [-Tau N] [-ATLASNAME NAME] [-FEATURE FEATURE] [-GRIDON GRIDON]"
  echo "Defaults: Tau=1, ATLASNAME='116', GROUP='Active', TASK='rest', SES='01', FEATURE='corrMatrixICA24Motion', GRIDON=1"
  exit 1
}

while [ $# -gt 0 ]; do
  case "$1" in
    -group)     GROUP="$2"; shift 2;;
    -TASK)      TASK="$2"; shift 2;;
    -SES)       SES="$2"; shift 2;;
    -Tau)       Tau="$2"; shift 2;;
    -GRIDON)    GRIDON="$2"; shift 2;;
    -FEATURE)    FEATURE="$2"; shift 2;;
    -ATLASNAME) ATLASNAME="$2"; shift 2;;
    -h|--help)  usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

# 2) Log file name
LOGFILE="logs/hopf_${GROUP}_${TASK}_${SES}_tau-${Tau}_atlas-${ATLASNAME}_feature-${FEATURE}_grid-${GRIDON}.log"
mkdir -p "$(dirname "$LOGFILE")"

echo "=========================================================="
echo "Launching Hopf model scripts with these vars:"
echo "  GROUP     = $GROUP"
echo "  TASK      = $TASK"
echo "  SES       = $SES"
echo "  ATLAS     = $ATLASNAME"
echo "  FEATURE   = $FEATURE"
echo "  Tau       = $Tau"
echo "  GRIDON    = $GRIDON"
echo "Logging to: $LOGFILE"
echo "=========================================================="

# %%% run('run_hopf_grid_group_psipet.m'); \
# %%% run('run_hopf_subject_refinement_psipet.m'); \

# 3) Launch MATLAB
nohup "$MATLABCMD" -nodisplay -nosplash -r "\
try; \
    GROUP='$GROUP'; \
    ATLASNAME='$ATLASNAME'; \
    Tau=$Tau; \
    TASK='$TASK'; \
    SES='$SES'; \
    GRIDON=$GRIDON; \
    FEATURE='$FEATURE'; \
    fprintf('Running MATLAB scripts for: GROUP=%s, TASK=%s, SES=%s, Tau=%d, ATLAS=%s, FEATURE=%s\n', GROUP, TASK, SES, Tau, ATLASNAME, FEATURE); \
    run('run_hopf_grid_group_psipet.m'); \
    run('run_hopf_subject_refinement_psipet.m'); \
catch ME; \
    disp('Error in MATLAB script:'); \
    disp(getReport(ME)); \
end; \
exit;" \
> "$LOGFILE" 2>&1 &

echo "MATLAB job submitted in background (PID: $!)."
echo "Use 'tail -f $LOGFILE' to watch progress."
