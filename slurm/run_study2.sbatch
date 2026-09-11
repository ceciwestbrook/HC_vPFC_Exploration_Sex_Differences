#!/bin/bash
#SBATCH --job-name=study2_models
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=180G
# Log paths are set below from $DATA_ROOT rather than as #SBATCH directives,
# because SBATCH directives cannot use variables.

# Adjust partition/cluster to your allocation, e.g.:
# #SBATCH --cluster=smp
# #SBATCH --partition=smp

### Submits analysis/03_study2_models.R as a batch job.
###
###   sbatch slurm/run_study2.sbatch
###
### Submit from the repository root. Set DATA_ROOT if your data is not where
### config.R says it is; otherwise the two are expected to agree.

set -u

REPO="${REPO:-$(cd "$(dirname "$0")/.." && pwd)}"
DATA_ROOT="${DATA_ROOT:-/ix/cladouceur/DNPL_sex_differences_paper}"
LOGDIR="$DATA_ROOT/BSOC/fMRI/rerun_logs"
mkdir -p "$LOGDIR"
exec >>"$LOGDIR/study2_${SLURM_JOB_ID:-manual}.out" 2>&1

module purge
# R 4.6.0 — must match the version your personal library was built against.
# Verify the exact module name with:  module avail r
module load gcc/12.2.0 r/4.6.0

# Personal library holding fmri.pipeline and the rest. Set R_LIBS_USER in your
# environment, or edit this line; the default is R's own user library.
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.6}"

export NCORES="${SLURM_CPUS_PER_TASK:-20}"
# export FORCE=1        # overwrite outputs already written today
# export SAVE_DATA=1    # also save the model frames (large)

# sbatch inherits the submitting shell's environment, so a stray VERIFY_ONLY
# from an interactive session would silently stop the job after the gates.
unset VERIFY_ONLY

SCRIPT="$REPO/analysis/03_study2_models.R"
if [ ! -f "$SCRIPT" ]; then echo "missing: $SCRIPT"; exit 1; fi
if [ ! -f "$REPO/config.R" ]; then echo "missing: $REPO/config.R"; exit 1; fi

cd "$REPO/analysis" || exit 1      # so the config lookup finds ../config.R
echo "repo:      $REPO"
echo "data root: $DATA_ROOT"
echo "cores:     $NCORES"
Rscript "$SCRIPT"
status=$?
echo "Rscript exit status: $status"
exit $status
