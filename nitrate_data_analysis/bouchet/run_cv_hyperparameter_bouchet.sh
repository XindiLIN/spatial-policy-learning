#!/bin/bash
# -c 16/--mem 128G validated against a real 9-region, 1-threshold run
# (job 20941012): 8h07m wall-clock, 42.77% CPU efficiency, 45.08GB memory
# (35.21% of 128G) -- see sacct/seff output and git history for the full
# diagnosis. cv_thresholds now covers 3 thresholds (log(2), log(5), log(10))
# instead of 1, tripling the (area, threshold, m, combo) task count to
# ~3240 -- CPU-time should scale roughly 3x too (~167 core-hours), so even
# at the same efficiency, wall-clock could plausibly exceed 12h; bumped -t
# accordingly. --mem should NOT need to grow: k_raw (the dominant per-task
# memory cost) depends only on (area, m), not threshold, so it isn't
# duplicated per threshold. Checkpointing means a walltime kill isn't
# catastrophic either way -- just resubmit to pick up where it left off.
#SBATCH -J cv_hyperparameter
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o slurm-cv-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xindi.lin@yale.edu

cd "$SLURM_SUBMIT_DIR"

module reset
module load R/4.4.2-gfbf-2024a

# OPENBLAS is already this R module's default FlexiBLAS backend (confirmed
# via sessionInfo() -- NETLIB is available but not the default), so this
# export is just for explicitness/reproducibility, not a fix by itself. Note
# the name is "OPENBLAS", not "OPENBLAS-OPENMP" -- `flexiblas list` displays
# the latter (matching its library filename), but that's not the name the
# FLEXIBLAS env var actually accepts; using it silently fell back to the
# (identical, correct) default with a warning on stderr. The actual fix is
# OMP_NUM_THREADS=1: OpenBLAS multithreads internally by default, and
# cv_hyperparameter.R separately parallelizes across many forked R processes
# via mclapply (both across areas and across (area, m, combo) tasks) --
# without this, each of those processes would ALSO try to spawn multiple
# OpenBLAS threads, oversubscribing the core allocation (N workers x M BLAS
# threads each) rather than speeding anything up.
export FLEXIBLAS=OPENBLAS
export OMP_NUM_THREADS=1
Rscript -e 'writeLines(grep("BLAS", capture.output(sessionInfo()), value = TRUE))'

Rscript nitrate_data_analysis/cv_hyperparameter.R
