#!/bin/bash
# -c 16/--mem 128G validated against a real run: cv_hyperparameter.R with
# cv_areas <- c("North East") only (120 (area, m, combo) tasks) completed in
# 7m21s (sacct), most of that setup rather than the CV grid itself -- see the
# git history for the diagnosis (mclapply parallelized only at the (area, m)
# level, which underused this allocation, plus a since-confirmed-unnecessary
# BLAS-backend theory) and the fixes (combo-level task granularity,
# OMP_NUM_THREADS=1 below, and parallelizing per-area setup across areas).
# cv_areas now covers all 9 regions (up to 1080 (area, m, combo) tasks) --
# if --mem ever proves insufficient for the larger regions' Gram matrices
# (e.g. Central's n_train ~4244 vs North East's ~1161), check sacct/seff
# after a run and scale --mem up from there, same as phase ABC's tuning.
#SBATCH -J cv_hyperparameter
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 12:00:00
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
