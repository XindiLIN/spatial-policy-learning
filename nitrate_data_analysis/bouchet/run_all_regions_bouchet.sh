#!/bin/bash
# A previous sequential-across-regions version of this script ran for real
# (job 20985635): 4h06m wall-clock but only 6.40% CPU efficiency (9.19GB
# memory, 7.18% of 128G) -- confirming piece (iv)'s single-threaded DC
# algorithm left most of a 16-core allocation idle through each region's
# entire DC phase, one region at a time. The script has since been
# restructured to parallelize across (area, threshold) pairs instead (up to
# 16 concurrent, each getting a share of the remaining cores for piece
# (iii)'s own internal parallelism), which should fix that -- but this
# specific restructuring hasn't been run on Bouchet yet, so -c/--mem/-t here
# are still an estimate. threshold_vals now covers 3 thresholds (up to 27
# (area, threshold) tasks total, vs. 9 before); rough math (CPU-time scaling
# ~3x from the old run, now spread across up to 16 concurrent tasks instead
# of ~1 effective core) suggests wall-clock could plausibly drop to 1.5-3.5h
# even with 3x the work -- kept -t at 12:00:00 as headroom, not a tight
# estimate. Memory may be higher than the old 9.19GB now that multiple
# regions' Gram matrices are in memory simultaneously instead of one at a
# time -- check sacct/seff after this run and adjust accordingly.
#SBATCH -J run_all_regions
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 12:00:00
#SBATCH -o slurm-allregions-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xindi.lin@yale.edu

cd "$SLURM_SUBMIT_DIR"

module reset
module load R/4.4.2-gfbf-2024a

# See run_cv_hyperparameter_bouchet.sh for why both of these matter: OPENBLAS
# is already this module's default FlexiBLAS backend (this export is just for
# explicitness), and OMP_NUM_THREADS=1 prevents OpenBLAS's own internal
# multithreading from oversubscribing the core allocation on top of this
# script's own mclapply-based parallelism (piece iii's per-observation warm
# start, via init_mc.cores).
export FLEXIBLAS=OPENBLAS
export OMP_NUM_THREADS=1
Rscript -e 'writeLines(grep("BLAS", capture.output(sessionInfo()), value = TRUE))'

Rscript nitrate_data_analysis/run_all_regions.R
