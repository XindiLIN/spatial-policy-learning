#!/bin/bash
# -c/--mem mirror run_cv_hyperparameter_bouchet.sh's validated allocation as a
# starting point -- this script hasn't itself been run on Bouchet yet, so
# treat -t 12:00:00 as generous headroom, not a validated estimate. Regions
# are processed sequentially (see run_all_regions.R's main loop), each one
# using the full -c allocation for piece (iii)'s per-observation DC warm
# start, so memory needs should if anything be lower than the CV job's (which
# held Gram matrices for every (area, m) pair in kraw_by_area_m at once,
# vs. one region/one m here at a time). Check sacct/seff after a real run and
# adjust -c/--mem/-t the same way phase ABC's and the CV job's were tuned.
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
