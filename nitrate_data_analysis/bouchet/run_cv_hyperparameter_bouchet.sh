#!/bin/bash
# -c/--mem below are an INITIAL, UNVALIDATED estimate -- unlike
# simulation/bouchet/run_nonparametric_phaseABC_bouchet.sh's -c 32/--mem=128G,
# which were tuned against real sacct/seff OOM history, this script hasn't
# been run on Bouchet yet. cv_hyperparameter.R currently defaults to
# cv_areas <- c("North East") (the smallest region), giving 6 (area, m) tasks
# total for the full m_grid -- mc.cores beyond ~6 won't add parallelism for
# that config yet. --mem is set generously ahead of need (128G, not 64G) for
# two reasons that will both matter once cv_areas expands beyond one region:
# (1) each of the up to 9 areas can run as its own concurrent (area, m) task
# once cv_areas has more than one entry -- more simultaneous workers than
# just the m_grid alone -- and (2) real nitrate regions are larger than the
# simulation's synthetic data (e.g. Central's n_train ~4244 vs the
# simulation's ~1750), so each worker's Gram matrix/smoothers footprint is
# bigger too. When cv_areas is later expanded toward all 9 regions (up to 54
# tasks), revisit -c/--mem the same way phase ABC's were tuned: check
# sacct/seff after a real run and adjust (scale --mem roughly proportionally
# with -c, not just in total, per that script's own note).
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
Rscript nitrate_data_analysis/cv_hyperparameter.R
