#!/bin/bash
#SBATCH -J nonparametric_sim_abc
#SBATCH -c 64
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o slurm-abc-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xindi.lin@yale.edu

cd "$SLURM_SUBMIT_DIR"

module reset
module load R/4.4.2-gfbf-2024a
Rscript simulation/simulation_nonparametric_run_phaseABC.R
