#!/bin/bash
#SBATCH -J nonparametric_sim
#SBATCH -c 20
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o slurm-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xindi.lin@yale.edu

cd "$SLURM_SUBMIT_DIR"

module reset
module load R/4.4.2-gfbf-2024a
Rscript simulation/simulation_nonparametric_run.R
