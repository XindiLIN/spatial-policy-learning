#!/bin/bash
#SBATCH -J nonparametric_sim
#SBATCH -c 20
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o slurm-%j.out

module reset
module load R/4.4.1-foss-2022b
Rscript simulation/simulation_nonparametric_run.R
