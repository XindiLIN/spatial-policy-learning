#!/bin/bash
#SBATCH -J install_r_packages
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH -o slurm-install-%j.out

module reset
module load R/4.4.2-gfbf-2024a
Rscript simulation/install_packages.R
