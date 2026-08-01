#!/bin/bash
# -c/--mem history: -c 64 --mem=64G and -c 64 --mem=128G both got OOM-killed during phase C
# (job 19064715's sacct showed MaxRSS hit ~128G, right at the cap). Dropping to -c 32 while
# keeping --mem=128G (job 19089256) finally avoided OOM: seff reported Memory Utilized
# 121.82 GB / 128.00 GB (95.17%) -- ~3.8 GB/worker actual usage, with exit code 0. So
# -c 32 --mem=128G is validated as workable; that run instead hit the 12h time limit
# (TIMEOUT) partway through phase C, which is fine since phase C is resumable (see
# hp_partial/ checkpointing in simulation_nonparametric_run_phaseABC.R) -- just resubmit.
# If -c is ever raised again, scale --mem up proportionally (~4 GB/core), not just in total.
#SBATCH -J nonparametric_sim_abc
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH -t 12:00:00
#SBATCH -o slurm-abc-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xindi.lin@yale.edu

cd "$SLURM_SUBMIT_DIR"

module reset
module load R/4.4.2-gfbf-2024a
Rscript simulation/simulation_nonparametric_run_phaseABC.R
