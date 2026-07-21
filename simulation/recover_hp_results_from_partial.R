library(dplyr)

# If simulation_nonparametric_run_phaseABC.R's Phase C crashes/OOMs partway through, this
# reconstructs hp_results_all from whatever per-task checkpoint files it managed to write to
# simulation/output/hp_partial/ before dying, so the completed work isn't lost. Requires
# simulation_nonparametric_thr_ctxs.rds, which phaseABC.R now saves right after Phase B --
# before Phase C runs -- specifically so it survives a Phase C failure.
#
# Note the result will only cover whichever (seed, threshold, m) tasks finished before the
# crash; it is not a substitute for a completed Phase C run, but it's real, useable data for
# whatever tasks it does cover, and can save you from redoing that work.

thr_ctxs_path <- 'simulation/output/simulation_nonparametric_thr_ctxs.rds'
hp_partial_dir <- 'simulation/output/hp_partial'

if (!file.exists(thr_ctxs_path)) {
  stop(sprintf("Missing %s -- can't recover without it (needed to attach seed/threshold_quantile).", thr_ctxs_path))
}

partial_files <- list.files(hp_partial_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(partial_files) == 0) {
  stop(sprintf("No partial checkpoint files found in %s.", hp_partial_dir))
}

thr_ctxs <- readRDS(thr_ctxs_path)

hp_results <- lapply(partial_files, readRDS)
hp_results_all <- bind_rows(hp_results)

hp_results_all$seed <- vapply(hp_results_all$thr_idx, function(i) thr_ctxs[[i]]$seed_value, numeric(1))
hp_results_all$threshold_quantile <- vapply(hp_results_all$thr_idx, function(i) thr_ctxs[[i]]$threshold_quantile, numeric(1))
hp_results_all$selected <- FALSE
for (idx in unique(hp_results_all$thr_idx)) {
  rows <- which(hp_results_all$thr_idx == idx)
  best <- rows[which.max(hp_results_all$mcc[rows])]
  hp_results_all$selected[best] <- TRUE
}

n_tasks_expected <- length(thr_ctxs) * 3  # 3 values of m

cat(sprintf("Recovered %d of an expected %d Phase C tasks (%d combo-level rows) from %s.\n",
            length(partial_files), n_tasks_expected, nrow(hp_results_all), hp_partial_dir))
cat("NOTE: 'selected' (best combo per seed/threshold) is only meaningful for (seed, threshold)\n")
cat("pairs where all 3 m-groups finished -- check which thr_idx values have all 3 m's present\n")
cat("before trusting 'selected' for downstream use (e.g. feeding into phaseD.R).\n")

dir.create('simulation/output', recursive = TRUE, showWarnings = FALSE)
saveRDS(hp_results_all, file = 'simulation/output/simulation_nonparametric_hp_results.rds')
