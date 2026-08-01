# Non-interactive package install for the nitrate_data_analysis pipeline -- explicitly creates
# and uses a personal library so this works standalone as a batch job (there's no terminal to
# answer an interactive "create a personal library?" prompt when run via Rscript in an sbatch
# job). Mirrors simulation/install_packages.R but covers the nitrate pipeline's own dependency
# list (region_pipeline_funcs.R, functions/miscellaneous.r, cv_hyperparameter.R,
# run_all_regions.R, policy_visualization.r, dataset_visualization.R), which only partially
# overlaps with the simulation pipeline's -- e.g. WeightIt (CBPS propensity weights) is required
# here but not by the simulation code, and isn't in simulation/install_packages.R's list.
lib_dir <- Sys.getenv("R_LIBS_USER")
if (lib_dir == "") {
  lib_dir <- file.path(Sys.getenv("HOME"), "R",
                        paste0(R.version$platform, "-library"),
                        paste(R.version$major, sub("\\..*", "", R.version$minor), sep = "."))
}
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(lib_dir)

install.packages(
  c("GpGp", "Metrics", "WeightIt", "caret", "dplyr", "e1071", "fields", "ggplot2", "glmnet",
    "kernlab", "recipes", "sf", "yardstick", "emmeans", "patchwork", "tigris", "tidyr",
    "stringr"),
  lib = lib_dir,
  repos = "https://cloud.r-project.org"
)
