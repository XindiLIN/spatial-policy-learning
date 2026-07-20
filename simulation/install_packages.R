# Non-interactive package install -- explicitly creates and uses a personal library so this
# works standalone as a batch job (there's no terminal to answer an interactive "create a
# personal library?" prompt when run via Rscript in an sbatch job).
lib_dir <- Sys.getenv("R_LIBS_USER")
if (lib_dir == "") {
  lib_dir <- file.path(Sys.getenv("HOME"), "R",
                        paste0(R.version$platform, "-library"),
                        paste(R.version$major, sub("\\..*", "", R.version$minor), sep = "."))
}
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(lib_dir)

install.packages(
  c("tidymodels","tidyverse","SuperLearner","dplyr","ggplot2","GpGp","Metrics","yardstick","glmnet","kernlab",
    "fields","sp","gstat","caret","e1071","recipes"),
  lib = lib_dir,
  repos = "https://cloud.r-project.org"
)
