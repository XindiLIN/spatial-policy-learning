library(dplyr)
library(tidyr)

output_dir <- "nitrate_data_analysis/output"
area_results_dir <- file.path(output_dir, "area_results")

areas <- c("central", "east_central", "north_central", "north_east", "north_west",
           "south_central", "south_east", "south_west", "west_central")
threshold_labels <- c("log2", "log5", "log10")

# For each (area, threshold), pull the mean-of-regions metrics for direct/indirect x train/test.
metrics_long <- lapply(threshold_labels, function(threshold_label) {
  rows <- lapply(areas, function(area_key) {
    path <- file.path(area_results_dir, sprintf("region_result_%s_%s.rds", area_key, threshold_label))
    if (!file.exists(path)) return(NULL)
    r <- readRDS(path)
    bind_rows(
      as_tibble(r$direct_metrics_train)   %>% mutate(area_key = area_key, threshold_label = threshold_label, method = "Our Method",  split = "Training"),
      as_tibble(r$direct_metrics_test)    %>% mutate(area_key = area_key, threshold_label = threshold_label, method = "Our Method",  split = "Test"),
      as_tibble(r$indirect_metrics_train) %>% mutate(area_key = area_key, threshold_label = threshold_label, method = "Non-spatial", split = "Training"),
      as_tibble(r$indirect_metrics_test)  %>% mutate(area_key = area_key, threshold_label = threshold_label, method = "Non-spatial", split = "Test")
    )
  })
  bind_rows(rows)
}) %>% bind_rows()

# Average across regions (dropping NA MCC folds, e.g. indirect method's boundary-clamp
# collapse at log10 for some regions -- see East Central at log10 for the clearest case).
metrics_summary <- metrics_long %>%
  group_by(threshold_label, method, split) %>%
  summarise(
    n_regions_mcc = sum(!is.na(mcc)),
    acc = mean(acc, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE),
    two_sided_f1 = mean(two_sided_f1, na.rm = TRUE),
    .groups = "drop"
  )

print(metrics_summary, n = Inf)

# Wide table matching the paper's presentation: rows = metric, columns = threshold x method x split
metrics_wide <- metrics_summary %>%
  pivot_longer(cols = c(acc, mcc, two_sided_f1), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric, acc = "Accuracy", mcc = "MCC", two_sided_f1 = "Two-sided F1"),
    threshold_label = factor(threshold_label, levels = c("log2", "log5", "log10"))
  ) %>%
  select(threshold_label, metric, method, split, value) %>%
  pivot_wider(names_from = c(threshold_label, split, method), values_from = value)

print(metrics_wide)

saveRDS(metrics_summary, file.path(output_dir, "area_results", "policy_metrics_comparison_summary.rds"))
