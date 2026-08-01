library(tigris)
library(ggplot2)
library(tidyr)
library(tidymodels)
library(dplyr)
library(patchwork)
library(stringr)
library(emmeans)
library(sf)

source("functions/miscellaneous.r")


# Nonsp Indirect method v.s. Our method (acc/mcc)

policy_comparison_combine <- readRDS("nitrate_data_analysis/output/policy_comparison_combine.rds")

policy_comparison_combine_sf <- st_as_sf(policy_comparison_combine, coords = c("longitude", "latitude"), crs = 4326)
policy_comparison_combine_sf_test <- policy_comparison_combine_sf[policy_comparison_combine_sf$test == 1, ]

train_nonspatial <- calculate_acc_mcc_two_sided_f1(T = log(policy_comparison_combine_sf$WellDepth),
                                                   Y = log(policy_comparison_combine_sf$concentration_plus_median),
                                                   policy = policy_comparison_combine_sf$indirect_policy,
                                                   threshold_val = log(5))

test_nonspatial <- calculate_acc_mcc_two_sided_f1(T = log(policy_comparison_combine_sf_test$WellDepth),
                                                  Y = log(policy_comparison_combine_sf_test$concentration_plus_median),
                                                  policy = policy_comparison_combine_sf_test$indirect_policy,
                                                  threshold_val = log(5))

train_our_method <- calculate_acc_mcc_two_sided_f1(T = log(policy_comparison_combine_sf$WellDepth),
                                                   Y = log(policy_comparison_combine_sf$concentration_plus_median),
                                                   policy = policy_comparison_combine_sf$policy,
                                                   threshold_val = log(5))

test_our_method <- calculate_acc_mcc_two_sided_f1(T = log(policy_comparison_combine_sf_test$WellDepth),
                                                  Y = log(policy_comparison_combine_sf_test$concentration_plus_median),
                                                  policy = policy_comparison_combine_sf_test$policy,
                                                  threshold_val = log(5))

# assemble into one table: rows = metric, columns = method x train/test
metrics_comparison <- data.frame(
  Metric = c("Two-sided F1", "MCC", "Accuracy"),
  Training_Nonspatial = c(train_nonspatial$two_sided_f1, train_nonspatial$mcc, train_nonspatial$acc),
  Training_OurMethod = c(train_our_method$two_sided_f1, train_our_method$mcc, train_our_method$acc),
  Test_Nonspatial = c(test_nonspatial$two_sided_f1, test_nonspatial$mcc, test_nonspatial$acc),
  Test_OurMethod = c(test_our_method$two_sided_f1, test_our_method$mcc, test_our_method$acc)
)
metrics_comparison[, -1] <- round(metrics_comparison[, -1], 3)

print(metrics_comparison)

