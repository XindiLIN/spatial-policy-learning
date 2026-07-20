library(tidyverse)

source('functions/miscellaneous.r')

estimated_policy_lst <- readRDS('simulation/output/simulation_nonparametric_estimated_policy.rds')
data_test_lst <- readRDS('simulation/output/simulation_nonparametric_data_test.rds')

threshold_quantiles <- c(0.4, 0.5, 0.6)

# boxplot of the absolute error

estimated_policy_all <- bind_rows(estimated_policy_lst)

df_diff <- estimated_policy_all %>%
  # Pivot wide to get columns: direct, indirect, indirect_nonsp, optimal
  pivot_wider(names_from = method, values_from = value) %>%

  # Calculate squared differences
  mutate(
    diff_direct = abs(direct - optimal),
    diff_indirect = abs(indirect - optimal),
    diff_indirect_nonsp = abs(indirect_nonsp - optimal)
  ) %>%
  dplyr::select(threshold_quantile, starts_with("diff_")) %>%
  pivot_longer(
    cols = starts_with("diff_"),
    names_to = "method",
    values_to = "abs_error"
  ) %>%
  mutate(method = gsub("diff_", "", method))

ggplot(df_diff) +
  geom_boxplot(aes(x = threshold_quantile, y = abs_error, fill = method,
                   group = interaction(threshold_quantile, method))) +
  coord_cartesian(ylim = c(0, 6)) +
  scale_x_continuous(
    name = "Threshold Quantile", # Change the X-axis title here
    breaks = c(0.4,0.5,0.6))+
  scale_y_continuous(
    name = "ABsolute Error", # Change the X-axis title here
    breaks = c(0,2,6),
    limits = c(0,6))+
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(labels = c("direct" = "Direct", "indirect" = "Indirect", "indirect_nonsp" = "Non-spatial Indirect"),
                    values = c("direct" = "#F8766D",
                               "indirect" = "#00BFC4",
                               "indirect_nonsp" = "#FDBF6F"),
                    name = "Estimation Method")

# draw the binary classification metrics, for every method and every quantile, we get a metrics

data_test_all <- bind_rows(data_test_lst)



estimated_policy_all_wider <- estimated_policy_all %>%
  # Pivot wide to get columns: direct, indirect, indirect_nonsp, optimal
  pivot_wider(names_from = method, values_from = value)

estimated_policy_all_wider_combine <- left_join(estimated_policy_all_wider, data_test_all, by = c('seed','threshold_quantile','test_obs_id'))

performance_metrics <- expand.grid(method = c("direct", "indirect", "indirect_nonsp","optimal"), threshold_quantile = threshold_quantiles)
performance_metrics$mcc <- NA
performance_metrics$acc <- NA
performance_metrics$two_sided_f1 <- NA

for(i in 1:nrow(performance_metrics)){
  method <- performance_metrics[i,'method']
  threshold_quantile <- performance_metrics[i,'threshold_quantile']
  df_tmp <- estimated_policy_all_wider_combine[(estimated_policy_all_wider_combine$threshold_quantile == threshold_quantile), ]
  result <- calculate_acc_mcc_two_sided_f1(T = df_tmp[,'T'],
                                           Y = df_tmp[,'Y'],
                                           policy = df_tmp[,method],
                                           threshold_val = df_tmp[,'threshold_val'])
  performance_metrics[i,'mcc'] <- result$mcc
  performance_metrics[i,'acc'] <- result$acc
  performance_metrics[i,'two_sided_f1'] <- result$two_sided_f1
}

# 2. Reshape data to long format
df_long <- performance_metrics %>%
  pivot_longer(
    cols = c(mcc, acc, two_sided_f1),
    names_to = "metric",
    values_to = "value"
  )

# 3. Plot with facet_wrap

metric_labs <- c(
  acc = "Accuracy",
  two_sided_f1 = "Aggregated F1",
  mcc = "MCC"
)

ggplot(df_long, aes(x = threshold_quantile, y = value, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = metric_labs)) +  # 'scales = "free_y"' allows y-axis to vary per plot
  labs(title = "Performance Metrics by Threshold Quantile",
       x = "Threshold Quantile") +
  scale_color_manual(labels = c("direct" = "Direct", "indirect" = "Indirect", "indirect_nonsp" = "Non-spatial Indirect", "optimal" = "True"),
                   values = c("direct" = "#F8766D",
                              "indirect" = "#00BFC4",
                              "indirect_nonsp" = "#FDBF6F",
                              "optimal" = "purple"),
                   name = "Estimation Method") +
  theme_minimal()
