# Paired environment-level uncertainty for the nested global validation.
# Avoids treating line-environment observations as independent replicates.
suppressPackageStartupMessages(library(data.table))
set.seed(20261030)
d <- fread("results/nested_multienv_cv_by_environment.csv")
envmean <- d[, .(mean_correlation = mean(correlation, na.rm = TRUE),
  mean_rmse = mean(rmse), mean_top20_overlap = mean(top20_overlap),
  mean_selection_regret = mean(selection_regret)), by = .(EnvID, fraction, strategy)]
metrics <- c("mean_correlation", "mean_rmse", "mean_top20_overlap", "mean_selection_regret")
boot_n <- 10000L; out <- list(); z <- 0L
for (frac in sort(unique(envmean$fraction))) {
  a <- dcast(envmean[fraction == frac], EnvID ~ strategy, value.var = metrics)
  for (metric in metrics) {
    rand_col <- paste0(metric, "_random")
    strat_col <- paste0(metric, "_stratified_proportional")
    delta <- a[[strat_col]] - a[[rand_col]]
    delta <- delta[is.finite(delta)]
    if (length(delta) < 3L) next
    b <- replicate(boot_n, mean(sample(delta, length(delta), replace = TRUE)))
    z <- z + 1L
    out[[z]] <- data.table(fraction = frac, metric = metric, n_environments = length(delta),
      stratified_minus_random = mean(delta), ci95_low = unname(quantile(b, .025)),
      ci95_high = unname(quantile(b, .975)))
  }
}
contrasts <- rbindlist(out)
fwrite(contrasts, "reports/paired_environment_strategy_contrasts.csv")
message("Wrote paired differences with percentile intervals from 10,000 environment-level bootstrap samples.")
