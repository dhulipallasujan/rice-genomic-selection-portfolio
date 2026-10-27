# Paired sensitivity to duplicate line-environment-trait handling.
# Bootstrap resampling is at the environment level; the panel has only 15 trials.
suppressPackageStartupMessages(library(data.table))
set.seed(20261030)
mean_file <- "results/nested_multienv_cv_environment_summary.csv"
exclude_file <- "results/duplicate_exclusion_cv_environment_summary.csv"
stopifnot(file.exists(mean_file), file.exists(exclude_file))
a <- fread(mean_file); b <- fread(exclude_file)
metrics <- c("mean_correlation", "mean_rmse", "mean_baseline_rmse",
             "mean_rmse_difference", "mean_top20_overlap", "mean_selection_regret")
setnames(a, metrics, paste0(metrics, "_meanpolicy"))
setnames(b, metrics, paste0(metrics, "_excludepolicy"))
d <- merge(a, b, by = c("EnvID", "strategy", "fraction"), all = FALSE)
if (!nrow(d)) stop("No shared environment-strategy-fraction groups to compare.")
for (m in metrics) d[, (paste0(m, "_delta")) := get(paste0(m, "_excludepolicy")) - get(paste0(m, "_meanpolicy"))]
delta_names <- paste0(metrics, "_delta")
summary_rows <- list(); z <- 0L
for (s in unique(d$strategy)) for (f in sort(unique(d$fraction))) {
  x <- d[strategy == s & fraction == f]
  for (m in delta_names) {
    vals <- x[[m]]; vals <- vals[is.finite(vals)]
    boot <- replicate(10000L, mean(sample(vals, length(vals), replace = TRUE)))
    z <- z + 1L
    summary_rows[[z]] <- data.table(strategy = s, fraction = f, metric = sub("_delta$", "", m),
      environments = length(vals), mean_delta = mean(vals), ci_low = unname(quantile(boot, .025)),
      ci_high = unname(quantile(boot, .975)))
  }
}
fwrite(d, "results/duplicate_sensitivity_by_environment.csv")
fwrite(rbindlist(summary_rows), "reports/duplicate_sensitivity.csv")
message("Wrote paired duplicate-policy contrasts for ", uniqueN(d$EnvID), " environments.")
