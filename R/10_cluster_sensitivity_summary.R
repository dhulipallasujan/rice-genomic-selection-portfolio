# Summarize pre-specified cluster-count sensitivity runs from R/05.
suppressPackageStartupMessages(library(data.table))
ks <- c(4L, 6L, 8L, 10L, 12L)
files <- c("results/nested_multienv_cv_summary.csv",
           sprintf("results/cluster_k%d_cv_summary.csv", ks[ks != 8L]))
labels <- c(8L, ks[ks != 8L])
if (!all(file.exists(files))) stop("Run R/05_nested_multienv_models.R for k = 4, 6, 10, and 12 first.")
res <- rbindlist(lapply(seq_along(files), function(i) {
  x <- fread(files[i]); x[, cluster_k := labels[i]]; x
}), use.names = TRUE, fill = TRUE)
setcolorder(res, c("cluster_k", setdiff(names(res), "cluster_k")))
fwrite(res, "reports/cluster_count_sensitivity.csv")
env_files <- c("results/nested_multienv_cv_environment_summary.csv",
               sprintf("results/cluster_k%d_cv_environment_summary.csv", ks[ks != 8L]))
env <- rbindlist(lapply(seq_along(env_files), function(i) {
  x <- fread(env_files[i]); x[, cluster_k := labels[i]]; x
}), use.names = TRUE, fill = TRUE)
metrics <- c("mean_correlation", "mean_rmse", "mean_top20_overlap", "mean_selection_regret")
set.seed(20261031); nboot <- 10000L; contrasts <- list(); z <- 0L
for (k in labels) for (f in sort(unique(env$fraction))) {
  x <- env[cluster_k == k & fraction == f]
  for (m in metrics) {
    w <- dcast(x, EnvID ~ strategy, value.var = m)
    delta <- w$stratified_proportional - w$random
    delta <- delta[is.finite(delta)]
    boot <- replicate(nboot, mean(sample(delta, length(delta), replace = TRUE)))
    z <- z + 1L
    contrasts[[z]] <- data.table(cluster_k = k, fraction = f, metric = m,
      environments = length(delta), stratified_minus_random = mean(delta),
      ci95_low = unname(quantile(boot, .025)), ci95_high = unname(quantile(boot, .975)))
  }
}
fwrite(rbindlist(contrasts), "reports/cluster_count_strategy_contrasts.csv")
message("Wrote nested summaries and paired environment bootstrap contrasts for k = 4, 6, 8, 10, and 12.")
