# Repeated within-environment held-out-line baseline for grain yield.
# Question: can marker relationships predict a line's yield within an already
# represented environment? Each environment is fit separately, so a test line
# is excluded from that environment's training set; splits are not global
# across environments and this is not a multi-environment model.
# This is not new-environment validation and does not reproduce BGLR models.
suppressPackageStartupMessages(library(data.table))
set.seed(20261027)
n_repeats <- 30L
dir.create("results", showWarnings = FALSE)
Gtab <- fread("data/processed/genotype_dosage.csv", check.names = FALSE)
P <- fread("data/processed/phenotypes_wide.csv")
line_ids <- sub("^L", "", names(Gtab)[-1])
M <- t(as.matrix(Gtab[, -1]))
rownames(M) <- line_ids
if (anyNA(M)) {
  for (j in seq_len(ncol(M))) M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
}
# Center and scale marker columns, dropping any constant loci.
sdm <- apply(M, 2, sd); M <- M[, is.finite(sdm) & sdm > 0, drop = FALSE]
M <- scale(M)
K <- tcrossprod(M) / ncol(M)
rownames(K) <- colnames(K) <- line_ids
P[, LineID := as.character(LineID)]
P[, EnvID := as.character(EnvID)]
P <- P[!is.na(YLD)]
envs <- sort(unique(P$EnvID))
fractions <- c(.25, .50, .75)
strategies <- c("random", "stratified_proportional")
out <- list(); z <- 0L

# Cluster labels come only from genotype similarities, never phenotypes.
D <- 1 - K
D[D < 0] <- 0
diag(D) <- 0
hc <- hclust(as.dist(D), method = "average")
for (rep in seq_len(n_repeats)) {
  for (env in envs) {
    obs <- P[EnvID == env]
    ids <- intersect(obs$LineID, line_ids)
    if (length(ids) < 12L) next
    ntest <- max(2L, floor(length(ids) * .20))
    test_ids <- sample(ids, ntest)
    pool <- setdiff(ids, test_ids)
    yenv <- setNames(obs$YLD, obs$LineID)
    base <- mean(yenv[pool], na.rm = TRUE)
    for (frac in fractions) for (strategy in strategies) {
      ntrain <- max(5L, min(length(pool), floor(length(pool) * frac)))
      if (strategy == "random") train_ids <- sample(pool, ntrain)
      else {
        cl <- cutree(hc, k = min(8L, max(2L, floor(length(pool) / 4L))))
        groups <- sort(unique(cl[pool]))
        sizes <- vapply(groups, function(g) sum(cl[pool] == g), integer(1))
        quota_raw <- ntrain * sizes / sum(sizes); quota <- floor(quota_raw)
        remainder <- ntrain - sum(quota)
        if (remainder > 0L) quota[order(quota_raw - quota, decreasing = TRUE)[seq_len(remainder)]] <- quota[order(quota_raw - quota, decreasing = TRUE)[seq_len(remainder)]] + 1L
        train_ids <- unlist(lapply(seq_along(groups), function(i) {
          candidates <- pool[cl[pool] == groups[i]]
          candidates[sample.int(length(candidates), quota[i])]
        }), use.names = FALSE)
      }
      ytr <- yenv[train_ids]; mu <- mean(ytr, na.rm = TRUE)
      lam <- .1
      pred <- as.numeric(mu + K[test_ids, train_ids, drop = FALSE] %*%
                           solve(K[train_ids, train_ids, drop = FALSE] + diag(lam, length(train_ids)), ytr - mu))
      names(pred) <- test_ids
      truth <- yenv[test_ids]
      envpred <- rep(mu, length(truth))
      pickn <- max(1L, ceiling(length(truth) * .20))
      rankpred <- names(sort(pred, decreasing = TRUE))[seq_len(pickn)]
      ranktrue <- names(sort(truth, decreasing = TRUE))[seq_len(pickn)]
      z <- z + 1L
      out[[z]] <- data.table(replicate = rep, environment = env, n_observed = length(ids), n_train = ntrain,
        fraction = frac, strategy = strategy,
        cor = if (sd(pred) > 0 && sd(truth) > 0) cor(pred, truth) else NA_real_,
        rmse = sqrt(mean((pred - truth)^2)),
        baseline_rmse = sqrt(mean((envpred - truth)^2)),
        top20_overlap = length(intersect(rankpred, ranktrue)) / pickn,
        selected_yield = mean(truth[rankpred]), oracle_top_yield = mean(truth[ranktrue]),
        selected_regret = mean(truth[ranktrue]) - mean(truth[rankpred]))
    }
  }
}
res <- rbindlist(out)
fwrite(res, "results/baseline_cv_replicates.csv")
env_summary <- res[, .(repeats = .N,
  mean_correlation = mean(cor, na.rm = TRUE), mean_rmse = mean(rmse), mean_baseline_rmse = mean(baseline_rmse),
  mean_rmse_difference = mean(rmse - baseline_rmse), mean_top20_overlap = mean(top20_overlap),
  mean_selected_yield = mean(selected_yield), mean_selected_regret = mean(selected_regret)),
  by = .(environment, strategy, fraction)]
fwrite(env_summary, "results/baseline_cv_by_environment.csv")
# Give each environment equal weight in the overall summary; do not pool raw
# RMSE or yields across trials with different scales.
summary <- env_summary[, .(environments = .N,
  mean_correlation = mean(mean_correlation, na.rm = TRUE), mean_rmse = mean(mean_rmse),
  mean_baseline_rmse = mean(mean_baseline_rmse), mean_rmse_difference = mean(mean_rmse_difference),
  sd_rmse_difference_across_environments = sd(mean_rmse_difference),
  mean_top20_overlap = mean(mean_top20_overlap), mean_selected_regret = mean(mean_selected_regret)),
  by = .(strategy, fraction)]
fwrite(summary, "results/baseline_cv_summary.csv")
message("Completed ", nrow(res), " model configurations across ", length(envs), " environments and ", n_repeats, " repeats.")
