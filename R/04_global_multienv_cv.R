# Global held-out-line multi-environment prediction and training-sample comparison.
# Test lines are removed from every environment's training records for a repeat.
# Linear covariance components: environment (E), genomic main effect (G), and
# genomic-by-environment interaction (GxE). Sampling uses genotype data only.
suppressPackageStartupMessages(library(data.table))
set.seed(20261028)
n_repeats <- 10L
fractions <- c(.25, .50, .75)
strategies <- c("random", "stratified_proportional")
Gtab <- fread("data/processed/genotype_dosage.csv", check.names = FALSE)
P <- fread("data/processed/phenotypes_wide.csv")
line_ids <- sub("^L", "", names(Gtab)[-1])
M <- t(as.matrix(Gtab[, -1])); rownames(M) <- line_ids
for (j in seq_len(ncol(M))) if (anyNA(M[, j])) M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
sdm <- apply(M, 2, sd); M <- M[, is.finite(sdm) & sdm > 0, drop = FALSE]
M <- scale(M)
K <- tcrossprod(M) / ncol(M); dimnames(K) <- list(line_ids, line_ids)
D <- 1 - K; D[D < 0] <- 0; diag(D) <- 0
hc <- hclust(as.dist(D), method = "average")
cl <- cutree(hc, k = 8L)
P[, LineID := as.character(LineID)]; P[, EnvID := as.character(EnvID)]
P <- P[!is.na(YLD) & LineID %in% line_ids]
env_levels <- sort(unique(P$EnvID))
P[, e := match(EnvID, env_levels)]
P[, g := match(LineID, line_ids)]
P[, obs_id := .I]
lambda <- .10
out <- list(); z <- 0L

choose_stratified <- function(pool, n, cluster) {
  groups <- sort(unique(cluster[pool]))
  sizes <- vapply(groups, function(g) sum(cluster[pool] == g), integer(1))
  raw <- n * sizes / sum(sizes); quota <- floor(raw)
  rem <- n - sum(quota)
  if (rem > 0L) {
    ix <- order(raw - quota, decreasing = TRUE)[seq_len(rem)]
    quota[ix] <- quota[ix] + 1L
  }
  unlist(lapply(seq_along(groups), function(i) {
    candidates <- pool[cluster[pool] == groups[i]]
    candidates[sample.int(length(candidates), quota[i])]
  }), use.names = FALSE)
}

for (replicate in seq_len(n_repeats)) {
  observed_lines <- sort(unique(P$LineID))
  n_test <- max(2L, floor(.20 * length(observed_lines)))
  test_ids <- sample(observed_lines, n_test)
  candidate_ids <- setdiff(observed_lines, test_ids)
  test_rows <- P[LineID %in% test_ids]
  for (fraction in fractions) for (strategy in strategies) {
    n_train <- max(8L, min(length(candidate_ids), floor(length(candidate_ids) * fraction)))
    train_ids <- if (strategy == "random") sample(candidate_ids, n_train) else choose_stratified(candidate_ids, n_train, cl)
    train_rows <- P[LineID %in% train_ids]
    stopifnot(!length(intersect(train_rows$LineID, test_rows$LineID)))

    # Observation-level covariance matrices, with training-only centering.
    Kg <- K[train_rows$g, train_rows$g, drop = FALSE]
    Ke <- outer(train_rows$e, train_rows$e, "==") * 1
    Kge <- Kg * Ke
    Ktrain <- Kg + Ke + Kge
    Kg_cross <- K[test_rows$g, train_rows$g, drop = FALSE]
    Ke_cross <- outer(test_rows$e, train_rows$e, "==") * 1
    Kcross <- Kg_cross + Ke_cross + Kg_cross * Ke_cross
    grand_mean <- mean(train_rows$YLD)
    alpha <- solve(Ktrain + diag(lambda, nrow(Ktrain)), train_rows$YLD - grand_mean)
    pred <- as.numeric(grand_mean + Kcross %*% alpha)
    test_eval <- copy(test_rows)
    test_eval[, prediction := pred]
    envmeans <- train_rows[, .(env_mean = mean(YLD)), by = EnvID]
    test_eval <- merge(test_eval, envmeans, by = "EnvID", all.x = TRUE, sort = FALSE)
    test_eval[is.na(env_mean), env_mean := grand_mean]
    test_eval[, baseline_prediction := env_mean]

    env_metrics <- test_eval[, {
      k <- max(1L, ceiling(.20 * .N))
      chosen <- order(prediction, decreasing = TRUE)[seq_len(k)]
      oracle <- order(YLD, decreasing = TRUE)[seq_len(k)]
      .(n_test = .N,
        correlation = if (.N > 1L && sd(prediction) > 0 && sd(YLD) > 0) cor(prediction, YLD) else NA_real_,
        rmse = sqrt(mean((prediction - YLD)^2)), baseline_rmse = sqrt(mean((baseline_prediction - YLD)^2)),
        rmse_difference = sqrt(mean((prediction - YLD)^2)) - sqrt(mean((baseline_prediction - YLD)^2)),
        top20_overlap = length(intersect(chosen, oracle)) / k,
        selected_yield = mean(YLD[chosen]), oracle_top_yield = mean(YLD[oracle]),
        selection_regret = mean(YLD[oracle]) - mean(YLD[chosen]))
    }, by = EnvID]
    z <- z + 1L
    out[[z]] <- env_metrics[, `:=`(replicate = replicate, strategy = strategy, fraction = fraction,
                                      n_train_lines = n_train, n_train_records = nrow(train_rows))]
  }
}
detail <- rbindlist(out)
fwrite(detail, "results/global_multienv_cv_by_environment.csv")
env_summary <- detail[, .(repeats = .N,
  mean_correlation = mean(correlation, na.rm = TRUE), mean_rmse = mean(rmse), mean_baseline_rmse = mean(baseline_rmse),
  mean_rmse_difference = mean(rmse_difference), sd_rmse_difference = sd(rmse_difference),
  mean_top20_overlap = mean(top20_overlap), mean_selection_regret = mean(selection_regret)),
  by = .(EnvID, strategy, fraction)]
fwrite(env_summary, "results/global_multienv_cv_environment_summary.csv")
overall <- env_summary[, .(environments = .N,
  mean_correlation = mean(mean_correlation, na.rm = TRUE), mean_rmse = mean(mean_rmse),
  mean_baseline_rmse = mean(mean_baseline_rmse), mean_rmse_difference = mean(mean_rmse_difference),
  sd_rmse_difference_across_environments = sd(mean_rmse_difference),
  mean_top20_overlap = mean(mean_top20_overlap), mean_selection_regret = mean(mean_selection_regret)),
  by = .(strategy, fraction)]
fwrite(overall, "results/global_multienv_cv_summary.csv")
message("Completed ", nrow(detail), " environment-by-fold result groups across ", n_repeats,
        " global line-held-out repeats; test lines were excluded from all training environments.")
