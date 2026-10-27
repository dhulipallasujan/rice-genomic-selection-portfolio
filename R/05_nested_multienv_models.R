# Nested model/regularization selection for the global held-out-line design.
# All test lines are excluded from all training environments. Model type and
# ridge parameter are selected using inner grouped folds from training lines only.
suppressPackageStartupMessages(library(data.table))
set.seed(20261029)
n_outer <- 10L; n_inner <- 3L
fractions <- c(.25, .50, .75)
strategies <- c("random", "stratified_proportional")
models <- c("G+E", "G+E+GE")
lambdas <- c(.01, .1, 1)
cluster_k <- getOption("rice.cluster.k", 8L)
if (length(cluster_k) != 1L || cluster_k < 2L || cluster_k > 20L) stop("rice.cluster.k must be an integer from 2 to 20")
Gtab <- fread("data/processed/genotype_dosage.csv", check.names = FALSE)
duplicate_policy <- getOption("rice.duplicate.policy", "mean")
if (!duplicate_policy %in% c("mean", "exclude")) stop("rice.duplicate.policy must be mean or exclude")
phenotype_file <- if (duplicate_policy == "exclude") "data/processed/phenotypes_wide_exclude_dupkeys.csv" else "data/processed/phenotypes_wide.csv"
P <- fread(phenotype_file)
line_ids <- sub("^L", "", names(Gtab)[-1])
M <- t(as.matrix(Gtab[, -1])); rownames(M) <- line_ids
for (j in seq_len(ncol(M))) if (anyNA(M[, j])) M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
sdm <- apply(M, 2, sd); M <- M[, is.finite(sdm) & sdm > 0, drop = FALSE]
M <- scale(M); K <- tcrossprod(M) / ncol(M); dimnames(K) <- list(line_ids, line_ids)
D <- 1 - K; D[D < 0] <- 0; diag(D) <- 0
cl <- cutree(hclust(as.dist(D), method = "average"), k = cluster_k)
P[, LineID := as.character(LineID)]; P[, EnvID := as.character(EnvID)]
P <- P[!is.na(YLD) & LineID %in% line_ids]
env_levels <- sort(unique(P$EnvID)); P[, e := match(EnvID, env_levels)]; P[, g := match(LineID, line_ids)]

choose_stratified <- function(pool, n) {
  groups <- sort(unique(cl[pool]))
  sizes <- vapply(groups, function(g) sum(cl[pool] == g), integer(1))
  raw <- n * sizes / sum(sizes); quota <- floor(raw)
  rem <- n - sum(quota)
  if (rem > 0L) {
    ix <- order(raw - quota, decreasing = TRUE)[seq_len(rem)]
    quota[ix] <- quota[ix] + 1L
  }
  unlist(lapply(seq_along(groups), function(i) {
    cand <- pool[cl[pool] == groups[i]]
    cand[sample.int(length(cand), quota[i])]
  }), use.names = FALSE)
}
predict_kernel <- function(train, valid, model, lambda) {
  kg <- K[train$g, train$g, drop = FALSE]
  ke <- outer(train$e, train$e, "==") * 1
  kt <- kg + ke
  kx <- K[valid$g, train$g, drop = FALSE] + outer(valid$e, train$e, "==") * 1
  if (model == "G+E+GE") {
    kt <- kt + kg * ke
    kx <- kx + K[valid$g, train$g, drop = FALSE] * outer(valid$e, train$e, "==")
  }
  mu <- mean(train$YLD)
  as.numeric(mu + kx %*% solve(kt + diag(lambda, nrow(kt)), train$YLD - mu))
}
rmse_equal_env <- function(dat, pred) {
  tmp <- copy(dat); tmp[, pred__ := pred]
  mean(tmp[, sqrt(mean((YLD - pred__)^2)), by = EnvID]$V1)
}
out <- list(); z <- 0L
for (outer in seq_len(n_outer)) {
  all_ids <- sort(unique(P$LineID)); ntest <- max(2L, floor(.20 * length(all_ids)))
  test_ids <- sample(all_ids, ntest); pool <- setdiff(all_ids, test_ids)
  valid_outer <- P[LineID %in% test_ids]
  for (frac in fractions) for (strategy in strategies) {
    ntrain <- max(12L, floor(length(pool) * frac))
    selected <- if (strategy == "random") sample(pool, ntrain) else choose_stratified(pool, ntrain)
    # Inner folds are grouped by line, never by individual line-environment rows.
    inner_fold <- sample(rep(seq_len(n_inner), length.out = length(selected)))
    score_rows <- list(); q <- 0L
    for (model in models) for (lambda in lambdas) {
      fold_scores <- numeric()
      for (fold in seq_len(n_inner)) {
        inner_valid_ids <- selected[inner_fold == fold]
        inner_train_ids <- setdiff(selected, inner_valid_ids)
        tr <- P[LineID %in% inner_train_ids]; va <- P[LineID %in% inner_valid_ids]
        if (!nrow(tr) || !nrow(va)) next
        pv <- predict_kernel(tr, va, model, lambda)
        fold_scores <- c(fold_scores, rmse_equal_env(va, pv))
      }
      q <- q + 1L
      score_rows[[q]] <- data.table(model = model, lambda = lambda,
                                    inner_mean_rmse = mean(fold_scores), n_inner_folds = length(fold_scores))
    }
    tuning <- rbindlist(score_rows)
    best <- tuning[order(inner_mean_rmse)][1]
    tr <- P[LineID %in% selected]
    ptest <- predict_kernel(tr, valid_outer, best$model, best$lambda)
    eval <- copy(valid_outer); eval[, prediction := ptest]
    envmeans <- tr[, .(envmean = mean(YLD)), by = EnvID]
    eval <- merge(eval, envmeans, by = "EnvID", all.x = TRUE, sort = FALSE)
    grand <- mean(tr$YLD); eval[is.na(envmean), envmean := grand]
    ev <- eval[, {
      k <- max(1L, ceiling(.20 * .N)); chosen <- order(prediction, decreasing = TRUE)[seq_len(k)]; oracle <- order(YLD, decreasing = TRUE)[seq_len(k)]
      .(n_test = .N,
        correlation = if (.N > 1L && sd(prediction) > 0 && sd(YLD) > 0) cor(prediction, YLD) else NA_real_,
        rmse = sqrt(mean((prediction - YLD)^2)), baseline_rmse = sqrt(mean((envmean - YLD)^2)),
        rmse_difference = sqrt(mean((prediction - YLD)^2)) - sqrt(mean((envmean - YLD)^2)),
        top20_overlap = length(intersect(chosen, oracle)) / k,
        selection_regret = mean(YLD[oracle]) - mean(YLD[chosen]))
    }, by = EnvID]
    ev[, `:=`(outer_repeat = outer, strategy = strategy, fraction = frac, n_train_lines = ntrain,
              selected_model = best$model, selected_lambda = best$lambda,
              inner_tuning_rmse = best$inner_mean_rmse)]
    z <- z + 1L; out[[z]] <- ev
  }
}
detail <- rbindlist(out)
prefix <- if (duplicate_policy == "exclude") "duplicate_exclusion" else if (cluster_k == 8L) "nested_multienv" else sprintf("cluster_k%d", cluster_k)
fwrite(detail, sprintf("results/%s_cv_by_environment.csv", prefix))
model_choices <- unique(detail[, .(outer_repeat, strategy, fraction, selected_model, selected_lambda)])
model_choices <- model_choices[, .N, by = .(strategy, fraction, selected_model, selected_lambda)]
fwrite(model_choices, sprintf("results/%s_model_selection_counts.csv", prefix))
env_summary <- detail[, .(outer_repeats = .N,
  mean_correlation = mean(correlation, na.rm = TRUE), mean_rmse = mean(rmse),
  mean_baseline_rmse = mean(baseline_rmse), mean_rmse_difference = mean(rmse_difference),
  mean_top20_overlap = mean(top20_overlap), mean_selection_regret = mean(selection_regret)),
  by = .(EnvID, strategy, fraction)]
fwrite(env_summary, sprintf("results/%s_cv_environment_summary.csv", prefix))
overall <- env_summary[, .(environments = .N,
  mean_correlation = mean(mean_correlation, na.rm = TRUE), mean_rmse = mean(mean_rmse),
  mean_baseline_rmse = mean(mean_baseline_rmse), mean_rmse_difference = mean(mean_rmse_difference),
  sd_rmse_difference_across_environments = sd(mean_rmse_difference),
  mean_top20_overlap = mean(mean_top20_overlap), mean_selection_regret = mean(mean_selection_regret)),
  by = .(strategy, fraction)]
fwrite(overall, sprintf("results/%s_cv_summary.csv", prefix))
message("Completed nested selection for ", nrow(detail), " environment-fold groups using duplicate policy '", duplicate_policy, "' and k = ", cluster_k, "; model and penalty selected using training lines only.")
