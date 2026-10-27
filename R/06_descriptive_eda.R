# Compact descriptive analysis and publication-ready diagnostics.
suppressPackageStartupMessages(library(data.table))
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
P <- fread("data/processed/phenotypes_wide.csv")
Gtab <- fread("data/processed/genotype_dosage.csv", check.names = FALSE)
line_ids <- sub("^L", "", names(Gtab)[-1])
M <- t(as.matrix(Gtab[, -1])); rownames(M) <- line_ids
for (j in seq_len(ncol(M))) if (anyNA(M[, j])) M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
sdm <- apply(M, 2, sd); M <- M[, is.finite(sdm) & sdm > 0, drop = FALSE]
Msc <- scale(M)
P[, LineID := as.character(LineID)]; P[, EnvID := as.character(EnvID)]

coverage <- dcast(P[, .(observed = 1L), by = .(LineID, EnvID)], LineID ~ EnvID, value.var = "observed", fill = 0)
fwrite(P[, .(n_lines = uniqueN(LineID), yield_records = sum(!is.na(YLD)),
             mean_yield = mean(YLD, na.rm = TRUE), sd_yield = sd(YLD, na.rm = TRUE),
             min_yield = min(YLD, na.rm = TRUE), max_yield = max(YLD, na.rm = TRUE)), by = EnvID],
       "reports/yield_by_environment.csv")
fwrite(P[, .(n_lines = uniqueN(LineID), n_environments = uniqueN(EnvID),
             missing_DTF = sum(is.na(DTF)), missing_HT = sum(is.na(HT)), missing_YLD = sum(is.na(YLD))), by = LineID],
       "reports/trait_coverage_by_line.csv")

# Genotype PCA and unsupervised groups used by the diversity sampling baseline.
pc <- prcomp(Msc, center = FALSE, scale. = FALSE, rank. = 5)
pcs <- data.table(LineID = rownames(pc$x), PC1 = pc$x[,1], PC2 = pc$x[,2])
K <- tcrossprod(Msc) / ncol(Msc); dimnames(K) <- list(line_ids, line_ids)
D <- 1-K; D[D < 0] <- 0; diag(D) <- 0
pcs[, cluster := factor(cutree(hclust(as.dist(D), method = "average"), k = 8L)[LineID])]
fwrite(pcs, "results/genotype_pca_scores.csv")
png("results/figures/genotype_pca.png", width = 1800, height = 1350, res = 220)
cols <- grDevices::hcl.colors(nlevels(pcs$cluster), "Dark 3")
plot(pcs$PC1, pcs$PC2, col = cols[pcs$cluster], pch = 19, cex = .8,
     xlab = sprintf("PC1 (%.1f%%)", 100*summary(pc)$importance[2,1]),
     ylab = sprintf("PC2 (%.1f%%)", 100*summary(pc)$importance[2,2]),
     main = "Rice elite panel: genome-wide marker PCA")
legend("topright", legend = paste("Cluster", levels(pcs$cluster)), col = cols, pch = 19, bty = "n", cex = .8)
dev.off()

# Unbalanced phenotype coverage: rows are environments, columns are lines.
line_order <- pcs$LineID
env_order <- sort(unique(P$EnvID))
env_codes <- c("Bangladesh_Gazipur_2019_Wet" = "BD_GZ_19W", "Bangladesh_Nizmawna_2019_Wet" = "BD_NM_19W",
  "India_Cuttack_2019_Wet" = "IN_CU_19W", "India_Hyderabad_2018_Wet" = "IN_HY_18W",
  "India_Hyderabad_2019_Dry" = "IN_HY_19D", "India_Hyderabad_2019_Wet" = "IN_HY_19W",
  "India_Maruteru_2019_Wet" = "IN_MA_19W", "India_Raipur_2019_Wet" = "IN_RP_19W",
  "Kenya_Ahero_2019_Dry" = "KE_AH_19D", "Kenya_Mwea_2020_Wet" = "KE_MW_20W",
  "Mozambique_Chokwe_2020_Wet" = "MZ_CK_20W", "Mozambique_Maputo_2020_Wet" = "MZ_MP_20W",
  "Philippines_Los Banos_2019_Dry" = "PH_LB_19D", "Philippines_Los Banos_2019_Wet" = "PH_LB_19W",
  "Tanzania_Dakawa_2020_Wet" = "TZ_DK_20W")
fwrite(data.table(line_index = seq_along(line_order), LineID = line_order), "reports/line_index_key.csv")
zmat <- matrix(0, nrow = length(line_order), ncol = length(env_order),
               dimnames = list(line_order, env_order))
for (i in seq_len(nrow(P))) zmat[P$LineID[i], P$EnvID[i]] <- 1
png("results/figures/phenotype_coverage.png", width = 1800, height = 1050, res = 220)
par(mar = c(5, 8, 5, 2))
image(seq_along(line_order), seq_along(env_order), zmat,
      col = c("#F2F4F7", "#147D76"), axes = FALSE, xlab = "",
      ylab = "", main = "Observed yield records by line and environment")
axis(2, at = seq_along(env_order), labels = unname(env_codes[env_order]), las = 2, cex.axis = .9)
axis(1, at = seq(1, length(line_order), length.out = 6), labels = round(seq(1, length(line_order), length.out = 6)), cex.axis = .8)
mtext("Genotype line index (see reports/line_index_key.csv)", side = 1, line = 3)
mtext("Dark = observed; pale = missing phenotype", side = 3, line = .4, cex = .85)
dev.off()

# Environment-wise yield distributions, with points revealing sample sizes.
png("results/figures/yield_by_environment.png", width = 2050, height = 1250, res = 220)
par(mar = c(11, 5, 4, 1))
boxplot(YLD ~ factor(EnvID, levels = env_order), data = P, las = 2,
        col = "#DDEFEA", border = "#147D76", ylab = "Published adjusted grain yield (source units)",
        xlab = "", main = "Grain-yield distribution differs across trials", outline = FALSE, cex.axis = .65)
for (j in seq_along(env_order)) {
  yy <- P[EnvID == env_order[j] & !is.na(YLD), YLD]
  points(jitter(rep(j, length(yy)), amount = .14), yy, pch = 16, cex = .35, col = grDevices::adjustcolor("#245B70", alpha.f = .38))
}
dev.off()
message("Wrote environment summaries and three diagnostic figures.")
