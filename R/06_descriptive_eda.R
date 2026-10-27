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

# Yield coverage by environment. Count non-missing outcomes, not merely rows
# in the phenotype table: some line-environment rows have no yield value.
env_order <- sort(unique(P$EnvID))
env_codes <- c("Bangladesh_Gazipur_2019_Wet" = "BD_GZ_19W", "Bangladesh_Nizmawna_2019_Wet" = "BD_NM_19W",
  "India_Cuttack_2019_Wet" = "IN_CU_19W", "India_Hyderabad_2018_Wet" = "IN_HY_18W",
  "India_Hyderabad_2019_Dry" = "IN_HY_19D", "India_Hyderabad_2019_Wet" = "IN_HY_19W",
  "India_Maruteru_2019_Wet" = "IN_MA_19W", "India_Raipur_2019_Wet" = "IN_RP_19W",
  "Kenya_Ahero_2019_Dry" = "KE_AH_19D", "Kenya_Mwea_2020_Wet" = "KE_MW_20W",
  "Mozambique_Chokwe_2020_Wet" = "MZ_CK_20W", "Mozambique_Maputo_2020_Wet" = "MZ_MP_20W",
  "Philippines_Los Banos_2019_Dry" = "PH_LB_19D", "Philippines_Los Banos_2019_Wet" = "PH_LB_19W",
  "Tanzania_Dakawa_2020_Wet" = "TZ_DK_20W")
yield_coverage <- P[, .(observed = sum(!is.na(YLD))), by = EnvID][match(env_order, EnvID)]
yield_coverage[, label := unname(env_codes[EnvID])]
png("results/figures/phenotype_coverage.png", width = 2200, height = 1350, res = 220)
par(mar = c(5, 9, 5, 4), xaxs = "i", yaxs = "i", las = 1,
    family = "sans", col.axis = "#344054", col.lab = "#344054")
ord <- order(yield_coverage$observed)
yy <- barplot(yield_coverage$observed[ord], horiz = TRUE,
              names.arg = yield_coverage$label[ord],
              xlim = c(0, 120), xaxt = "n", xlab = "Lines with an observed yield value (of 107)",
              col = ifelse(yield_coverage$observed[ord] == max(yield_coverage$observed),
                           "#147D76", "#E39A55"), border = NA, space = .28,
              main = "Yield data coverage varies across trials",
              cex.names = 1.05, cex.lab = 1.05, cex.main = 1.35)
axis(1, at = seq(0, 100, 20), labels = seq(0, 100, 20), col = "#98A2B3", col.ticks = "#98A2B3")
for (x in seq(0, 100, 20)) abline(v = x, col = "#E8EDF2", lwd = 1)
# Redraw bars over the subtle count grid so the values remain easy to compare.
barplot(yield_coverage$observed[ord], horiz = TRUE,
        names.arg = rep("", nrow(yield_coverage)), add = TRUE,
        col = ifelse(yield_coverage$observed[ord] == max(yield_coverage$observed),
                     "#147D76", "#E39A55"), border = NA, space = .28,
        axes = FALSE)
abline(v = max(yield_coverage$observed), lty = 2, lwd = 1.5, col = "#667085")
text(yield_coverage$observed[ord] + 2, yy,
     labels = sprintf("%d / 107", yield_coverage$observed[ord]),
     adj = 0, cex = .9, col = "#344054", xpd = NA)
mtext("Counts exclude missing yield values; full environment names are in reports/yield_by_environment.csv",
      side = 3, line = .6, cex = .9, col = "#667085")
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
