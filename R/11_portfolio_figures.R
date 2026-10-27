# Portfolio figure: model error versus the environment-mean benchmark by cluster cut.
suppressPackageStartupMessages(library(data.table))
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
d <- fread("reports/cluster_count_sensitivity.csv")
png("results/figures/cluster_sensitivity_rmse.png", width = 1800, height = 1200, res = 200)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 1), oma = c(0, 0, 2, 0))
cols <- c(random = "#245B70", stratified_proportional = "#C05A35")
for (f in sort(unique(d$fraction))) {
  x <- d[fraction == f]
  yr <- range(x$mean_rmse_difference, 0)
  plot(NA, xlim = range(x$cluster_k), ylim = yr, xaxt = "n",
       xlab = "Number of genotype clusters", ylab = "RMSE - environment-mean RMSE",
       main = paste0("Training budget: ", 100*f, "%"))
  axis(1, at = sort(unique(x$cluster_k)))
  abline(h = 0, lty = 2, col = "#777777")
  for (s in names(cols)) {
    y <- x[strategy == s][order(cluster_k)]
    lines(y$cluster_k, y$mean_rmse_difference, type = "b", pch = if (s == "random") 16 else 17,
          lwd = 2, col = cols[s])
  }
  if (f == min(d$fraction))
    legend("bottomright", legend = c("Random", "Proportional stratified"),
           col = cols, pch = c(16, 17), lty = 1, bty = "white", cex = .8)
}
mtext("The apparent value of stratification changes with the cluster cut", outer = TRUE, cex = 1.2, font = 2)
dev.off()
message("Wrote portfolio cluster-sensitivity figure.")
