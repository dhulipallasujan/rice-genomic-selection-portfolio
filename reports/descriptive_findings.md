# Descriptive findings

The analysis-ready file contains 1,061 unique line-environment phenotype rows across 107 lines and 15 environments. Per-environment line counts range from 37 to 107, while observed yield counts range from 31 to 107 because some rows have missing yield values. The figure `results/figures/phenotype_coverage.png` now counts only non-missing yield values and labels each bar with its numerator and the full panel size. Full trial names and counts are in `reports/yield_by_environment.csv`.

All 882 retained SNPs were polymorphic in the panel and had no missing calls after the source HapMap was parsed. Genomic PCA explains 11.1% of marker variation on PC1 and 6.3% on PC2 (17.4% combined). The first eight-cluster cut has group sizes 11, 23, 3, 21, 30, 15, 2, and 2; these are algorithmic sampling strata, not asserted biological populations. Cluster-count sensitivity was run at 4, 6, 8, 10, and 12 groups; see `reports/cluster_count_sensitivity.md`. The outcomes vary with the cut.

Published adjusted yield distributions differ substantially across environments. Therefore, pooled raw-scale accuracy can be dominated by environment means. The validation reports compare environment-specific error and average across environments equally. See `results/figures/yield_by_environment.png` and `reports/yield_by_environment.csv` for the distribution and summary per trial.

Figures and tables are reproducible with `R/06_descriptive_eda.R`. The marker PCA and coverage summary are diagnostics; they do not imply biological groups or explain causal mechanisms.

