# Descriptive findings

The analysis-ready file contains 1,061 unique line-environment phenotype rows across 107 lines and 15 environments. Per-environment line counts range from 37 to 107; two Philippines environments include all 107 lines, while the Hyderabad 2018 wet trial includes 37. The distribution of observed line-environment records is shown in `results/figures/phenotype_coverage.png`; `reports/line_index_key.csv` maps the plotted line indices to source IDs.

All 882 retained SNPs were polymorphic in the panel and had no missing calls after the source HapMap was parsed. Genomic PCA explains 11.1% of marker variation on PC1 and 6.3% on PC2 (17.4% combined). The first eight-cluster cut has group sizes 11, 23, 3, 21, 30, 15, 2, and 2; these are algorithmic sampling strata, not asserted biological populations. Cluster-count sensitivity was run at 4, 6, 8, 10, and 12 groups; see `reports/cluster_count_sensitivity.md`. The outcomes vary with the cut.

Published adjusted yield distributions differ substantially across environments. Therefore, pooled raw-scale accuracy can be dominated by environment means. The validation reports compare environment-specific error and average across environments equally. See `results/figures/yield_by_environment.png` and `reports/yield_by_environment.csv` for the distribution and summary per trial.

Figures and tables are reproducible with `R/06_descriptive_eda.R`. The marker PCA and heatmap are diagnostics; they do not imply biological groups or explain causal mechanisms.

