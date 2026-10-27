# Cluster-count sensitivity

## Purpose

The proportional stratified rule depends on how the marker relationship tree is cut into groups. The primary analysis uses eight groups, but this choice has no biological interpretation. To assess dependence on this choice, the nested global line-holdout pipeline was rerun with 4, 6, 10, and 12 groups. The runs used the same seeds, outer line splits, budget fractions, model grid, and training-only inner validation; only the genotype-cluster partition changed. Intervals resample the 15 environments.

## Finding

The sampling comparison is sensitive to cluster count. At k=8, proportional stratification had lower top-20% overlap and higher regret than random sampling at the 75% budget. With k=4, stratification improved 25% top-set overlap by 0.046 (95% interval 0.005 to 0.084) and reduced regret by 0.065 (−0.117 to −0.013). With k=6, it improved 50% top-set overlap by 0.075 (0.023 to 0.126). At k=10, the 75% stratified model had lower RMSE than random by 0.009 (−0.019 to −0.0004); at k=12, the 50% difference was −0.020 (−0.035 to −0.007). Positive RMSE differences favor random; negative differences favor stratification. These are several exploratory comparisons across cluster counts, budgets, and metrics, so intervals are not adjusted for multiple comparisons and should not be read as confirmatory evidence.

Across cluster counts, the preferred strategy is not stable for all decision outcomes. The defensible result is that training-set design can affect retrospective selection outcomes in this panel, but the estimated benefit depends on how diversity is operationalized. A supervisor-led, pre-specified sampling criterion and independent validation are necessary before recommending one rule.

## Files and reproduction

`R/05_nested_multienv_models.R` accepts `options(rice.cluster.k = ...)` (default 8). `R/10_cluster_sensitivity_summary.R` combines the five cluster-count summaries and calculates paired environment bootstrap intervals. Machine-readable results are `reports/cluster_count_sensitivity.csv` and `reports/cluster_count_strategy_contrasts.csv`.
