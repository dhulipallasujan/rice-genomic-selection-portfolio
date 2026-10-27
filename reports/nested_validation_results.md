# Nested, globally held-out-line validation

## Design and run

Executed in RStudio on 2026-10-27. Ten outer repeats held out 20% of observed genotype lines globally, removing them from every training environment. For each repeat, random or proportional stratified training samples were drawn at 25%, 50%, or 75% of the remaining candidate lines. Model type (`G+E` or `G+E+GE`) and ridge penalty (0.01, 0.1, 1) were selected by three-fold inner validation grouped by line using training data only. The outer test phenotypes were never used for sample, model, or penalty selection. The environment-specific mean was the primary baseline. Environment-level metrics were averaged equally in the overall summary.

The cluster assignments use a fixed eight-group hierarchical clustering of marker relationship distances. The PCA diagnostic shows that two PCs explain 17.4% of marker variance and the eight groups are uneven; they are computational strata, not biological labels. Alternative cuts at k = 4, 6, 10, and 12 are now evaluated in `reports/cluster_count_sensitivity.md`; the strategy contrast changes with the cut.

## Results

`RMSE difference` is nested-model RMSE minus the environment-mean baseline, so positive values favor the mean baseline. Higher correlation and top-set overlap are favorable; lower regret is favorable.

| Sampling | Candidate fraction | Correlation | Nested model RMSE | Environment-mean RMSE | RMSE difference | Top-20% overlap | Yield regret |
|---|---:|---:|---:|---:|---:|---:|---:|
| Random | 25% | 0.076 | 0.609 | 0.595 | +0.014 | 0.273 | 0.622 |
| Proportional stratified | 25% | 0.024 | 0.616 | 0.602 | +0.013 | 0.250 | 0.633 |
| Random | 50% | 0.059 | 0.594 | 0.592 | +0.003 | 0.268 | 0.626 |
| Proportional stratified | 50% | 0.048 | 0.591 | 0.584 | +0.007 | 0.264 | 0.646 |
| Random | 75% | 0.115 | 0.578 | 0.577 | +0.001 | 0.331 | 0.572 |
| Proportional stratified | 75% | 0.087 | 0.582 | 0.579 | +0.002 | 0.281 | 0.647 |

The inner procedure selected the `G+E+GE` model and penalty 1 most often at 50% and 75%; at 25%, selection was mixed. Full frequencies are in `results/nested_multienv_model_selection_counts.csv`.

## Paired environment-level strategy contrast

I averaged the ten outer repeats within each environment, then bootstrapped the 15 environments 10,000 times. Values below are proportional-stratified minus random, with percentile 95% intervals. Negative top-set overlap differences and positive regret differences favor random sampling.

| Candidate fraction | Top-20% overlap difference (95% interval) | Regret difference (95% interval) |
|---|---:|---:|
| 25% | −0.023 (−0.069, 0.024) | +0.011 (−0.014, 0.036) |
| 50% | −0.004 (−0.045, 0.038) | +0.020 (−0.041, 0.088) |
| 75% | −0.049 (−0.092, −0.009) | +0.076 (+0.006, +0.168) |

At the 75% budget, random training selection had better top-set overlap and lower regret in this bootstrap. At 25% and 50%, intervals include zero. With only 15 environments, these intervals are uncertain and should be treated as exploratory rather than definitive.

## Interpretation

After nested tuning, the genomic model did not improve equal-environment RMSE over the environment-mean baseline at any budget. At the largest budget, the gap was small. Random sampling had better average top-20% overlap and lower regret than proportional stratification at 50% and 75%; the two methods were similarly weak at 25%. These outcomes do not support a claim that diversity sampling improves this panel's selection decisions. The cluster-count sensitivity shows that some stratified settings outperform the baseline at selected budgets, while others do not; the eight-cluster result should not be generalized. See `reports/cluster_count_sensitivity.md`.

This is a historical retrospective exercise with 107 lines, mostly adjusted BLUP targets, 15 uneven environments, and ten outer repeats. Top-20% outcomes are noisy where a test environment contains few lines. It is not a prospective estimate of genetic gain or evidence about any company's current breeding program.

## Reproducible outputs

- `R/05_nested_multienv_models.R`: complete nested pipeline.
- `results/nested_multienv_cv_by_environment.csv`: environment-level results by outer repeat.
- `results/nested_multienv_cv_environment_summary.csv`: repeat averages for each environment and strategy.
- `results/nested_multienv_cv_summary.csv`: equal-environment summary.
- `results/nested_multienv_model_selection_counts.csv`: inner-selected model and penalty frequencies.
- `R/08_uncertainty_analysis.R`: paired environment bootstrap.
- `reports/paired_environment_strategy_contrasts.csv`: strategy contrasts and bootstrap intervals.


