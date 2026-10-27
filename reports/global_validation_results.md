# Globally held-out-line multi-environment validation

**Run status:** Completed in RStudio on 2026-10-27. Ten outer repeats were used; in each repeat, 20% of observed lines were withheld from every environment. For each training budget and strategy, phenotype records from only the selected training lines were used. Results comprise 900 environment-by-fold summaries (15 environments × 10 repeats × 3 fractions × 2 strategies).

## Models

This fixed-kernel screen used an environment mean and a genomic-by-environment kernel (`G+E+GE`) with penalty 0.1. It did not tune model type or penalty. Training-set composition was selected from genotype data alone: random sampling or proportional stratified sampling over eight genomic clusters. Test lines were absent from all training environments. The supplied outcomes are adjusted line-environment values, primarily BLUPs. A nested, training-only model-selection analysis is reported separately in `reports/nested_validation_results.md`.

## Equal-environment summary

`RMSE difference` is selected-model RMSE minus the environment-mean RMSE; a positive number means the model was worse on average. Yield regret is in the source yield units.

| Sampling strategy | Fraction of candidate lines | Mean correlation | Model RMSE | Environment-mean RMSE | RMSE difference | Top-20% overlap | Yield regret |
|---|---:|---:|---:|---:|---:|---:|---:|
| Random | 25% | 0.161 | 0.583 | 0.589 | −0.006 | 0.287 | 0.586 |
| Proportional stratified | 25% | 0.193 | 0.594 | 0.601 | −0.008 | 0.312 | 0.598 |
| Random | 50% | 0.170 | 0.582 | 0.580 | +0.002 | 0.310 | 0.575 |
| Proportional stratified | 50% | 0.177 | 0.590 | 0.590 | −0.000 | 0.314 | 0.583 |
| Random | 75% | 0.199 | 0.583 | 0.579 | +0.004 | 0.323 | 0.535 |
| Proportional stratified | 75% | 0.205 | 0.580 | 0.574 | +0.006 | 0.307 | 0.555 |

## Interpretation

The fixed kernel slightly improved mean RMSE over the environment mean at 25%, was approximately tied for proportional stratification at 50%, and was worse at 75%. Proportional stratification had modestly higher correlation in this screen, but random sampling had better top-set overlap and lower regret at 50% and 75%. These results are a screen only; strategy effects should be interpreted from the nested validation and environment-paired uncertainty analysis.

The small panel, a 20% held-out set, only ten outer repeats, and short inner tuning make this a strong diagnostic rather than a definitive optimization benchmark. Treat each environment as the main unit for uncertainty. Examine environment-specific results and the selected model/penalty frequencies before thesis conclusions. There is no prospective evaluation and no evidence here about current commercial breeding germplasm.

## Machine-readable outputs

- `results/global_multienv_cv_by_environment.csv`: environment-level metrics for each outer repeat, budget, and strategy.
- `results/global_multienv_cv_environment_summary.csv`: repeat-averaged environment results.
- `results/global_multienv_cv_summary.csv`: equal-environment overall summary.



