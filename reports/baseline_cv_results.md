# First baseline validation results

**Status:** Executed successfully in RStudio on 2026-10-27. The final run used seed 20261027 and produced 2,700 model configurations: 15 environments × 30 repeats × 3 training fractions × 2 sampling strategies. The full replicate-level results are in `results/baseline_cv_replicates.csv`; environment-specific summaries are in `results/baseline_cv_by_environment.csv`.

## Validation question

For grain yield, can marker relationships predict lines held out within an already observed environment? Each environment is modeled separately. For a given prediction, the held-out line is absent from that environment’s training set. These are not globally grouped line splits across all environments, and the predictor is not a full multi-environment model. The environment mean is the baseline. Marker relationships use a simple kernel ridge predictor with fixed regularization (0.1); random and proportional stratified training samples are compared at 25%, 50%, and 75% of the available training pool.

## Equal-environment averages

RMSE difference is predictor RMSE minus environment-mean RMSE; positive values mean the marker predictor was worse on average.

| Training strategy | Pool fraction | Mean Pearson r | Marker RMSE | Mean baseline RMSE | RMSE difference | Top-20% overlap | Selection regret |
|---|---:|---:|---:|---:|---:|---:|---:|
| Random | 25% | 0.091 | 0.601 | 0.589 | +0.013 | 0.304 | 0.562 |
| Proportional stratified | 25% | 0.125 | 0.592 | 0.588 | +0.004 | 0.302 | 0.560 |
| Random | 50% | 0.124 | 0.596 | 0.576 | +0.021 | 0.313 | 0.534 |
| Proportional stratified | 50% | 0.148 | 0.594 | 0.576 | +0.018 | 0.317 | 0.537 |
| Random | 75% | 0.158 | 0.586 | 0.569 | +0.017 | 0.333 | 0.527 |
| Proportional stratified | 75% | 0.161 | 0.587 | 0.568 | +0.019 | 0.347 | 0.507 |

Regret is in the published grain-yield units. Top-set metrics are averages across environment-specific held-out sets and should be interpreted cautiously: some environments have few held-out candidates.

## What this first run suggests

The marker model did not improve mean RMSE over the environment mean at any tested fraction. Proportional stratification was close to random sampling; at 75% it had slightly higher top-set overlap and lower regret, but only marginally higher correlation and worse RMSE. The pattern does not support a broad claim that stratified sampling is superior. The globally held-out-line analysis provides the more deployment-aligned comparison.

## Limitations and next work

- This is a baseline diagnostic, not a complete replication of the source paper’s Bayesian multi-environment models.
- Training fraction is a fraction of the non-test pool within each environment; absolute training counts vary by environment.
- A separate line can be in training for one environment model and held out for another. This is valid for the stated within-environment prediction use, but does not test completely unseen lines across all environments.
- Regularization was fixed rather than tuned in nested validation. Compare against a tuned ridge baseline with tuning confined to training lines.
- Thirty repeats improve split stability, but uncertainty remains substantial across the 15 environments; report paired environment-level differences rather than treating individual line predictions as independent.
- Compare duplicate averaging with exclusion, verify phenotype scale and source documentation, and test a globally held-out-line design. Do not claim unseen-environment performance from this analysis.
