# Rice genomic selection portfolio case study

## Decision question

Can choosing a smaller genotype-diverse training population improve grain-yield prediction or selection decisions for unphenotyped lines across the environments represented in this published rice panel?

## Data and quality controls

The source release contains 107 genotyped elite lines and 882 polymorphic SNPs, with adjusted phenotypes for days to flowering, plant height, and grain yield across 15 location-year-season environments. There are 1,061 unique line-environment rows, with uneven trial coverage. Genotype and phenotype IDs matched completely. The marker data were encoded as allele dosage from 0 to 1 using the published allele order; missing marker calls were absent in the source panel. The phenotype supplements are mostly adjusted BLUPs and should not be mistaken for raw plot measurements.

The source table contains 15 duplicated line-environment-trait keys for line 4016512 across five environments and three traits. The primary analysis follows the authors' preparation rule and averages duplicate values. A sensitivity analysis excludes all values for duplicate keys. The results are similar under both policies; the duplicate policy does not change the central conclusion. Full key-level and paired outputs are retained.

## Validation and analysis

The primary design repeats a global 20% line holdout ten times. Each held-out line is absent from every training environment. Among remaining candidates, random and proportional genotype-cluster-stratified training sets are compared at 25%, 50%, and 75% budgets. The model and ridge penalty are selected in three inner folds grouped by line, using training lines only. Performance is summarized equally across environments rather than weighting the largest trials most heavily.

Outcomes are Pearson correlation, RMSE relative to an environment-mean benchmark, top-20% overlap with the observed best set, and selection regret. The latter two are historical retrospective indicators based on the published adjusted targets; they are not realized genetic gain. Strategy contrasts use paired environment averages with percentile bootstrap intervals over 15 environments.

## Results

| Training selection | Budget | Correlation | Model RMSE | Environment-mean RMSE | RMSE difference | Top-20% overlap | Regret |
|---|---:|---:|---:|---:|---:|---:|---:|
| Random | 25% | 0.076 | 0.609 | 0.595 | +0.014 | 0.273 | 0.622 |
| Proportional stratified | 25% | 0.024 | 0.616 | 0.602 | +0.013 | 0.250 | 0.633 |
| Random | 50% | 0.059 | 0.594 | 0.592 | +0.003 | 0.268 | 0.626 |
| Proportional stratified | 50% | 0.048 | 0.591 | 0.584 | +0.007 | 0.264 | 0.646 |
| Random | 75% | 0.115 | 0.578 | 0.577 | +0.001 | 0.331 | 0.572 |
| Proportional stratified | 75% | 0.087 | 0.582 | 0.579 | +0.002 | 0.281 | 0.647 |

Positive RMSE difference means the genomic model has larger error than the environment-mean benchmark. At the primary eight-cluster setting, the tuned kernel did not improve equal-environment RMSE at any budget. The environment-bootstrap interval for stratified-minus-random top-set overlap at 75% was −0.093 to −0.009; for regret it was +0.003 to +0.169. At 25% and 50%, intervals for these strategy contrasts included zero. However, the cluster-count sensitivity changes the picture: k=4, 6, 10, and 12 produced favorable stratified comparisons for different budgets and outcomes, including lower regret or RMSE in selected settings. These exploratory intervals are not adjusted for multiple comparisons. The data therefore do not support a universal winner; the estimated effect depends on how genotype diversity is defined. See `reports/cluster_count_sensitivity.md` and its paired results table.

Excluding duplicated keys changed average nested model RMSE by less than 0.003 at every strategy and budget. It changed average top-20% overlap and regret by at most about 0.007 and 0.005, respectively. The main comparison is therefore not driven by the duplicate averaging decision.

The simpler fixed-kernel screen and within-environment baseline are included as diagnostics, not competing primary claims. Results and definitions are in `reports/global_validation_results.md` and `reports/baseline_cv_results.md`.

![Equal-environment RMSE difference across cluster counts; zero is the environment-mean benchmark.](../results/figures/cluster_sensitivity_rmse.png)

## Interpretation for a breeding analytics audience

In this dataset, a plausible genomic modeling pipeline was not automatically better than a simple environment-specific mean. Proportional genotype-cluster sampling did not provide a stable advantage across the tested cluster cuts, budgets, and decision outcomes. This is useful evidence: it shows why a breeding decision should be evaluated using a deployment-matched holdout, a simple benchmark, and selection outcomes alongside prediction correlation. It also shows the value of recording uncertainty and data-cleaning sensitivity rather than reporting a single best score.

For an applied breeding team, the transferable skills demonstrated here are reproducible data ingestion and QC, marker/phenotype alignment, grouped cross-validation, training-only tuning, genotype-based training-set design, multi-environment summaries, breeder-facing selection metrics, and communication of limitations. The current public case study does not establish performance on a company's proprietary germplasm or current breeding pipeline.

## Limitations

- The panel is small (107 lines) and represents one historical elite panel; test estimates may not transfer to other breeding populations.
- Only 15 environments support the main uncertainty analysis, and several trials have substantially fewer observed lines than others.
- Most phenotype targets are adjusted BLUPs. They are not raw plot data, and no plot-level heritability, spatial effects, or prospective selection gain can be estimated here.
- The stratification rule is proportional sampling from a hierarchical clustering. Sensitivity at k=4, 6, 8, 10, and 12 shows that strategy contrasts vary with the cut; alternative relationship-based criteria remain to be assessed.
- Ten outer repeats and a short inner tuning grid make this a rigorous portfolio demonstration, not a final operational recommendation.
- The literature positioning is focused rather than systematic; the source paper already evaluates this exact panel under other validation scenarios.

## Reproducibility and files

Start with `README.md`. The source supplements and checksums are in `data/raw/` and `data/source_manifest.csv`. The main analysis is `R/05_nested_multienv_models.R`; outputs are `results/nested_multienv_cv_*`. The duplicate sensitivity is `R/09_duplicate_sensitivity.R`; its paired environment results are `results/duplicate_sensitivity_by_environment.csv` and `reports/duplicate_sensitivity.csv`. Methods, audit, literature positioning, figures, and supporting model screens are retained in their project folders.

## Next steps to make this thesis-ready

1. Agree the deployment target with a supervisor and expand the literature review systematically, including how this panel study differs from prior rice training-population optimization.
2. Test alternative relationship-based training-set criteria and a pre-specified marker-QC sensitivity, maintaining identical outer held-out lines for paired comparisons.
3. Test a leave-one-environment-out model that uses the supplied environmental covariates, with all preprocessing and tuning inside training environments and no target-environment categorical effect.
4. Add a genuinely independent current breeding population or prospective season before making operational recommendations.
5. Rerun from a clean project copy, review reuse conditions, and convert this report into the thesis methods/results chapters with supervisor feedback.

## Sources

- Nguyen et al. (2023), [Multi-environment genomic selection in rice elite breeding lines](https://doi.org/10.1186/s12284-023-00623-6), source panel and data.
- Isidro et al. (2015), [Training set optimization under population structure in genomic selection](https://doi.org/10.1007/s00122-014-2418-4), rice/wheat training-population methods.
- Fernández-González et al. (2023), [A comparison of methods for training population optimization in genomic selection](https://doi.org/10.1007/s00122-023-04265-6), comparative multi-species benchmark.
- de Verdal et al. (2023), [Optimization of multi-generation multi-location genomic prediction models for recurrent genomic selection in an upland rice population](https://doi.org/10.1186/s12284-023-00661-0), a different population and deployment setting.




