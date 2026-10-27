# Source data audit and preparation record

Audit and preparation were executed in RStudio on 2026-10-27 using the six supplements for Nguyen et al. (2023). The original OneDrive project was left untouched; this is the reproducible working copy.

## Source inventory and checks

| Source | Shape / contents | Result |
|---|---:|---|
| HapMap genotypes (MOESM3) | 882 marker rows × 107 genotype lines plus metadata | All 882 markers are polymorphic; no missing calls were observed. |
| Phenotypes (MOESM4) | 3,186 long records | Three traits across 15 environments and 107 lines; mostly adjusted BLUPs, not raw plot measurements. |
| Environmental covariates (MOESM5) | 15 environments × 24 covariates | No missing values in the supplied table. |
| Supplementary workbook (MOESM1) | Trial and trait documentation | Retained for interpretation and later supervisor review. |
| ZIP (MOESM6) | Five author R scripts | Inspected as text, not blindly executed. The authors average duplicate cells in their wide pivot. |

All 107 phenotype IDs match genotype IDs. There are 15 duplicated line-environment-trait keys for line 4016512 across five environments and all three traits. The primary dataset averages these records to follow the published preparation rule. The sensitivity dataset excludes each duplicated key entirely. The duplicate keys and paired results are retained in reports and machine-readable output.

## Reproducible preparation

`R/01_download_source_data.R` downloads the supplements and records URL, date, file size, and SHA-256. `R/00_audit_inputs.R` checks file integrity, dimensions, missingness, IDs, environments, and duplicates without executing archived code. `R/02_prepare_data.R` converts each biallelic call to dosage 0, 0.5, or 1 according to the HapMap allele order, preserves missing calls as NA, and writes the marker map, analysis-ready phenotype tables, and preparation QC.

The source files remain in `data/raw/`. Consult the publisher's reuse terms before redistributing them.

## Analysis completed

- A 30-repeat within-environment whole-line holdout baseline.
- A fixed genomic-by-environment kernel screen with ten global line-holdout repeats.
- A nested globally held-out-line comparison: ten outer repeats, three inner grouped folds, model and penalty selection on training lines only, three training budgets, and random versus proportional genotype-stratified sampling.
- Environment-level bootstrap uncertainty for paired strategy contrasts.
- Duplicate-averaging versus duplicated-key exclusion sensitivity.
- Cluster-count sensitivity at k = 4, 6, 8, 10, and 12 with paired environment bootstrap intervals.
- Descriptive phenotype coverage, yield distribution, and genotype PCA figures.

See `reports/portfolio_case_study.md` for the synthesis and `README.md` for reproduction steps. The central result is qualified: the eight-cluster primary setting does not beat the environment-mean RMSE benchmark at any tested budget, but some alternative cluster cuts produce favorable comparisons. The preferred training strategy therefore depends on how diversity is defined.

## Limitations and further work

The panel has only 107 lines and 15 environments, the trial layout is uneven, and most phenotype targets are adjusted BLUPs. The work supports a transparent retrospective case study, not plot-level accuracy, realized prospective gain, or operational claims for a commercial program. Further thesis work should expand the literature search, evaluate alternative relationship-based training-population criteria, test an environmental-covariate model under leave-one-environment-out validation, and seek independent or prospective validation.
