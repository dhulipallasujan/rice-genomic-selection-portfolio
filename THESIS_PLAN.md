# Thesis plan

## Proposed title

**Choosing the training population, not just the prediction model: phenotyping-budget strategies for genomic selection in elite rice across target environments**

## Why this question

Genomic prediction is useful to breeders only when it helps them decide which lines to advance. Prediction accuracy by itself does not tell a program how many lines it must phenotype, whether its training set represents the candidates, or whether the selected lines remain useful across target environments.

The public IRRI elite-core-panel resource offers a realistic case: a starting panel of 111 elite lines, 15 trials, 12 locations, three years, and an unbalanced genotype-by-environment table. After source-study marker quality control, 107 lines and 882 SNPs remained. Supplementary data also include three adjusted traits and environmental covariates. The original article already evaluated genomic prediction and environmental effects. The proposed thesis focuses on a complementary decision question: the cost and benefit of alternative training-line sampling strategies, measured by selection outcomes under strict validation.

## Primary research question

At different phenotyping-budget levels, which training-line sampling strategy best preserves genomic selection gain for grain yield when predicting untested lines across target environments?

## Objectives

1. Reconstruct and document a clean, auditable genotype–phenotype–environment dataset from the published supplementary files.
2. Estimate a transparent baseline for grain-yield prediction using environment means and genomic main effects.
3. Compare random, genomic-diversity, and prediction-oriented training-line selection at multiple training-set sizes.
4. Test performance under validation designs that represent distinct breeding uses: new lines in represented environments and unobserved line-by-environment combinations; hold out whole environments only where data support a defensible evaluation.
5. Evaluate breeder-facing outcomes: mean yield of the selected top fraction, selection regret, overlap of selected lines across environments, and uncertainty around those outcomes.
6. Examine whether flowering time and plant height provide useful auxiliary information for yield prediction without allowing target-environment phenotypes to leak into training.

## Working hypotheses

- **H1:** At a fixed phenotyping budget, genomic-diversity or prediction-oriented sampling can improve yield selection outcomes over purely random sampling, but the benefit will vary by target environment.
- **H2:** Ranking methods by correlation alone can produce a different conclusion from ranking them by top-fraction yield gain and selection regret.
- **H3:** Training-set size will show diminishing returns, and the point at which gains level off will depend on the validation scenario and environment relatedness.
- **H4:** Auxiliary traits may help only when their genetic and phenotypic relationships to yield are sufficiently strong in the training data; any benefit will be reported by environment rather than assumed to generalize.

These are hypotheses to test, not expected or promised results.

## Data and units

The source article reports an elite panel of 111 irrigated-rice lines evaluated in 15 location-year-season environments across 12 sites during 2018–2020. Its marker QC retained 107 lines and 882 SNPs. Supplementary materials include adjusted phenotypes (days to flowering, plant height, grain yield), environment identities and trial details, and environmental covariates. The core observation is the line × environment record. Before analysis, the exact numbers of lines, markers, environments, replications, missing cells, and matched identifiers must be recomputed from the files.

The adjusted phenotype table is not raw plot data. The thesis must describe that limitation and must not claim to re-estimate plot-level heritability or spatial field effects from BLUPs alone.

## Analysis design

### Data audit and quality control

- Verify source file checksums or record file sizes and dates; preserve a source manifest.
- Check duplicate identifiers, duplicated line–environment cells, missing trait values, marker coding, allele frequencies, marker missingness, and line overlap.
- Report population structure and genomic relatedness before defining grouped splits.
- Confirm whether environmental covariates are available for all target trials and inspect missingness before proposing covariate models.
- Save the cleaned analysis tables and a data dictionary; keep downloaded source files out of version control unless redistribution terms clearly allow it.

### Prediction benchmarks

Use a small, interpretable model set:

1. Environment mean baseline.
2. Genomic main-effect ridge/GBLUP model.
3. One multi-environment genomic model, selected after checking the source's model definitions and the available R packages. A reaction-norm or genomic × environment model is a candidate, not a commitment.
4. A multi-trait model only if training-only trait relationships and the incomplete data layout support a fair comparison.

Avoid deep learning: the number of lines and environments is too small to justify a large neural model, and the thesis should prioritize design quality and breeding interpretation.

### Training-population strategies

- Random sampling, repeated across seeds.
- Genomic-diversity sampling, such as a maximin or cluster-balanced selection, with the rule fixed before test outcomes are viewed.
- Prediction-oriented sampling such as CDmean only if its implementation can be independently checked and the required target genotypes are available without using target phenotypes.
- Training fractions of 20%, 40%, 60%, and 80%, adjusted if a fraction leaves too few lines or makes a validation split unstable.

### Validation scenarios

Each scenario answers a separate question and will be labeled clearly:

- **New-line validation:** held-out lines are absent from training across the observed environments; related lines are grouped where sample size permits.
- **Sparse testing:** predict missing line × environment cells when information from other environments is available for some lines.
- **New-environment validation:** hold out an entire environment and predict it only if a model has a principled way to handle an unseen environment (for example, measured environmental covariates). A categorical target-environment effect cannot be estimated from data that were fully withheld.
- **Temporal validation:** use earlier years to predict later trials only if the trial layout and number of environments provide enough support. If not, report this scenario as infeasible rather than forcing it.

All hyperparameter selection and any marker filtering must be nested within the training split. The same held-out rows must be used for all compared strategies.

### Primary outcomes

- Mean observed grain yield among the predicted top 10% and top 20% of candidates, compared with the mean yield of the same-sized random selections and the observed top set.
- Selection regret: gap between the observed mean of the true top set and the mean of the model-selected set.
- Rank correlation / Pearson correlation and RMSE as secondary prediction metrics.
- Selection stability: overlap of selected lines across repeated splits and target environments.
- Uncertainty intervals from repeated sampling and environment-level summaries; do not treat thousands of line–environment records as independent replicates.

### Sensitivity and limitations

Repeat key results across reasonable marker QC thresholds, selection fractions, and grouped split definitions. With only 111 lines and 15 trials, estimates may be noisy. This is a public historical panel, not a prospective test of current ANGRAU material. The main result will be a transparent case study and workflow, not a claim that one strategy is universally best for rice breeding.

## Twelve-month deliverables

- Month 2: literature matrix, data dictionary, reproducibility and reuse check, final supervisor-approved protocol.
- Month 4: audited analysis-ready tables and descriptive report.
- Month 6: baseline genomic prediction pipeline and validation plan.
- Month 8: completed sampling-strategy comparisons.
- Month 10: selection-gain, stability, and uncertainty results.
- Month 11: full thesis draft and manuscript draft.
- Month 12: clean rerun, final thesis, defense slides, and a public code release if permitted by the university and data terms.

## Stop/go checks

Before committing to the main analysis, confirm that the supplementary data are downloadable and useable, that the genotype and phenotype line IDs match, and that the trial layout supports the proposed splits. Ask the supervisor to review novelty against the source paper and current literature. If line sampling cannot be simulated without leakage or the data are too sparse for one validation scenario, narrow that objective while retaining the central selection-budget question.


## Completed analysis update (2026-10-27)

The working portfolio now includes audited and prepared data, a within-environment baseline, a global line-held-out fixed-kernel screen, a nested globally grouped validation with training-only tuning, paired environment bootstrap uncertainty, duplicate-key sensitivity, five cluster-count settings, and reproducible figures. The current primary analysis uses 10 outer repeats, 20% globally held-out lines, three inner grouped folds, three candidate budgets, two sampling rules, two model structures, and three ridge penalties. See `reports/portfolio_case_study.md` for the results narrative and `README.md` for reproduction steps.

The results refine the original hypotheses. At the primary eight-cluster setting, the tuned model did not beat the environment-mean RMSE benchmark at any budget. The relative value of proportional stratification changes across cluster counts and decision outcomes, so no single strategy is supported as a universal winner. The duplicate-handling sensitivity is small. The main remaining thesis work is expanding the review, evaluating alternative relationship-based sampling and a properly validated environmental-covariate leave-one-environment-out model, then obtaining independent or prospective validation before operational recommendations.
