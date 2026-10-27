# Portfolio pitch: genomic selection under uneven multi-environment data

## Resume-ready project description

**Rice Genomic Selection and Training-Set Design (R, genomic prediction, reproducible analytics).** Built an end-to-end pipeline around six published rice trial supplements: audited provenance and data structure, aligned 107 genotyped lines with 882 SNPs and 15 uneven environments, encoded marker dosages, and documented duplicate phenotype handling. Compared random and genotype-stratified training populations at three phenotyping budgets using globally held-out-line validation, grouped nested tuning, and a simple environment-mean benchmark. Evaluated prediction error and breeder-facing selection overlap/regret, quantified environment-level uncertainty, and tested sensitivity to duplicate handling and genotype cluster count. The project found that the preferred sampling strategy depends on how diversity is defined, demonstrating the need for deployment-matched validation and decision metrics rather than correlation alone.

## Short interview version

“I wanted to treat genomic selection as a breeding decision problem, not just a model leaderboard. I built a reproducible R workflow from the published inputs, then held out whole genotype lines across every environment so the same line could not leak into training elsewhere. Model and penalty selection happened inside grouped training folds. I compared three phenotyping budgets and evaluated RMSE against an environment baseline, plus top-set overlap and selection regret. The surprising result was that the diversity-sampling effect changed when I changed the cluster cut, so I reported that sensitivity instead of claiming one strategy was best. The practical lesson is to match validation to the deployment question and check whether the decision stays stable across reasonable design choices.”

## Skills shown in the repository

- RStudio workflow and reproducible analysis scripts
- Data ingestion, provenance checks, genotype/phenotype matching, and phenotype duplicate sensitivity
- SNP dosage encoding and genomic relationship kernels
- Grouped outer validation and nested training-only tuning to control leakage
- Multi-environment analysis with equal-environment summaries
- Breeder-oriented selection metrics: top-fraction overlap and selection regret
- Environment-level bootstrap uncertainty and design sensitivity
- Clear reporting of negative, conditional, and non-transferable findings

## How to present it accurately

Describe this as a public-data portfolio case study using one historical rice panel. Do not present it as an analysis of Syngenta data, a prospective trial, or proof that a result transfers to a commercial program. A strong next interview discussion is how you would adapt the workflow to a breeding team's target population, trial design, compute constraints, and decision costs after confirming data access and evaluation requirements.
