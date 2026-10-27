# Literature and novelty positioning

**Scope note:** This is a focused, source-led positioning review, not a systematic review. Search terms and forward/backward citation chasing should be documented and expanded with a librarian or supervisor before describing the thesis as a comprehensive review.

| Study | Population and question | Main finding relevant here | Implication for this thesis |
|---|---|---|---|
| Isidro et al. (2015), *Theoretical and Applied Genetics*, [DOI](https://doi.org/10.1007/s00122-014-2418-4) | Rice and wheat genomic selection; compared random, stratified, CDmean, PEVmean and stratified CDmean under population structure | The best training-set criterion depended on the population and trait; stratified sampling performed best in their strongly structured rice dataset. | Diversity sampling has a clear precedent. Report panel structure first; do not assume one sampling strategy should win. |
| Nguyen et al. (2023), *Rice*, [DOI](https://doi.org/10.1186/s12284-023-00623-6) | This exact 107-line, 882-marker elite core panel; multi-environment prediction of DTF, HT and yield | Evaluated seven model structures and three validation scenarios. Yield predictability varied widely by scenario; G×E models did not consistently improve on main effects. | The data and baseline question are not new. Differentiate the thesis through globally grouped line holdout, training-budget decisions, strict tuning, and decision outcomes. |
| Fernández-González et al. (2023), *Theoretical and Applied Genetics*, [open article](https://doi.org/10.1007/s00122-023-04265-6) | Seven multi-species datasets, including a 327-line rice dataset with four environments; compared many training population optimization criteria and genomic models using 40 repeats | Optimization benefits varied across datasets; the rice dataset showed relatively large gains. CDmean was strong when test-set information was available, while untargeted results depended on criterion and population structure. | Broad claims that rice training-population optimization is unstudied are not defensible. Emphasize the distinct elite panel, deployment scenario, candidate budget and breeder decision measures. |
| de Verdal et al. (2023), *Rice*, [open article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10533757/) | Multi-generation, multi-location genomic prediction in an upland rice recurrent-selection population | Evaluated prediction across generations and locations, addressing a different population and deployment problem. | Separate elite-panel within-program prediction from transfer across generations or populations. |

## Defensible contribution statement (provisional)

This project is a reproducible, decision-oriented case study of how training-line composition and budget affect grain-yield prediction and selection in one small, unbalanced elite rice panel. It compares random and genotype-diversity sampling under globally held-out-line validation, with model and penalty selection confined to training data, and reports breeder-facing ranking outcomes alongside prediction error. It is an extension and practical re-evaluation of established methods, not the invention of training-population optimization.

## Important caveats for the literature review

- The 2023 multi-species benchmark included a different rice panel and a different validation design; it does not settle which method works in this elite core panel.
- Isidro et al.'s rice result under strong population structure cannot be transferred automatically to this ECP. The source article characterized its panel as a narrow indica elite group, so this project must quantify its own structure.
- Targeted training-set selection may use test-line genotypes (available before phenotyping) but must not use test phenotypes. If implemented later, it must use a training-only prediction criterion and an inner split; do not call it untargeted.
- Retrospective BLUP-based accuracy is not prospective genetic gain, and results from 107 historical lines cannot establish general performance at Syngenta or any current breeding program.

