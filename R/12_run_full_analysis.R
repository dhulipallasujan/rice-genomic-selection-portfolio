# Rebuild the complete analysis in order from the downloaded source supplements.
# Run from the RStudio project root. Source data are not downloaded by this script.
required <- c("data/raw/12284_2023_623_MOESM1_ESM.xlsx",
              "data/raw/12284_2023_623_MOESM2_ESM.pdf",
              "data/raw/12284_2023_623_MOESM3_ESM.txt",
              "data/raw/12284_2023_623_MOESM4_ESM.csv",
              "data/raw/12284_2023_623_MOESM5_ESM.csv",
              "data/raw/12284_2023_623_MOESM6_ESM.zip")
if (!all(file.exists(required))) stop("Source supplements are incomplete. Run R/01_download_source_data.R first.")

source("R/00_audit_inputs.R")
source("R/02_prepare_data.R")
source("R/03_baseline_cv.R")
source("R/04_global_multienv_cv.R")
options(rice.duplicate.policy = "mean", rice.cluster.k = 8L)
source("R/05_nested_multienv_models.R")
source("R/06_descriptive_eda.R")
source("R/08_uncertainty_analysis.R")
options(rice.duplicate.policy = "exclude", rice.cluster.k = 8L)
source("R/05_nested_multienv_models.R")
options(rice.duplicate.policy = "mean", rice.cluster.k = 8L)
source("R/09_duplicate_sensitivity.R")
for (kk in c(4L, 6L, 10L, 12L)) {
  options(rice.cluster.k = kk)
  source("R/05_nested_multienv_models.R")
}
options(rice.duplicate.policy = "mean", rice.cluster.k = 8L)
source("R/10_cluster_sensitivity_summary.R")
source("R/11_portfolio_figures.R")
message("Full analysis and sensitivity pipeline completed.")
