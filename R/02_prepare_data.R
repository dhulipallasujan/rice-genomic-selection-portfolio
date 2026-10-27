# Prepare the published rice panel for reproducible genomic prediction.
# Genotypes are encoded 0/0.5/1 for first-allele homozygote/heterozygote/
# second-allele homozygote, matching the numerical convention used in the
# source paper. Ambiguous and missing calls remain NA.
suppressPackageStartupMessages({library(data.table)})
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("reports", showWarnings = FALSE, recursive = TRUE)

hapmap <- fread("data/raw/12284_2023_623_MOESM3_ESM.txt", sep = "\t", na.strings = c("NA", "N", "-", "NN"), check.names = FALSE)
pheno <- fread("data/raw/12284_2023_623_MOESM4_ESM.csv", na.strings = c("", "NA"))
stopifnot(all(c("rs#", "alleles", "chrom", "pos") %in% names(hapmap)))
line_ids <- setdiff(names(hapmap), c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID", "assayLSID", "panelLSID", "QCcode"))
stopifnot(length(line_ids) == 107L, nrow(hapmap) == 882L)

decode_call <- function(call, allele_pair) {
  if (is.na(call) || call %in% c("N", "-", "NN", "NA", "./.")) return(NA_real_)
  alleles <- strsplit(gsub("[^ACGT]", "", allele_pair), "", fixed = TRUE)[[1]]
  if (length(alleles) != 2L || alleles[1] == alleles[2]) return(NA_real_)
  call <- gsub("[/|]", "", call)
  if (nchar(call) == 1L) {
    if (call == alleles[1]) return(0)
    if (call == alleles[2]) return(1)
    return(NA_real_)
  }
  if (nchar(call) != 2L || any(!strsplit(call, "", fixed = TRUE)[[1]] %in% alleles)) return(NA_real_)
  sum(strsplit(call, "", fixed = TRUE)[[1]] == alleles[2]) / 2
}
G <- matrix(NA_real_, nrow = nrow(hapmap), ncol = length(line_ids), dimnames = list(hapmap[["rs#"]], line_ids))
for (j in seq_along(line_ids)) G[, j] <- mapply(decode_call, hapmap[[line_ids[j]]], hapmap$alleles)
keep_marker <- apply(G, 1, function(x) any(!is.na(x)) && length(unique(x[!is.na(x)])) > 1L)
G <- G[keep_marker, , drop = FALSE]
if (anyDuplicated(rownames(G))) stop("Duplicate marker IDs; resolve before modeling.")

# Preserve the source-preparation rule: average duplicate GID x environment x trait values.
pheno[, EnvID := paste(Country, Location, Year, Season, sep = "_")]
pheno <- pheno[!is.na(BLUPs)]
dup <- pheno[, .N, by = .(GID, EnvID, Trait)][N > 1]
fwrite(dup, "reports/duplicate_phenotype_keys.csv")
pheno_wide <- dcast(pheno, GID + Country + Location + Year + Season + EnvID ~ Trait,
                    value.var = "BLUPs", fun.aggregate = mean, fill = NA_real_)
setnames(pheno_wide, "GID", "LineID")
pheno_wide[, LineID := as.character(LineID)]
if (!setequal(pheno_wide$LineID, line_ids)) stop("Phenotype and genotype line IDs do not match.")
dup_keys <- pheno[, .N, by = .(GID, Country, Location, Year, Season, Trait)][N > 1L,
  .(GID, Country, Location, Year, Season, Trait)]
pheno_exclude <- pheno[!dup_keys, on = .(GID, Country, Location, Year, Season, Trait)]
pheno_wide_exclude <- dcast(pheno_exclude, GID + Country + Location + Year + Season + EnvID ~ Trait,
                            value.var = "BLUPs", fill = NA_real_)
setnames(pheno_wide_exclude, "GID", "LineID")
pheno_wide_exclude[, LineID := as.character(LineID)]
Gout <- G
colnames(Gout) <- paste0("L", line_ids) # safe, stable names; original IDs remain explicit
fwrite(as.data.table(Gout, keep.rownames = "MarkerID"), "data/processed/genotype_dosage.csv")
fwrite(pheno_wide, "data/processed/phenotypes_wide.csv")
fwrite(pheno_wide_exclude, "data/processed/phenotypes_wide_exclude_dupkeys.csv")
fwrite(data.table(MarkerID = rownames(G), Chromosome = hapmap$chrom[keep_marker], Position = hapmap$pos[keep_marker], Alleles = hapmap$alleles[keep_marker]), "data/processed/marker_map.csv")
qc <- data.table(metric = c("source_markers", "retained_polymorphic_markers", "lines", "phenotype_rows", "phenotype_duplicate_keys_averaged", "genotype_missing_fraction"),
                 value = c(nrow(hapmap), nrow(G), length(line_ids), nrow(pheno), nrow(dup), mean(is.na(G))))
fwrite(qc, "reports/preparation_qc.csv")
message("Prepared ", nrow(G), " polymorphic markers across ", length(line_ids), " lines; genotype missing fraction = ", signif(mean(is.na(G)), 4))
