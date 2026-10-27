# Inventory source files and audit the published data tables.
# Run from the project root: Rscript R/00_audit_inputs.R

required_packages <- c("data.table", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
}

files <- list.files("data/raw", recursive = TRUE, full.names = TRUE,
                    all.files = FALSE, include.dirs = FALSE)
files <- files[!grepl("README\\.md$", files, ignore.case = TRUE)]
if (!length(files)) {
  stop("No source files found in data/raw. See data/raw/README.md for the source list.")
}

# Read bytes directly before hashing. digest(file = TRUE) relies on
# file.access(), which can report unreadable for readable files on Windows.
sha256_file <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con))
  digest::digest(readBin(con, what = "raw", n = as.integer(file.info(path)[["size"]])),
                 algo = "sha256", serialize = FALSE)
}

manifest <- data.frame(
  file = basename(files),
  source_url = NA_character_,
  access_date = NA_character_,
  bytes = file.info(files)$size,
  modified_utc = format(file.info(files)$mtime, tz = "UTC", usetz = TRUE),
  sha256 = vapply(files, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
if (file.exists("data/source_manifest.csv")) {
  prior <- tryCatch(read.csv("data/source_manifest.csv", stringsAsFactors = FALSE),
                    error = function(e) NULL)
  if (!is.null(prior) && all(c("file", "source_url", "access_date") %in% names(prior))) {
    idx <- match(manifest$file, prior$file)
    manifest$source_url <- prior$source_url[idx]
    manifest$access_date <- prior$access_date[idx]
  }
}
write.csv(manifest, "data/source_manifest.csv", row.names = FALSE)
print(manifest, row.names = FALSE)

text_tables <- list()
summary_rows <- list()

for (path in files) {
  extension <- tolower(tools::file_ext(path))
  cat("\n---", path, "---\n")
  if (extension %in% c("csv", "txt", "tsv")) {
    # The publisher's .txt genotype file is HapMap, which is tab-delimited.
    sep <- if (extension %in% c("txt", "tsv")) "\t" else ","
    tab <- tryCatch(
      data.table::fread(path, sep = sep, data.table = FALSE,
                        check.names = FALSE, fill = TRUE, na.strings = c("", "NA", ".")),
      error = function(e) {
        message("Could not parse as delimited text: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(tab)) {
      cat("Full dimensions:", nrow(tab), "rows x", ncol(tab), "columns\n")
      cat("Column names:\n")
      print(names(tab))
      cat("Missing values by column:\n")
      print(sort(colSums(is.na(tab)), decreasing = TRUE))
      text_tables[[basename(path)]] <- tab
      summary_rows[[basename(path)]] <- data.frame(
        file = basename(path), rows = nrow(tab), columns = ncol(tab),
        missing_cells = sum(is.na(tab)), stringsAsFactors = FALSE
      )
    }
  } else if (extension == "zip") {
    cat("Archive contents (not extracted or executed):\n")
    print(utils::unzip(path, list = TRUE)[, c("Name", "Length")], row.names = FALSE)
  } else {
    cat("Non-tabular source; retained for manual reference.\n")
  }
}

if (all(c("12284_2023_623_MOESM3_ESM.txt", "12284_2023_623_MOESM4_ESM.csv") %in% names(text_tables))) {
  geno <- text_tables[["12284_2023_623_MOESM3_ESM.txt"]]
  pheno <- text_tables[["12284_2023_623_MOESM4_ESM.csv"]]
  hapmap_metadata <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                       "center", "protLSID", "assayLSID", "panelLSID", "QCcode")
  marker_columns <- intersect(hapmap_metadata, names(geno))
  line_ids <- setdiff(names(geno), hapmap_metadata)
  phenotype_ids <- unique(as.character(pheno$GID))
  cat("\n--- Cross-file checks ---\n")
  cat("HapMap markers:", nrow(geno), "\n")
  cat("HapMap line columns:", length(line_ids), "\n")
  cat("Phenotype unique lines:", length(phenotype_ids), "\n")
  cat("Phenotype rows:", nrow(pheno), "\n")
  cat("Phenotype environments (Location × Year × Season):",
      length(unique(paste(pheno$Location, pheno$Year, pheno$Season, sep = "|"))), "\n")
  cat("Phenotype traits:\n")
  print(sort(unique(pheno$Trait)))
  cat("Lines with phenotype but no genotype:\n")
  print(setdiff(phenotype_ids, line_ids))
  cat("Genotyped lines with no phenotype:\n")
  print(setdiff(line_ids, phenotype_ids))
  cat("Duplicate line-environment-trait records:",
      sum(duplicated(pheno[c("GID", "Location", "Year", "Season", "Trait")])), "\n")
}

if (length(summary_rows)) {
  dir.create("reports", showWarnings = FALSE)
  write.csv(do.call(rbind, summary_rows), "reports/input_audit_summary.csv", row.names = FALSE)
}

cat("\nAudit complete. No cleaning, model fitting, or outcome calculation is performed here.\n")

