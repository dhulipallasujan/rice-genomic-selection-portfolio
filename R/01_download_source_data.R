# Download the original supplementary files for Nguyen et al. (2023).
# Run from the project root: Rscript R/01_download_source_data.R

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Install the 'digest' package first: install.packages('digest')")
}

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
doi_path <- "art%3A10.1186%2Fs12284-023-00623-6/MediaObjects/"
publisher_root <- paste0("https://static-content.springer.com/esm/", doi_path)
files <- c(
  "12284_2023_623_MOESM1_ESM.xlsx",
  "12284_2023_623_MOESM2_ESM.pdf",
  "12284_2023_623_MOESM3_ESM.txt",
  "12284_2023_623_MOESM4_ESM.csv",
  "12284_2023_623_MOESM5_ESM.csv",
  "12284_2023_623_MOESM6_ESM.zip"
)

manifest_path <- "data/source_manifest.csv"
if (file.exists(manifest_path)) {
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
} else {
  manifest <- data.frame(
    file = character(), source_url = character(), access_date = character(),
    bytes = numeric(), sha256 = character(), stringsAsFactors = FALSE
  )
}

for (name in files) {
  destination <- file.path("data/raw", name)
  if (!file.exists(destination)) {
    url <- paste0(publisher_root, utils::URLencode(name, reserved = TRUE))
    message("Downloading ", name)
    status <- tryCatch(
      utils::download.file(url, destination, mode = "wb", quiet = FALSE),
      error = function(e) stop("Download failed for ", name, ": ", conditionMessage(e))
    )
    if (!identical(status, 0L) || !file.exists(destination) || file.info(destination)$size == 0) {
      stop("Download did not produce a valid file: ", name)
    }
  }

  url <- paste0(publisher_root, utils::URLencode(name, reserved = TRUE))
  record <- data.frame(
    file = name, source_url = url, access_date = as.character(Sys.Date()),
    bytes = unname(file.info(destination)$size),
    sha256 = digest::digest(destination, algo = "sha256", file = TRUE),
    stringsAsFactors = FALSE
  )
  manifest <- manifest[manifest$file != name, , drop = FALSE]
  manifest <- rbind(manifest, record)
}

manifest <- manifest[match(files, manifest$file), , drop = FALSE]
write.csv(manifest, manifest_path, row.names = FALSE, na = "")
print(manifest, row.names = FALSE)
message("Download and provenance manifest complete.")
