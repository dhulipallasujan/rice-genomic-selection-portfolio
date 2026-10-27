# Source data

The candidate source is Nguyen et al. (2023), *Multi-environment Genomic Selection in Rice Elite Breeding Lines*, DOI [10.1186/s12284-023-00623-6](https://doi.org/10.1186/s12284-023-00623-6).

The article is open access under CC BY 4.0. The article's supplementary files include:

- Additional file 1 (XLSX): panel, trial design, environmental covariates, and published prediction summaries.
- Additional file 3 (TXT): 882-SNP genotypes for the elite breeding lines.
- Additional file 4 (CSV): adjusted phenotypes for three traits across 15 environments.
- Additional file 5 (CSV): environmental covariates.
- Additional file 6 (ZIP): scripts used in the source study.

Run `Rscript R/01_download_source_data.R` from the project root to download the original supplementary files using the publisher's stable supplementary-file paths. The script does not overwrite existing files. It writes `data/source_manifest.csv` with the source URL, download date, byte size, and SHA-256 checksum. The source data are not included in this repository scaffold.

Retain attribution to the original authors and cite the article. If the supplementary data include any separately credited third-party material, follow its stated terms. The project will avoid committing raw files by default and will instead preserve download instructions and checksums.

