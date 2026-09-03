# install_packages.R -- installs everything the app and pipeline need.
# Run once with: Rscript install_packages.R

required <- c("shiny", "DT", "bslib", "openxlsx", "ggplot2", "xml2", "xfun", "shinyFiles")
optional <- c("enviPat", "rmarkdown")

installed <- rownames(installed.packages())

to_install <- setdiff(required, installed)
if (length(to_install) > 0) {
  cat("Installing required packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  cat("All required packages already installed.\n")
}

missing_optional <- setdiff(optional, installed)
if (length(missing_optional) > 0) {
  cat("\nOptional packages not installed:", paste(missing_optional, collapse = ", "), "\n")
  cat("  enviPat   -- higher-accuracy isotope patterns (mass_isotope.R falls back\n")
  cat("               to a built-in convolution method if absent).\n")
  cat("  rmarkdown -- required to render HTML/PDF reports from build_report.R.\n")
  ans <- if (interactive()) readline("Install optional packages now? [y/N] ") else "n"
  if (tolower(ans) == "y") install.packages(missing_optional, repos = "https://cloud.r-project.org")
}

# Bioconductor packages (Spectra + mzR) -- robust mzML/mzXML reading
# (streaming, numpress, indexed-mzML) via R/spectra_io.R's read_ms_file().
# Optional: the app falls back to a built-in xml2-based parser
# (parse_mzml() in R/ms_matching.R) if these aren't installed, so a user who
# can't get Bioconductor/mzR working can still run everything else.
bioc_optional <- c("Spectra", "mzR")
missing_bioc <- setdiff(bioc_optional, installed)
if (length(missing_bioc) > 0) {
  cat("\nOptional Bioconductor packages not installed:", paste(missing_bioc, collapse = ", "), "\n")
  cat("  Spectra/mzR -- faster, more robust mzML/mzXML (and, via ProteoWizard\n")
  cat("                 msconvert, vendor .raw) reading. mzR is a compiled,\n")
  cat("                 platform-sensitive dependency -- installation can take\n")
  cat("                 a while and occasionally needs system libraries; the\n")
  cat("                 app works fine without it, just with a slower/less\n")
  cat("                 robust reader for large or unusual mzML files.\n")
  ans <- if (interactive()) readline("Install optional Bioconductor packages now? [y/N] ") else "n"
  if (tolower(ans) == "y") {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
  }
}

cat("\nDone.\n")
cat("\nBatch/parallel raw-file processing is a separate, optional Python\n")
cat("component (needed only for that workflow): pip install -r inst/python/requirements.txt\n")
