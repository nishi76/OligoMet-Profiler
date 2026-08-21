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

cat("\nDone. Check for python3 on PATH separately -- it is used for mzML\n")
cat("binary decoding (falls back to a pure-R decoder if unavailable).\n")
