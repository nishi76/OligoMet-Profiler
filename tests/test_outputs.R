# test_outputs.R -- test workbook + report generation for the inotersen reference case

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
source(file.path(.pkg_root, "R", "chemistry_dict.R"))
source(file.path(.pkg_root, "R", "oligo_io.R"))
source(file.path(.pkg_root, "R", "metabolites.R"))
source(file.path(.pkg_root, "R", "mass_isotope.R"))
source(file.path(.pkg_root, "R", "fragments.R"))
source(file.path(.pkg_root, "R", "ms_matching.R"))
source(file.path(.pkg_root, "R", "build_workbook.R"))
source(file.path(.pkg_root, "R", "build_report.R"))

cat("==== Testing workbook + report generation ====\n\n")

# Parse the inotersen reference case
spec <- parse_input(INOTERSEN_TRIPLET)
validation <- validate_reference(verbose = FALSE)

# Generate metabolites
mets <- generate_metabolites(spec, opts = list(
  oligo_name = "inotersen", max_3p = 10, max_5p = 10, endo = TRUE,
  endo_sites = "all", min_frag_len = 3))
cat("Generated", length(mets), "metabolites\n")

# Build options
opts <- list(
  z_range = 3:12, n_iso = 5, max_oxid = 6, h_offset = 0,
  use_envipat = FALSE, include_internal = TRUE,
  max_3p = 10, max_5p = 10, endo = TRUE,
  adducts = c("H", "Na", "K", "NH4"), ppm_tol = 10
)

# ---- Build workbook ----
cat("\n--- Building Excel workbook ---\n")
wb_path <- build_workbook(spec, mets, STANDARD_DICT, validation,
                           ms_results = NULL, ms_info = NULL, opts = opts,
                           output_file = "inotersen_metabolite_library.xlsx")
cat("Workbook:", wb_path, "\n")
cat("File size:", file.size(wb_path), "bytes\n")

# Copy to the results dir (override with OLIGOMET_RESULTS_DIR; defaults to
# a session temp dir so this runs on any machine)
results_dir <- Sys.getenv("OLIGOMET_RESULTS_DIR", file.path(tempdir(), "results"))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
results_path <- file.path(results_dir, "inotersen_metabolite_library.xlsx")
file.copy(wb_path, results_path, overwrite = TRUE)
# R's file.copy may produce 0-byte on S3 mount; use cp instead
if (file.size(results_path) == 0) {
  system2("cp", args = c(wb_path, results_path))
}
cat("Copied to results:", results_path, "\n")
cat("Results file size:", file.size(results_path), "bytes\n")

# Verify workbook contents
cat("\n--- Verifying workbook ---\n")
library(openxlsx)
wb <- loadWorkbook(results_path)
cat("Sheets:", paste(names(wb), collapse = ", "), "\n")
for (sh in names(wb)) {
  df <- read.xlsx(results_path, sheet = sh, rows = 1:3, colNames = TRUE)
  cat(sprintf("  %s: %d cols, first row: %s\n", sh, ncol(df),
              paste(colnames(df)[1:min(4, ncol(df))], collapse=", ")))
}

# ---- Build report ----
cat("\n--- Building HTML report ---\n")
plot_dir <- file.path(tempdir(), "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
report_path <- build_report(spec, mets, STANDARD_DICT, validation,
                             ms_results = NULL, ms_info = NULL, opts = opts,
                             output_file = "inotersen_metabolite_report",
                             output_format = "html",
                             plot_dir = plot_dir)
cat("Report:", report_path, "\n")
if (file.exists(report_path)) {
  cat("Report file size:", file.size(report_path), "bytes\n")
  # Copy to results
  results_report <- file.path(results_dir, "inotersen_metabolite_report.html")
  file.copy(report_path, results_report, overwrite = TRUE)
  if (file.size(results_report) == 0) {
    system2("cp", args = c(report_path, results_report))
  }
  cat("Copied to results:", results_report, "\n")
}

# Copy plots to results too
for (f in list.files(plot_dir, pattern = "^plot_.*\\.png$", full.names = TRUE)) {
  dest <- file.path(results_dir, basename(f))
  file.copy(f, dest, overwrite = TRUE)
  if (file.size(dest) == 0) system2("cp", args = c(f, dest))
}
cat("\nPlots copied to results\n")

cat("\n==== Output generation complete ====\n")
