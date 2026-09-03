# =============================================================================
# run_batch_example.R -- OligoMet Profiler
# Ready-to-run bundled example: 6 synthetic mzML files (3 "control" + 3
# "treated" replicates) for the inotersen reference sequence, shipped
# under inst/extdata/batch_example/ (see generate_example.py there for how
# they were built). No CONFIG editing needed -- just run:
#
#   Rscript inst/examples/run_batch_example.R
#
# from the repository root (or, if the package is installed,
# OligoMetProfiler::run_app() plus this file's own driver logic below
# still works via system.file()). Requires Python 3.9+ with
# inst/python/requirements.txt installed -- see README.md.
#
# This is the runnable counterpart to run_batch_ms.R: same pipeline, but
# pre-filled so a fresh checkout can exercise the whole batch workflow
# (parallel deconvolution -> MS1 matching -> retained unmatched peaks ->
# MS2 confirmation -> two-group statistics -> workbook) with a single
# command and no external data of your own.
# =============================================================================

## ---- Bootstrap: find modules and the bundled example data ------------------
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else "."
}, error = function(e) ".")

.find_repo_root <- function(start) {
  candidates <- c(start, file.path(start, "..", ".."), getwd(), file.path(getwd(), "..", ".."))
  for (d in candidates) if (file.exists(file.path(d, "R", "chemistry_dict.R"))) return(normalizePath(d))
  NULL
}
repo_root <- .find_repo_root(script_dir)
if (is.null(repo_root)) {
  stop("Cannot locate pipeline modules. Run this script from the repository root: ",
       "Rscript inst/examples/run_batch_example.R")
}

for (.f in c("about.R", "progress_utils.R", "chemistry_dict.R", "oligo_io.R",
             "metabolites.R", "mass_isotope.R", "fragments.R", "ms_matching.R",
             "batch_ms_processing.R", "statistics.R", "build_workbook.R")) {
  source(file.path(repo_root, "R", .f))
}

example_dir <- tryCatch(
  system.file("extdata", "batch_example", package = "OligoMetProfiler"),
  error = function(e) ""
)
if (!nzchar(example_dir) || !dir.exists(example_dir)) {
  example_dir <- file.path(repo_root, "inst", "extdata", "batch_example")
}
if (!dir.exists(example_dir)) {
  stop("Cannot find the bundled example data (inst/extdata/batch_example/).")
}

## ===========================================================================
## CONFIG -- pre-filled for the bundled example (nothing to edit)
## ===========================================================================

MY_TRIPLET <- INOTERSEN_TRIPLET
custom_overrides <- list()

BATCH_FILES <- list.files(example_dir, pattern = "\\.mzML$", full.names = TRUE)
SAMPLE_META <- utils::read.csv(file.path(example_dir, "sample_meta.csv"), stringsAsFactors = FALSE)

ANALYSIS_MODE <- "two_group"
GROUP_A <- "control"
GROUP_B <- "treated"
RUN_MS2_CONFIRMATION <- TRUE

PARAMS <- list(
  oligo_name = "inotersen_example", max_3p = 3, max_5p = 3, endo = FALSE, min_frag_len = 3,
  z_range = 3:12, n_iso = 0, max_oxid = 0, h_offset = 0, use_envipat = FALSE,
  ppm_tol = 10, adducts = c("H"), frag_tol_ppm = 15, frag_z_range = 1:2,
  deconv_z_range = 3:12, deconv_ppm_tol = 20, roi_ppm = 15, rt_tol = 0.15,
  min_intensity = 1e4, min_scans = 3, max_gap_scans = 2, ms2_watch_ppm = 50,
  n_workers = 2, output_prefix = "inotersen_batch_example", results_dir = "results_batch_example"
)

## ===========================================================================
## END CONFIG -- pipeline execution below is identical to run_batch_ms.R
## ===========================================================================

run_batch_pipeline <- function() {
  cat("\n=============================================================\n")
  cat("  OligoMet Profiler -- bundled batch example\n")
  cat("  (inotersen, 3 control + 3 treated replicates)\n")
  cat("=============================================================\n")

  if (length(BATCH_FILES) == 0) stop("No mzML files found under: ", example_dir)
  if (!dir.exists(PARAMS$results_dir)) dir.create(PARAMS$results_dir, recursive = TRUE)

  dict <- build_dictionary(overrides = custom_overrides)
  spec <- parse_input(MY_TRIPLET, dict = dict)
  cat("  ", format_spec(spec), "\n")
  mets <- generate_metabolites(spec, opts = list(
    oligo_name = PARAMS$oligo_name, max_3p = PARAMS$max_3p, max_5p = PARAMS$max_5p,
    endo = PARAMS$endo, endo_sites = "all", min_frag_len = PARAMS$min_frag_len))
  cat("  Generated", length(mets), "metabolites\n")

  watchlist_path <- NULL
  if (RUN_MS2_CONFIRMATION) {
    watchlist_path <- tempfile(fileext = ".txt")
    write_precursor_watchlist(mets, dict, z_range = PARAMS$z_range,
                               max_oxid = PARAMS$max_oxid, h_offset = PARAMS$h_offset,
                               out_path = watchlist_path)
    cat("  Precursor watch-list written for targeted MS2 capture\n")
  }

  cat("\n--- Running parallel deconvolution on", length(BATCH_FILES), "file(s) ---\n")
  deconv <- run_batch_deconvolution(
    BATCH_FILES, output_dir = PARAMS$results_dir,
    output_file = paste0(PARAMS$output_prefix, "_features.tsv"),
    precursor_watchlist = watchlist_path, ms2_watch_ppm = PARAMS$ms2_watch_ppm,
    roi_ppm = PARAMS$roi_ppm, rt_tol = PARAMS$rt_tol, mass_tol_ppm = PARAMS$deconv_ppm_tol,
    z_range = PARAMS$deconv_z_range, min_intensity = PARAMS$min_intensity,
    min_scans = PARAMS$min_scans, max_gap_scans = PARAMS$max_gap_scans,
    n_workers = PARAMS$n_workers, progress = function(m) cat(" ", m, "\n"))
  features <- read_batch_features(deconv$features_path)
  ms2_by_sample <- if (!is.null(deconv$ms2_path)) read_batch_ms2(deconv$ms2_path) else NULL
  cat("  Features extracted:", nrow(features), "across", length(unique(features$sample)), "samples\n")
  if (!is.null(deconv$profile_mode_files) && nrow(deconv$profile_mode_files) > 0) {
    cat("  WARNING:", nrow(deconv$profile_mode_files), "file(s) appear to be PROFILE mode",
        "(not centroided):", paste(deconv$profile_mode_files$sample, collapse = ", "), "\n")
  }

  cat("\n--- Matching + MS2 confirmation ---\n")
  batch_results <- annotate_metabolites_batch(
    mets, features, ms2_by_sample, dict = dict, ppm_tol = PARAMS$ppm_tol,
    z_range = PARAMS$z_range, adducts = PARAMS$adducts, max_oxid = PARAMS$max_oxid,
    h_offset = PARAMS$h_offset, n_iso = PARAMS$n_iso, use_envipat = PARAMS$use_envipat,
    frag_tol_ppm = PARAMS$frag_tol_ppm, frag_z_range = PARAMS$frag_z_range)
  cat("  MS1 matches:", nrow(batch_results$ms1_matches), "\n")
  cat("  Unmatched (retained) peaks:", nrow(batch_results$unmatched), "\n")
  cat("  MS2 confirmations:", nrow(batch_results$ms2_confirmations), "\n")

  unmatched_csv <- file.path(PARAMS$results_dir, paste0(PARAMS$output_prefix, "_unmatched_peaks.csv"))
  utils::write.csv(batch_results$unmatched, unmatched_csv, row.names = FALSE)
  cat("  Unmatched peaks:", unmatched_csv, "\n")

  cat("\n--- Statistical comparison (", ANALYSIS_MODE, ") ---\n", sep = "")
  abund <- build_abundance_matrix(batch_results$ms1_matches)
  long <- abundance_long(abund, SAMPLE_META)
  stats_results <- list(mode = "two_group", result = compare_two_groups(long, GROUP_A, GROUP_B))
  print(stats_results$result[, c("met_id", "met_name", "mean_a", "mean_b", "log2fc", "p_value", "p_adj")])
  stats_csv <- file.path(PARAMS$results_dir, paste0(PARAMS$output_prefix, "_statistics.csv"))
  utils::write.csv(stats_results$result, stats_csv, row.names = FALSE)
  cat("  Statistics table:", stats_csv, "\n")

  cat("\n--- Building Excel workbook ---\n")
  build_opts <- list(z_range = PARAMS$z_range, n_iso = PARAMS$n_iso,
                      max_oxid = PARAMS$max_oxid, h_offset = PARAMS$h_offset,
                      use_envipat = PARAMS$use_envipat, include_internal = FALSE,
                      max_3p = PARAMS$max_3p, max_5p = PARAMS$max_5p, endo = PARAMS$endo,
                      adducts = PARAMS$adducts, ppm_tol = PARAMS$ppm_tol)
  wb_file <- paste0(PARAMS$output_prefix, "_library.xlsx")
  wb_path <- build_workbook(spec, mets, dict, NULL, NULL, NULL, build_opts, wb_file,
                             output_dir = PARAMS$results_dir,
                             batch_ms_results = batch_results, stats_results = stats_results)
  cat("  Workbook saved:", wb_path, "\n")

  cat("\n=============================================================\n")
  cat("  Bundled Batch Example Complete\n")
  cat("=============================================================\n")
  cat("  Files processed:", length(BATCH_FILES), "(", paste(SAMPLE_META$group, collapse = ", "), ")\n")
  cat("  MS1 matches:", nrow(batch_results$ms1_matches), "\n")
  cat("  Results directory:", normalizePath(PARAMS$results_dir), "\n")
  cat("=============================================================\n")

  invisible(list(spec = spec, mets = mets, dict = dict, features = features,
                 batch_results = batch_results, stats_results = stats_results,
                 workbook = wb_path))
}

if (!interactive()) {
  run_batch_pipeline()
}
