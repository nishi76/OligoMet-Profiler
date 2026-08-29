# =============================================================================
# run_batch_ms.R -- OligoMet Profiler
# Batch entry point: process many mzML/mzXML raw files in parallel (via the
# bundled Python deconvolution pipeline, inst/python/oligomet_deconv/),
# match the results against a theoretical metabolite library, confirm hits
# with MS2 fragment matching, and run a statistical comparison across
# experimental groups or a time series.
#
# Copy this file, edit the CONFIG block, and run:
#
#   Rscript run_batch_ms.R
#
# Requires Python 3.9+ on PATH with inst/python/requirements.txt installed
# (see that file, or README.md's "Batch processing" section).
#
# For a single sequence / single MS file, see run_custom_oligo.R instead.
# For the Shiny dashboard, see app.R.
#
# Author and developer: Nishikant Wase, PhD <nishikant.wase@gmail.com>
# Research Scientist, Thermo Fisher Scientific. An independent personal
# project, not a Thermo Fisher Scientific product.
# https://github.com/nishi76/OligoMet-Profiler -- MIT licence.
#
# FOR RESEARCH USE ONLY. Not for diagnostic, clinical, or regulatory
# submission use. Everything this pipeline reports is a computed
# prediction, not a measurement, and must be confirmed experimentally.
# Provided without warranty; the author accepts no liability for its use.
# Run oligomet_about(), or see DISCLAIMER.md, for the full statement.
# =============================================================================

## ---- Bootstrap: find modules and source them --------------------------------
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else "."
}, error = function(e) ".")

if (!file.exists(file.path(script_dir, "R", "chemistry_dict.R"))) {
  stop("Cannot locate pipeline modules (chemistry_dict.R not found). ",
       "Run this script from the package root, e.g. Rscript run_batch_ms.R.")
}

source(file.path(script_dir, "R", "about.R"))
source(file.path(script_dir, "R", "progress_utils.R"))
source(file.path(script_dir, "R", "chemistry_dict.R"))
source(file.path(script_dir, "R", "oligo_io.R"))
source(file.path(script_dir, "R", "metabolites.R"))
source(file.path(script_dir, "R", "mass_isotope.R"))
source(file.path(script_dir, "R", "fragments.R"))
source(file.path(script_dir, "R", "ms_matching.R"))
source(file.path(script_dir, "R", "batch_ms_processing.R"))
source(file.path(script_dir, "R", "statistics.R"))
source(file.path(script_dir, "R", "build_workbook.R"))

## ===========================================================================
## CONFIG -- edit this block for your experiment
## ===========================================================================

# --- Oligonucleotide (same notations as run_custom_oligo.R) -----------------
MY_TRIPLET <- "Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm"
custom_overrides <- list()

# --- Raw files ----------------------------------------------------------------
# mzML/mzXML files to process in parallel. One row of SAMPLE_META per file;
# `sample` must match the file's basename without extension.
BATCH_FILES <- list.files("raw_data", pattern = "\\.(mzML|mzXML)$",
                           full.names = TRUE, ignore.case = TRUE)

# --- Experimental design ------------------------------------------------------
# Provide EITHER a `group` column (two-group or multi-group comparison) OR a
# `timepoint` column (time-series trend), matching `sample` to BATCH_FILES.
# Example two-group design:
SAMPLE_META <- data.frame(
  sample = tools::file_path_sans_ext(basename(BATCH_FILES)),
  group = NA_character_,      # e.g. "control", "control", "treated", "treated"
  stringsAsFactors = FALSE
)
# Example time-series design instead:
# SAMPLE_META <- data.frame(
#   sample = tools::file_path_sans_ext(basename(BATCH_FILES)),
#   timepoint = c(0, 0, 1, 1, 2, 2),
#   stringsAsFactors = FALSE
# )

ANALYSIS_MODE <- "two_group"   # "two_group" | "multi_group" | "time_series"
GROUP_A <- "control"           # only used when ANALYSIS_MODE == "two_group"
GROUP_B <- "treated"

RUN_MS2_CONFIRMATION <- TRUE   # confirm MS1 hits with experimental MS2 spectra

# --- Pipeline parameters ------------------------------------------------------
PARAMS <- list(
  oligo_name    = "my_oligo",
  max_3p        = 10, max_5p = 10, endo = TRUE, min_frag_len = 3,
  z_range       = 3:12,          # theoretical charge envelope (matching)
  n_iso         = 5, max_oxid = 6, h_offset = 0, use_envipat = TRUE,
  ppm_tol       = 10,             # MS1 matching tolerance (ppm)
  adducts       = c("H", "Na", "K", "NH4"),
  frag_tol_ppm  = 25, frag_z_range = 1:3,   # MS2 confirmation tolerance
  # Python deconvolution parameters (see inst/python/oligomet_deconv/cli.py):
  deconv_z_range = 3:20, deconv_ppm_tol = 20, roi_ppm = 15, rt_tol = 0.15,
  min_intensity = 1e4, min_scans = 3, max_gap_scans = 2, ms2_watch_ppm = 50,
  n_workers     = NULL,           # NULL = os.cpu_count() - 1
  output_prefix = "my_oligo_batch",
  results_dir   = "results_batch"
)

## ===========================================================================
## END CONFIG -- below is the pipeline execution (no edits needed)
## ===========================================================================

run_batch_pipeline <- function() {
  cat("\n=============================================================\n")
  cat("  OligoMet Profiler\n")
  cat("  (batch MS processing mode)\n")
  cat("=============================================================\n")

  if (length(BATCH_FILES) == 0) stop("No mzML/mzXML files found -- check BATCH_FILES.")
  if (!dir.exists(PARAMS$results_dir)) dir.create(PARAMS$results_dir, recursive = TRUE)

  # Step 1: theoretical library
  dict <- build_dictionary(overrides = custom_overrides)
  spec <- parse_input(MY_TRIPLET, dict = dict)
  cat("  ", format_spec(spec), "\n")
  mets <- generate_metabolites(spec, opts = list(
    oligo_name = PARAMS$oligo_name, max_3p = PARAMS$max_3p, max_5p = PARAMS$max_5p,
    endo = PARAMS$endo, endo_sites = "all", min_frag_len = PARAMS$min_frag_len))
  cat("  Generated", length(mets), "metabolites\n")

  # Step 2: targeted MS2 watch-list (theoretical precursor m/z candidates)
  watchlist_path <- NULL
  if (RUN_MS2_CONFIRMATION) {
    watchlist_path <- tempfile(fileext = ".txt")
    write_precursor_watchlist(mets, dict, z_range = PARAMS$z_range,
                               max_oxid = PARAMS$max_oxid, h_offset = PARAMS$h_offset,
                               out_path = watchlist_path)
    cat("  Precursor watch-list written for targeted MS2 capture\n")
  }

  # Step 3: parallel deconvolution (Python)
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

  # Step 4: match against theoretical library + MS2 confirmation
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

  # Step 5: statistical comparison
  stats_results <- NULL
  has_meta <- nrow(batch_results$ms1_matches) > 0 &&
    all(SAMPLE_META$sample %in% unique(c(features$sample)))
  if (has_meta) {
    cat("\n--- Statistical comparison (", ANALYSIS_MODE, ") ---\n", sep = "")
    abund <- build_abundance_matrix(batch_results$ms1_matches)
    long <- abundance_long(abund, SAMPLE_META)
    stats_results <- switch(ANALYSIS_MODE,
      two_group = list(mode = "two_group", result = compare_two_groups(long, GROUP_A, GROUP_B)),
      multi_group = list(mode = "multi_group", result = compare_multi_groups(long)),
      time_series = list(mode = "time_series", result = compare_time_series(long)),
      stop("Unknown ANALYSIS_MODE: ", ANALYSIS_MODE))
    stats_csv <- file.path(PARAMS$results_dir, paste0(PARAMS$output_prefix, "_statistics.csv"))
    utils::write.csv(if (ANALYSIS_MODE == "multi_group") stats_results$result$omnibus else stats_results$result,
                      stats_csv, row.names = FALSE)
    cat("  Statistics table:", stats_csv, "\n")
  } else {
    cat("\n  Skipping statistics (fill in SAMPLE_META's group/timepoint column to enable)\n")
  }

  # Step 6: extend the workbook with batch/statistics sheets
  cat("\n--- Building Excel workbook ---\n")
  build_opts <- list(z_range = PARAMS$z_range, n_iso = PARAMS$n_iso,
                      max_oxid = PARAMS$max_oxid, h_offset = PARAMS$h_offset,
                      use_envipat = PARAMS$use_envipat, include_internal = TRUE,
                      max_3p = PARAMS$max_3p, max_5p = PARAMS$max_5p, endo = PARAMS$endo,
                      adducts = PARAMS$adducts, ppm_tol = PARAMS$ppm_tol)
  wb_file <- paste0(PARAMS$output_prefix, "_library.xlsx")
  wb_path <- build_workbook(spec, mets, dict, NULL, NULL, NULL, build_opts, wb_file,
                             output_dir = PARAMS$results_dir,
                             batch_ms_results = batch_results, stats_results = stats_results)
  cat("  Workbook saved:", wb_path, "\n")

  cat("\n=============================================================\n")
  cat("  Batch Pipeline Complete\n")
  cat("=============================================================\n")
  cat("  Files processed:", length(BATCH_FILES), "\n")
  cat("  Metabolites:", length(mets), "\n")
  cat("  MS1 matches:", nrow(batch_results$ms1_matches), "\n")
  cat("  Unmatched peaks:", nrow(batch_results$unmatched), "\n")
  cat("  Workbook:", wb_path, "\n")
  cat("=============================================================\n")

  invisible(list(spec = spec, mets = mets, dict = dict, features = features,
                 batch_results = batch_results, stats_results = stats_results,
                 workbook = wb_path))
}

## ---- Run --------------------------------------------------------------------
if (!interactive()) {
  run_batch_pipeline()
}
