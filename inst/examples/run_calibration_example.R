# =============================================================================
# run_calibration_example.R -- OligoMet Profiler
# Ready-to-run bundled example: a 5-level calibration curve (1, 5, 25, 100,
# 500 ng/mL), n=2 replicate injections per level, for the SAME inotersen
# reference sequence as inst/extdata/batch_example/ -- shipped under
# inst/extdata/calibration_example/ (see generate_calibration_example.py
# there for how they were built). No CONFIG editing needed -- just run:
#
#   Rscript inst/examples/run_calibration_example.R
#
# from the repository root. Requires Python 3.9+ with
# inst/python/requirements.txt installed -- see README.md.
#
# This exercises the full "Sample Type = standard + concentration ->
# calibration curve" path (fit_calibration_curve()/quantify_absolute() in
# R/statistics.R) end to end on a known-good synthetic dataset, so it
# doubles as a reference for what correctly set-up calibration standard
# metadata looks like, and as a way to confirm the pipeline itself is
# working before troubleshooting a real dataset that isn't producing a
# curve.
#
# Deconvolution parameters below are deliberately NOT the pipeline's own
# defaults (z_range = 3:12, not 3:20; min_charge_states = 3, not 2) --
# this dataset's parent envelope is exactly z=5,6,7,8 by construction, a
# consecutive small-integer charge set. Sweeping as wide as z=20 against
# only 4 closely related real peaks produces spurious "confirmed"
# 2-charge-state groups from simple-ratio harmonics (e.g. a z=6 peak
# reinterpreted at z=3, and a z=8 peak reinterpreted at z=4, coincide
# exactly since 6:3 and 8:4 are both 2:1) -- a real, general limitation of
# the sweep-based charge grouping (group_charge_states() in
# charge_group.py already documents the broader coincidental-collision
# issue), not a defect in this dataset. Narrowing z_range to the range
# actually expected for this known oligo, and requiring min_charge_states
# = 3 (comfortably below the real envelope's 4, comfortably above every
# harmonic collision's 2), cleanly recovers exactly one real feature per
# sample -- confirmed against this exact fixture before this script was
# written. A real, unknown analyte would instead want the pipeline's own
# wider defaults, tightening only if this same coincidental-collision
# pattern shows up in "Unidentified Peaks".
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
       "Rscript inst/examples/run_calibration_example.R")
}

for (.f in c("about.R", "progress_utils.R", "chemistry_dict.R", "oligo_io.R",
             "metabolites.R", "mass_isotope.R", "fragments.R", "ms_matching.R",
             "batch_ms_processing.R", "degradation.R", "statistics.R")) {
  source(file.path(repo_root, "R", .f))
}

example_dir <- tryCatch(
  system.file("extdata", "calibration_example", package = "OligoMetProfiler"),
  error = function(e) ""
)
if (!nzchar(example_dir) || !dir.exists(example_dir)) {
  example_dir <- file.path(repo_root, "inst", "extdata", "calibration_example")
}
if (!dir.exists(example_dir)) {
  stop("Cannot find the bundled example data (inst/extdata/calibration_example/).")
}

## ===========================================================================
## CONFIG -- pre-filled for the bundled example (nothing to edit)
## ===========================================================================

MY_TRIPLET <- INOTERSEN_TRIPLET  # same reference oligo as batch_example

BATCH_FILES <- list.files(example_dir, pattern = "\\.mzML$", full.names = TRUE)
SAMPLE_META <- utils::read.csv(file.path(example_dir, "sample_meta.csv"), stringsAsFactors = FALSE,
                                colClasses = "character")

CAL_WEIGHTING <- "1/x2"  # constant RELATIVE noise at every level -- see generate_calibration_example.py

PARAMS <- list(
  oligo_name = "inotersen_example", max_3p = 3, max_5p = 3, endo = FALSE, min_frag_len = 3,
  z_range = 3:12, n_iso = 0, max_oxid = 0, h_offset = 0, use_envipat = FALSE,
  ppm_tol = 10, adducts = c("H"),
  # See the module header comment above for why these two deconvolution
  # parameters are narrower than the pipeline's own defaults.
  deconv_z_range = 3:12, deconv_min_charge_states = 3,
  deconv_ppm_tol = 10, roi_ppm = 15, rt_tol = 0.15,
  min_intensity = 5000, min_scans = 3, max_gap_scans = 2,
  n_workers = 2, output_prefix = "inotersen_calibration_example", results_dir = "results_calibration_example"
)

## ===========================================================================
## END CONFIG
## ===========================================================================

run_calibration_pipeline <- function() {
  cat("\n=============================================================\n")
  cat("  OligoMet Profiler -- bundled calibration-curve example\n")
  cat("  (inotersen, 5 levels x n=2 replicates)\n")
  cat("=============================================================\n")

  if (length(BATCH_FILES) == 0) stop("No mzML files found under: ", example_dir)
  if (!dir.exists(PARAMS$results_dir)) dir.create(PARAMS$results_dir, recursive = TRUE)

  dict <- build_dictionary()
  spec <- parse_input(MY_TRIPLET, dict = dict)
  cat("  ", format_spec(spec), "\n")
  mets <- generate_metabolites(spec, opts = list(
    oligo_name = PARAMS$oligo_name, max_3p = PARAMS$max_3p, max_5p = PARAMS$max_5p,
    endo = PARAMS$endo, min_frag_len = PARAMS$min_frag_len),
    dict = dict)
  parent_id <- mets[[which(vapply(mets, function(m) identical(m$kind, "parent"), logical(1)))[1]]]$id
  cat("  Generated", length(mets), "metabolites; parent met_id =", parent_id, "\n")

  cat("\n--- Running parallel deconvolution on", length(BATCH_FILES), "file(s) ---\n")
  deconv <- run_batch_deconvolution(
    BATCH_FILES, output_dir = PARAMS$results_dir,
    output_file = paste0(PARAMS$output_prefix, "_features.tsv"),
    roi_ppm = PARAMS$roi_ppm, rt_tol = PARAMS$rt_tol, mass_tol_ppm = PARAMS$deconv_ppm_tol,
    z_range = PARAMS$deconv_z_range, min_intensity = PARAMS$min_intensity,
    min_scans = PARAMS$min_scans, max_gap_scans = PARAMS$max_gap_scans,
    min_charge_states = PARAMS$deconv_min_charge_states,
    n_workers = PARAMS$n_workers, progress = function(m) cat(" ", m, "\n"))
  features <- read_batch_features(deconv$features_path)
  cat("  Features extracted:", nrow(features), "across", length(unique(features$sample)), "samples\n")

  cat("\n--- MS1 matching ---\n")
  ms1_matches <- match_ms1_batch(mets, features, dict = dict, ppm_tol = PARAMS$ppm_tol,
                                  z_range = PARAMS$z_range, adducts = PARAMS$adducts,
                                  max_oxid = PARAMS$max_oxid, h_offset = PARAMS$h_offset)
  cat("  MS1 matches:", nrow(ms1_matches), "\n")

  cat("\n--- Calibration curve (", parent_id, ", weighting = ", CAL_WEIGHTING, ") ---\n", sep = "")
  curve <- fit_calibration_curve(ms1_matches, parent_id, SAMPLE_META, weighting = CAL_WEIGHTING)
  if (is.null(curve$model)) {
    cat("  Curve did not fit:", curve$note, "\n")
  } else {
    cat("  n_points:", curve$n_points, "\n")
    cat("  slope:", format(curve$slope, scientific = TRUE), "  intercept:", format(curve$intercept, scientific = TRUE), "\n")
    cat("  r_squared:", round(curve$r_squared, 5), "\n")
    cat("  concentration range:", paste(curve$conc_range, collapse = " - "), "ng/mL\n\n")
    print(curve$points[, c("sample", "concentration", "signal", "back_calc_concentration", "percent_re")])
  }

  cat("\n--- Absolute quantification (back-calculated concentrations) ---\n")
  quant <- quantify_absolute(ms1_matches, parent_id, SAMPLE_META, weighting = CAL_WEIGHTING)
  quant_csv <- file.path(PARAMS$results_dir, paste0(PARAMS$output_prefix, "_quantification.csv"))
  utils::write.csv(quant$quant, quant_csv, row.names = FALSE)
  cat("  Quantification table:", quant_csv, "\n")

  cat("\n=============================================================\n")
  cat("  Bundled Calibration Example Complete\n")
  cat("=============================================================\n")
  cat("  Files processed:", length(BATCH_FILES), "(5 levels x n=2)\n")
  cat("  Results directory:", normalizePath(PARAMS$results_dir), "\n")
  cat("=============================================================\n")

  invisible(list(spec = spec, mets = mets, dict = dict, features = features,
                 ms1_matches = ms1_matches, curve = curve, quant = quant))
}

if (!interactive()) {
  run_calibration_pipeline()
}
