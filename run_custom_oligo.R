# =============================================================================
# run_custom_oligo.R -- OligoMet Profiler
# General-purpose entry point for running the pipeline on any oligonucleotide.
# Works out of the box with standard chemistry (DNA/RNA/2'OMe/2'F/MOE/cEt/LNA,
# PO/PS linkages) and supports custom chemistry overrides for anything else.
# Copy this file, edit the CONFIG block, and run:
#
#   Rscript run_custom_oligo.R
#
# This demonstrates three things:
#   1. Standard chemistry (uses built-in dictionary codes)
#   2. Custom chemistry overrides (add your own base/sugar/linkage/conjugate)
#   3. All three input notations (triplet, OligoDistiller, structured)
#
# For the Shiny dashboard instead, see app.R. To self-test that the
# formula engine still matches a known mass, run validate_reference() from R
# or the validation scripts under tests/.
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
# When run via Rscript, detect the script's own directory.
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else "."
}, error = function(e) ".")

if (!file.exists(file.path(script_dir, "R", "chemistry_dict.R"))) {
  stop("Cannot locate pipeline modules (chemistry_dict.R not found). ",
       "Run this script from the package root, e.g. ",
       "Rscript run_custom_oligo.R.")
}

source(file.path(script_dir, "R", "about.R"))
source(file.path(script_dir, "R", "default_params.R"))
source(file.path(script_dir, "R", "progress_utils.R"))
source(file.path(script_dir, "R", "chemistry_dict.R"))
source(file.path(script_dir, "R", "oligo_io.R"))
source(file.path(script_dir, "R", "metabolites.R"))
source(file.path(script_dir, "R", "mass_isotope.R"))
source(file.path(script_dir, "R", "fragments.R"))
source(file.path(script_dir, "R", "ms_matching.R"))
source(file.path(script_dir, "R", "build_workbook.R"))
source(file.path(script_dir, "R", "build_report.R"))
source(file.path(script_dir, "R", "export_acquisition.R"))
source(file.path(script_dir, "R", "export_spectral.R"))

## ===========================================================================
## CONFIG -- edit this block for your oligonucleotide
## ===========================================================================

# --- Option A: Triplet notation ---------------------------------------------
# Format: [linkage][base][sugar] per token, dash-separated.
#   First token has NO linkage prefix (5' terminal).
#   Linkage prefix on token i = incoming bond (from position i-1).
#   Bases:  A G C T U S(m5C) D(2,6-diaminopurine) I(inosine)
#   Sugars: d(DNA) r(RNA) m(2'OMe) f(2'F) e or MOE (2'-MOE) cEt LNA
#   Linkages: o or p (PO), s (PS), u (PS stereochem variant)
#
# Example: 20-mer gapmer ASO, 2'OMe wings (positions 1-5, 16-20),
#          DNA gap (positions 6-15), all PS linkages.
#          Sequence: G m T m C m T m C m  T d C d T d C d T d T d  C m T m C m T m G m
#          Triplet:  Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm

MY_TRIPLET <- "Gm-sTm-sCm-sTm-sCm-sTd-sCd-sTd-sCd-sTd-sTd-sCm-sTm-sCm-sTm-sGm"

# --- Option B: OligoDistiller notation --------------------------------------
# Format: OH-<base><sugar>[*]-...-<base><sugar>-OH
#   * = PS outgoing bond; no * = PO.
#   Sugar suffix: m=2'OMe, f=2'F, d=DNA, r=RNA, e=MOE
#
# Same gapmer in OligoDistiller notation (all PS, so every token gets *):
# MY_SEQ <- "OH-Gm*-Tm*-Cm*-Tm*-Cm*-Td*-Cd*-Td*-Cd*-Td*-Td*-Cm*-Tm*-Cm*-Tm*-Gm*-OH"

# --- Option C: Structured input (most flexible) -----------------------------
# Use this when you need conjugates or non-standard codes.
#
# MY_SPEC <- list(
#   bases    = c("G","T","C","T","C", "T","C","T","C","T", "T","C","T","C","T","G"),
#   sugars   = c("m","m","m","m","m", "d","d","d","d","d", "d","m","m","m","m","m"),
#   linkages = c("s","s","s","s","s","s","s","s","s","s","s","s","s","s","s", NA),
#   conj5    = "none",
#   conj3    = "none"      # or "GalNAc", "cholesterol", "C6", "TEG", etc.
# )

# --- Custom chemistry overrides ---------------------------------------------
# Add or replace dictionary entries. Each override is a named list with
# 'formula' (compact named vector or formula string) and optional 'name'.
#
# Examples:
#   custom_overrides <- list(
#     # Add a custom base
#     X = list(formula = c(C=6,H=7,N=5,O=1), name = "custom modified base X"),
#     # Add a custom sugar
#     Q = list(formula = c(C=10,H=20,O=7), name = "custom 2'-mod sugar Q"),
#     # Override an existing entry with a corrected formula
#     cEt = list(formula = "C7H12O4", name = "constrained ethyl (confirmed)"),
#     # Add a custom conjugate
#     myLinker = list(formula = c(C=12,H=24,N=2,O=4), name = "custom linker",
#                     attach = "replace_H")
#   )
#
# Leave empty if using only standard dictionary codes:
custom_overrides <- list()

# --- Pipeline parameters ----------------------------------------------------
# Built on top of DEFAULT_PIPELINE_PARAMS (R/default_params.R) -- the same
# defaults the Shiny dashboard's sidebar uses -- so the two interfaces can't
# quietly drift apart. Override only what this run needs to change; anything
# left out keeps the shared default.
PARAMS <- modifyList(DEFAULT_PIPELINE_PARAMS, list(
  oligo_name     = "my_oligo",          # used in output filenames and labels
  z_range        = with(DEFAULT_PIPELINE_PARAMS, z_min:z_max),  # charge state range for envelopes
  output_prefix  = "my_oligo_metabolite",
  results_dir    = "results"            # created under the working directory
))

# --- MS data (optional) -----------------------------------------------------
# Set to a file path to run MS matching, or NULL for library-only mode.
# Accepted: .mzML, .mzXML, or peak-list text files (mz, intensity, charge).
MS_FILE <- NULL
# MS_FILE <- "/path/to/data.mzML"

## ===========================================================================
## END CONFIG -- below is the pipeline execution (no edits needed)
## ===========================================================================

run_pipeline <- function() {
  cat("\n=============================================================\n")
  cat("  OligoMet Profiler\n")
  cat("  (custom oligo mode)\n")
  cat("=============================================================\n")

  if (!dir.exists(PARAMS$results_dir)) dir.create(PARAMS$results_dir, recursive = TRUE)
  plot_dir <- tempdir()

  prog <- progress_tracker(c(
    "Building chemistry dictionary" = 1,
    "Parsing input" = 1,
    "Parent formula and mass" = 1,
    "Generating metabolite library" = 2,
    "Generating McLuckey MS/MS fragment ions" = 2,
    "Generating PRM inclusion list" = 1,
    "MS data import and matching" = 3,
    "Building Excel workbook" = 15,
    "Building report" = 4,
    "Building Orbitrap Exploris acquisition method lists" = 2
  ))

  # Step 1: Build dictionary with custom overrides
  progress_next(prog)
  dict <- build_dictionary(overrides = custom_overrides)
  n_custom <- length(custom_overrides)
  cat("  Dictionary entries:", length(dict), "\n")
  if (n_custom > 0) {
    cat("  Custom overrides applied:", n_custom, "\n")
    for (nm in names(custom_overrides)) {
      e <- dict[[nm]]
      cat("    ", nm, " -> ", format_formula(e$formula),
          " (", e$kind, ")\n", sep = "")
    }
  }

  # Step 2: Parse input (auto-detects notation)
  progress_next(prog)
  # Use whichever input is configured above:
  if (exists("MY_SPEC") && !is.null(MY_SPEC)) {
    spec <- parse_input(MY_SPEC, dict = dict)
  } else {
    seq_input <- MY_TRIPLET %||% ""
    spec <- parse_input(seq_input, dict = dict)
  }
  cat("  ", format_spec(spec), "\n")

  # Step 3: Compute parent formula and mass (sanity check)
  progress_next(prog)
  parent_f <- assemble_oligo_formula(spec$bases, spec$sugars, spec$linkages,
                                      spec$conj5, spec$conj3, dict = dict)
  parent_mass <- formula_mass(parent_f, mono = TRUE)
  parent_avg <- formula_mass(parent_f, mono = FALSE)
  cat("  Formula:", format_formula(parent_f), "\n")
  cat("  Mono mass:", sprintf("%.6f Da\n", parent_mass))
  cat("  Avg mass:", sprintf("%.4f Da\n", parent_avg))
  cat("  Length:", spec$n, "nucleotides\n")

  # Step 4: Generate metabolite library
  progress_next(prog)
  met_opts <- list(
    oligo_name = PARAMS$oligo_name,
    max_3p = PARAMS$max_3p,
    max_5p = PARAMS$max_5p,
    endo = PARAMS$endo,
    endo_sites = PARAMS$endo_sites,
    min_frag_len = PARAMS$min_frag_len
  )
  mets <- generate_metabolites(spec, opts = met_opts)
  cat("  Generated", length(mets), "metabolites\n")
  cat("  Parent:", mets[[1]]$name, "\n")
  cat("  3' truncations:", sum(sapply(mets, function(m) m$kind == "exo_3p")), "\n")
  cat("  5' truncations:", sum(sapply(mets, function(m) m$kind == "exo_5p")), "\n")
  cat("  Endonuclease:", sum(sapply(mets, function(m) grepl("endo", m$kind))), "\n")

  # Step 5: Generate fragment ions
  progress_next(prog)
  parent_frags <- generate_fragments(mets[[1]], dict = dict)
  parent_ifrags <- generate_internal_fragments(mets[[1]], dict = dict)
  cat("  Terminal fragments:", length(parent_frags), "\n")
  cat("  Internal fragments:", length(parent_ifrags), "\n")

  # Step 6: Generate PRM inclusion list
  progress_next(prog)
  prm <- prm_inclusion_list(mets, z_range = PARAMS$z_range,
                             max_oxid = PARAMS$max_oxid,
                             h_offset = PARAMS$h_offset, dict = dict)
  cat("  PRM entries:", nrow(prm), "\n")
  cat("  m/z range:", sprintf("%.2f - %.2f\n",
      min(prm$precursor_mz), max(prm$precursor_mz)))

  # Step 7: Optional MS matching
  progress_next(prog)
  ms_results <- NULL
  ms_info <- NULL
  if (!is.null(MS_FILE) && file.exists(MS_FILE)) {
    if (grepl("\\.mzML$|\\.mzml$|\\.mzXML$", MS_FILE, ignore.case = TRUE)) {
      ms_data <- parse_mzml(MS_FILE)
      ms_info <- ms_data$info
      ms1_features <- extract_ms1_features(ms_data$ms1, ppm = PARAMS$ppm_tol)
      ms_results <- annotate_metabolites(mets, ms1_features, ms_data$ms2,
        dict = dict, ppm_tol = PARAMS$ppm_tol, z_range = PARAMS$z_range,
        adducts = PARAMS$adducts, max_oxid = PARAMS$max_oxid,
        h_offset = PARAMS$h_offset, n_iso = PARAMS$n_iso,
        use_envipat = PARAMS$use_envipat)
      cat("  MS1 matches:", nrow(ms_results$ms1_matches), "\n")
    } else {
      ms1_features <- import_peak_list(MS_FILE)
      ms_results <- annotate_metabolites(mets, ms1_features, NULL,
        dict = dict, ppm_tol = PARAMS$ppm_tol, z_range = PARAMS$z_range,
        adducts = PARAMS$adducts, max_oxid = PARAMS$max_oxid,
        h_offset = PARAMS$h_offset, n_iso = PARAMS$n_iso,
        use_envipat = PARAMS$use_envipat)
      cat("  MS1 matches:", nrow(ms_results$ms1_matches), "\n")
    }
  } else {
    cat("  No MS data provided (library generation mode)\n")
  }

  # Step 8: Build Excel workbook
  progress_next(prog)
  build_opts <- list(
    z_range = PARAMS$z_range, n_iso = PARAMS$n_iso,
    max_oxid = PARAMS$max_oxid, h_offset = PARAMS$h_offset,
    use_envipat = PARAMS$use_envipat, include_internal = FALSE,
    max_3p = PARAMS$max_3p, max_5p = PARAMS$max_5p,
    endo = PARAMS$endo, adducts = PARAMS$adducts, ppm_tol = PARAMS$ppm_tol
  )
  wb_file <- paste0(PARAMS$output_prefix, "_library.xlsx")
  wb_path <- build_workbook(spec, mets, dict, NULL, ms_results, ms_info,
                             build_opts, wb_file, console_tracker = prog)
  results_wb <- file.path(PARAMS$results_dir, wb_file)
  system2("cp", args = c(wb_path, results_wb))
  cat("  Workbook saved:", results_wb, "\n")

  # Step 9: Build report
  progress_next(prog)
  report_file <- paste0(PARAMS$output_prefix, "_report")
  report_path <- build_report(spec, mets, dict, NULL, ms_results, ms_info,
                               build_opts, output_file = report_file,
                               output_format = "html", plot_dir = plot_dir)
  results_report <- file.path(PARAMS$results_dir, paste0(report_file, ".html"))
  if (file.exists(report_path)) {
    system2("cp", args = c(report_path, results_report))
    cat("  Report saved:", results_report, "\n")
  }
  for (f in list.files(plot_dir, pattern = "^plot_.*\\.png$", full.names = TRUE)) {
    system2("cp", args = c(f, file.path(PARAMS$results_dir, basename(f))))
  }
  cat("  Plots saved to results\n")

  # Step 10: Orbitrap Exploris acquisition method exports
  progress_next(prog)
  ms1_file <- file.path(PARAMS$results_dir, paste0(PARAMS$output_prefix, "_MS1_inclusion_list.csv"))
  ms2_file <- file.path(PARAMS$results_dir, paste0(PARAMS$output_prefix, "_MS2_PRM_target_list.csv"))
  frag_file <- file.path(PARAMS$results_dir, paste0(PARAMS$output_prefix, "_MS2_fragment_reference.csv"))
  ms2_z_range <- PARAMS$ms2_z_min:PARAMS$ms2_z_max
  utils::write.csv(
    thermo_ms1_inclusion_list(mets, dict, z_range = PARAMS$z_range,
      h_offset = PARAMS$h_offset, max_oxid = PARAMS$max_oxid,
      rt_end = PARAMS$method_length, max_targets = PARAMS$ms1_target_cap),
    ms1_file, row.names = FALSE)
  utils::write.csv(
    thermo_ms2_prm_target_list(mets, dict, z_range = ms2_z_range,
      h_offset = PARAMS$h_offset, max_oxid = PARAMS$max_oxid,
      rt_end = PARAMS$method_length, nce = PARAMS$hcd_nce,
      max_targets = PARAMS$ms2_target_cap),
    ms2_file, row.names = FALSE)
  utils::write.csv(ms2_fragment_reference(mets, dict, z_range = ms2_z_range),
                   frag_file, row.names = FALSE)
  cat("  MS1 inclusion list:", ms1_file, "\n")
  cat("  MS2 PRM target list:", ms2_file, "\n")
  cat("  MS2 fragment reference:", frag_file, "\n")

  progress_finish(prog)

  # Summary
  cat("\n=============================================================\n")
  cat("  Pipeline Complete\n")
  cat("=============================================================\n")
  cat("  Oligo:", PARAMS$oligo_name, "\n")
  cat("  Formula:", format_formula(parent_f), "\n")
  cat("  Mono mass:", sprintf("%.6f Da\n", parent_mass))
  cat("  Metabolites:", length(mets), "\n")
  cat("  Fragment ions:", length(parent_frags) + length(parent_ifrags), "\n")
  cat("  PRM entries:", nrow(prm), "\n")
  cat("  Workbook:", results_wb, "\n")
  cat("  Report:", results_report, "\n")
  cat("=============================================================\n")

  invisible(list(spec = spec, mets = mets, dict = dict, prm = prm,
                 ms_results = ms_results, workbook = results_wb,
                 report = results_report, ms1_inclusion = ms1_file,
                 ms2_prm_targets = ms2_file, fragment_reference = frag_file))
}

## ---- Run ------------------------------------------------------------------
if (!interactive()) {
  run_pipeline()
}
