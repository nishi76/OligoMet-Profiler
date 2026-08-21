# =============================================================================
# main_ION337.R -- OligoMet Profiler
# Reference/validation driver: reproduces the published ION337 metabolite
# library end to end and checks it against a known reference mass and
# (optionally) a reference workbook. Useful mainly as a self-test that the
# formula engine still matches a known result after code changes.
#
# For running the pipeline on your own sequence, use run_custom_oligo.R (a
# template you edit directly) or the Shiny app (app.R) instead -- this
# script's defaults are specific to the ION337 case. Passing --seq here
# still works (it swaps the input sequence) but Step 2's validation always
# checks against the ION337 reference mass, which will correctly report
# MISMATCH for any other sequence; that is expected, not a bug.
#
# Usage:
#   Rscript main_ION337.R                              # ION337 reference run
#   Rscript main_ION337.R --seq "Gm-sTm-sCm-..."        # different sequence
#                                                        #  (Step 2 will show MISMATCH -- expected)
#   Rscript main_ION337.R --ms data.mzML               # with MS data
#   Rscript main_ION337.R --format pdf                 # PDF report
#
# Grounded in:
#   - SynONIM (Lippens et al., JASMS 2024)
#   - Eluforsen metabolite profiling (Kim et al., Mol Ther Nucleic Acids 2019)
#   - OligoDistiller (Liu et al., Anal Chem 2025)
#   - FMVS automatic metabolite ID (Ye et al., J Chromatogr B 2025)
# =============================================================================

## ---- Module paths ---------------------------------------------------------
MODULE_DIR <- if (exists("MODULE_DIR")) MODULE_DIR else {
  # Try to find the module directory relative to this script
  script_dir <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else "."
  }, error = function(e) ".")
  script_dir
}
if (!file.exists(file.path(MODULE_DIR, "chemistry_dict.R"))) {
  stop("Cannot locate pipeline modules (chemistry_dict.R not found). ",
       "Run this script from the package root, e.g. ",
       "Rscript main_ION337.R, or set MODULE_DIR before sourcing.")
}

## ---- Source all modules ---------------------------------------------------
cat("Loading modules from", MODULE_DIR, "...\n")
source(file.path(MODULE_DIR, "progress_utils.R"))
source(file.path(MODULE_DIR, "chemistry_dict.R"))
source(file.path(MODULE_DIR, "oligo_io.R"))
source(file.path(MODULE_DIR, "metabolites.R"))
source(file.path(MODULE_DIR, "mass_isotope.R"))
source(file.path(MODULE_DIR, "fragments.R"))
source(file.path(MODULE_DIR, "ms_matching.R"))
source(file.path(MODULE_DIR, "build_workbook.R"))
source(file.path(MODULE_DIR, "build_report.R"))
source(file.path(MODULE_DIR, "export_acquisition.R"))

## ---- Parse command-line arguments -----------------------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    seq = NULL,           # oligonucleotide sequence (triplet/OD notation)
    ms_file = NULL,       # MS data file (mzML or peak list)
    ms_type = "ms1",      # MS data type: ms1, ms2, or both
    format = "html",      # report format: html or pdf
    output_prefix = NULL,  # output file prefix (default depends on --seq; see main())
    max_3p = 10,          # max 3' exonuclease truncations
    max_5p = 10,          # max 5' exonuclease truncations
    endo = TRUE,          # include endonuclease fragments
    max_oxid = 6,         # max PS->PO oxidation events
    z_range_min = 3,      # min charge state
    z_range_max = 12,     # max charge state
    ppm_tol = 10,         # mass tolerance in ppm
    adducts = "H,Na,K,NH4",  # adducts to check
    h_offset = 0,         # charge envelope offset (0=standard, 3.0046=ION337 WB)
    n_iso = 5,            # number of isotope peaks
    use_envipat = TRUE,   # use enviPat for isotope calculation
    validate_wb = NULL,   # path to reference workbook for validation
    results_dir = "results"           # created under the working directory
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else ""
    switch(key,
      "--seq" = { opts$seq <- val; i <- i + 2 },
      "--ms" = { opts$ms_file <- val; i <- i + 2 },
      "--ms-type" = { opts$ms_type <- val; i <- i + 2 },
      "--format" = { opts$format <- val; i <- i + 2 },
      "--output" = { opts$output_prefix <- val; i <- i + 2 },
      "--max-3p" = { opts$max_3p <- as.integer(val); i <- i + 2 },
      "--max-5p" = { opts$max_5p <- as.integer(val); i <- i + 2 },
      "--no-endo" = { opts$endo <- FALSE; i <- i + 1 },
      "--max-oxid" = { opts$max_oxid <- as.integer(val); i <- i + 2 },
      "--z-min" = { opts$z_range_min <- as.integer(val); i <- i + 2 },
      "--z-max" = { opts$z_range_max <- as.integer(val); i <- i + 2 },
      "--ppm" = { opts$ppm_tol <- as.numeric(val); i <- i + 2 },
      "--adducts" = { opts$adducts <- val; i <- i + 2 },
      "--h-offset" = { opts$h_offset <- as.numeric(val); i <- i + 2 },
      "--n-iso" = { opts$n_iso <- as.integer(val); i <- i + 2 },
      "--no-envipat" = { opts$use_envipat <- FALSE; i <- i + 1 },
      "--validate-wb" = { opts$validate_wb <- val; i <- i + 2 },
      "--results-dir" = { opts$results_dir <- val; i <- i + 2 },
      "--help" = {
        cat("Usage: Rscript main_ION337.R [options]\n",
            "  --seq <string>        Oligonucleotide sequence\n",
            "  --ms <file>            MS data file (mzML or peak list)\n",
            "  --format <html|pdf>    Report format\n",
            "  --output <prefix>      Output file prefix\n",
            "  --max-3p <int>         Max 3' truncations (default 10)\n",
            "  --max-5p <int>         Max 5' truncations (default 10)\n",
            "  --no-endo              Exclude endonuclease fragments\n",
            "  --max-oxid <int>       Max PS oxidation events (default 6)\n",
            "  --z-min <int>          Min charge state (default 3)\n",
            "  --z-max <int>          Max charge state (default 12)\n",
            "  --ppm <float>          Mass tolerance ppm (default 10)\n",
            "  --adducts <string>     Adducts (default H,Na,K,NH4)\n",
            "  --h-offset <float>     Envelope offset (default 0)\n",
            "  --n-iso <int>          Isotope peaks (default 5)\n",
            "  --no-envipat           Use built-in isotope calculation\n",
            "  --validate-wb <file>   Reference workbook for validation\n",
            sep = "")
        quit(status = 0)
      },
      { i <- i + 1 }  # skip unknown
    )
  }
  opts
}

## ---- End-to-end validation against ION337 workbook ------------------------
validate_vs_reference <- function(mets, dict, wb_path) {
  if (!file.exists(wb_path)) {
    cat("  Reference workbook not found:", wb_path, "\n")
    return(NULL)
  }
  library(openxlsx)
  ref <- read.xlsx(wb_path, sheet = 1, colNames = TRUE)
  # Find metabolite rows (Name not NA and Seq.Length not NA)
  name_rows <- which(!is.na(ref$Name) & !is.na(ref$Seq.Length))
  cat("  Reference workbook has", length(name_rows), "metabolites\n")

  results <- list()
  n_match <- 0; n_total <- 0; max_ppm <- 0

  for (i in name_rows) {
    ref_name <- as.character(ref$Name[i])
    ref_len <- as.integer(ref$Seq.Length[i])
    ref_formula <- as.character(ref$Formula[i])
    # Find matching metabolite by length and kind
    ref_kind <- if (grepl("^3 N_", ref_name)) "exo_3p" else
                if (grepl("^5 N_", ref_name)) "exo_5p" else
                if (grepl("Endo", ref_name)) "endo" else "parent"

    # Extract truncation number if present
    trunc_k <- NA
    if (grepl("N_(\\d+)", ref_name)) {
      trunc_k <- as.integer(sub(".*N_(\\d+).*", "\\1", ref_name))
    }

    # Find our metabolite
    our_met <- NULL
    if (ref_kind == "parent" && !is.na(ref_len)) {
      for (m in mets) {
        if (m$kind == "parent" && m$n == ref_len) { our_met <- m; break }
      }
    } else if (ref_kind == "exo_3p" && !is.na(trunc_k)) {
      for (m in mets) {
        if (m$kind == "exo_3p" && m$site == trunc_k) { our_met <- m; break }
      }
    } else if (ref_kind == "exo_5p" && !is.na(trunc_k)) {
      for (m in mets) {
        if (m$kind == "exo_5p" && m$site == trunc_k) { our_met <- m; break }
      }
    }

    if (is.null(our_met)) next
    n_total <- n_total + 1

    our_info <- metabolite_mass_info(our_met, dict)

    # Compare formula (by atom counts, order-independent)
    our_vec <- our_info$formula_vec
    ref_vec <- parse_formula(ref_formula)
    formula_match <- all(our_vec == ref_vec)

    # Compare mass (if available in workbook)
    # The workbook stores mass in charge-state rows, not the header row
    # So we compare our computed mass to the formula-derived mass
    ref_mass <- formula_mass(ref_vec, mono = TRUE)
    ppm <- abs(our_info$mono_mass - ref_mass) / ref_mass * 1e6

    if (formula_match && ppm < 1) n_match <- n_match + 1
    max_ppm <- max(max_ppm, ppm)

    results[[length(results) + 1]] <- data.frame(
      ref_name = ref_name, ref_len = ref_len,
      formula_match = formula_match, ppm = round(ppm, 4),
      our_mass = round(our_info$mono_mass, 4),
      ref_mass = round(ref_mass, 4),
      stringsAsFactors = FALSE)
  }

  res_df <- do.call(rbind, results)
  cat(sprintf("  Matched %d/%d metabolites (max ppm: %.4f)\n",
              n_match, n_total, max_ppm))
  list(results = res_df, n_match = n_match, n_total = n_total, max_ppm = max_ppm)
}

## ---- Main pipeline --------------------------------------------------------
main <- function() {
  opts <- parse_args()
  cat("\n")
  cat("=============================================================\n")
  cat("  OligoMet Profiler -- ION337 reference run\n")
  cat("=============================================================\n")

  prog <- progress_tracker(c(
    "Parsing input" = 1,
    "Chemistry engine self-test" = 1,
    "Generating metabolite library" = 2,
    "Computing masses and formulas" = 1,
    "Generating McLuckey MS/MS fragment ions" = 2,
    "Generating PRM inclusion list" = 1,
    "MS data import and matching" = 3,
    "End-to-end validation" = 1,
    "Building Excel workbook" = 15,
    "Building report" = 4,
    "Building Orbitrap Exploris acquisition method lists" = 2
  ))

  # Step 1: Parse input
  progress_next(prog)
  is_ion337 <- is.null(opts$seq)
  seq <- opts$seq %||% ION337_TRIPLET
  spec <- parse_input(seq)
  cat("  ", format_spec(spec), "\n")
  opts$output_prefix <- opts$output_prefix %||%
    (if (is_ion337) "ION337_metabolite" else "custom_oligo_metabolite")

  # Step 2: Chemistry engine self-test (always checks the ION337 reference
  # mass, regardless of the sequence being run here -- see header note)
  progress_next(prog)
  if (!is_ion337) {
    cat("  Note: --seq was supplied, so this check validates the formula\n",
        "  engine against ION337, not the sequence being run below.\n", sep = "")
  }
  validation <- validate_ION337(verbose = TRUE)

  # Step 3: Generate metabolite library
  progress_next(prog)
  met_opts <- list(
    oligo_name = if (is_ion337) "ION337" else opts$output_prefix,
    max_3p = opts$max_3p,
    max_5p = opts$max_5p,
    endo = opts$endo,
    endo_sites = "all",
    min_frag_len = 3
  )
  mets <- generate_metabolites(spec, opts = met_opts)
  cat("  Generated", length(mets), "metabolites\n")
  cat("  Parent:", mets[[1]]$name, "\n")
  cat("  3' truncations:", sum(sapply(mets, function(m) m$kind == "exo_3p")), "\n")
  cat("  5' truncations:", sum(sapply(mets, function(m) m$kind == "exo_5p")), "\n")
  cat("  Endonuclease:", sum(sapply(mets, function(m) grepl("endo", m$kind))), "\n")

  # Step 4: Compute masses and formulas
  progress_next(prog)
  parent_info <- metabolite_mass_info(mets[[1]])
  cat("  Parent formula:", parent_info$formula_str, "\n")
  cat("  Parent mono mass:", sprintf("%.6f Da", parent_info$mono_mass), "\n")
  cat("  Parent avg mass:", sprintf("%.4f Da", parent_info$avg_mass), "\n")

  # Step 5: Generate fragment ions
  progress_next(prog)
  parent_frags <- generate_fragments(mets[[1]])
  parent_ifrags <- generate_internal_fragments(mets[[1]])
  cat("  Terminal fragments:", length(parent_frags), "\n")
  cat("  Internal fragments:", length(parent_ifrags), "\n")

  # Step 6: Generate PRM inclusion list
  progress_next(prog)
  z_range <- opts$z_range_min:opts$z_range_max
  adducts <- strsplit(opts$adducts, ",")[[1]]
  prm <- prm_inclusion_list(mets, z_range = z_range, max_oxid = opts$max_oxid,
                             h_offset = opts$h_offset)
  cat("  PRM entries:", nrow(prm), "\n")
  cat("  m/z range:", sprintf("%.2f - %.2f\n", min(prm$precursor_mz), max(prm$precursor_mz)))

  # Step 7: Optional MS data import and matching
  progress_next(prog)
  ms_results <- NULL
  ms_info <- NULL
  if (!is.null(opts$ms_file) && file.exists(opts$ms_file)) {
    if (grepl("\\.mzML$|\\.mzml$|\\.mzXML$|\\.mzxml$", opts$ms_file, ignore.case = TRUE)) {
      ms_data <- parse_mzml(opts$ms_file)
      ms_info <- ms_data$info
      cat("  MS1 spectra:", ms_info$n_ms1, "\n")
      cat("  MS2 spectra:", ms_info$n_ms2, "\n")
      ms1_features <- extract_ms1_features(ms_data$ms1, ppm = opts$ppm_tol)
      cat("  MS1 features extracted:", nrow(ms1_features), "\n")
      ms_results <- annotate_metabolites(mets, ms1_features, ms_data$ms2,
        dict = ION337_DICT, ppm_tol = opts$ppm_tol, z_range = z_range,
        adducts = adducts, max_oxid = opts$max_oxid, h_offset = opts$h_offset,
        frag_tol_ppm = 25, frag_z_range = 1:3,
        n_iso = opts$n_iso, use_envipat = opts$use_envipat)
      cat("  MS1 matches:", nrow(ms_results$ms1_matches), "\n")
      cat("  Annotated metabolites:", nrow(ms_results$summary), "\n")
    } else {
      # Peak list
      ms1_features <- import_peak_list(opts$ms_file, type = opts$ms_type)
      cat("  Peak list features:", nrow(ms1_features), "\n")
      ms_results <- annotate_metabolites(mets, ms1_features, NULL,
        dict = ION337_DICT, ppm_tol = opts$ppm_tol, z_range = z_range,
        adducts = adducts, max_oxid = opts$max_oxid, h_offset = opts$h_offset,
        n_iso = opts$n_iso, use_envipat = opts$use_envipat)
      cat("  MS1 matches:", nrow(ms_results$ms1_matches), "\n")
    }
  } else {
    cat("  No MS data provided (library generation mode)\n")
  }

  # Step 8: Validate against reference workbook (only meaningful for the
  # ION337 case, or any workbook explicitly passed via --validate-wb; not
  # run automatically for other sequences since there's nothing to compare against)
  progress_next(prog)
  default_wb <- file.path(MODULE_DIR, "tests", "ION337_Workbook.xlsx")
  wb_path <- opts$validate_wb %||% Sys.getenv("ION337_WORKBOOK_PATH", default_wb)
  if (!is_ion337 && is.null(opts$validate_wb)) {
    cat("  Skipped (not the ION337 reference sequence; pass --validate-wb to compare anyway)\n")
    val_results <- NULL
  } else if (file.exists(wb_path)) {
    val_results <- validate_vs_reference(mets, ION337_DICT, wb_path)
  } else {
    cat("  No reference workbook found at", wb_path, "\n")
    val_results <- NULL
  }

  # Step 9: Build Excel workbook
  progress_next(prog)
  if (!dir.exists(opts$results_dir)) dir.create(opts$results_dir, recursive = TRUE)
  plot_dir <- tempdir()
  build_opts <- list(
    z_range = z_range, n_iso = opts$n_iso, max_oxid = opts$max_oxid,
    h_offset = opts$h_offset, use_envipat = opts$use_envipat,
    include_internal = TRUE, max_3p = opts$max_3p, max_5p = opts$max_5p,
    endo = opts$endo, adducts = adducts, ppm_tol = opts$ppm_tol
  )
  wb_file <- paste0(opts$output_prefix, "_library.xlsx")
  wb_path_out <- build_workbook(spec, mets, ION337_DICT, validation,
                                 ms_results, ms_info, build_opts, wb_file,
                                 console_tracker = prog)
  # Copy to results
  results_wb <- file.path(opts$results_dir, wb_file)
  system2("cp", args = c(wb_path_out, results_wb))
  cat("  Workbook saved:", results_wb, "\n")

  # Step 10: Build report
  progress_next(prog)
  report_file <- paste0(opts$output_prefix, "_report")
  report_path <- build_report(spec, mets, ION337_DICT, validation,
                               ms_results, ms_info, build_opts,
                               output_file = report_file,
                               output_format = opts$format,
                               plot_dir = plot_dir)
  # Copy to results
  ext <- if (opts$format == "pdf") ".pdf" else ".html"
  results_report <- file.path(opts$results_dir, paste0(report_file, ext))
  if (file.exists(report_path)) {
    system2("cp", args = c(report_path, results_report))
    cat("  Report saved:", results_report, "\n")
  }
  # Copy plots
  for (f in list.files(plot_dir, pattern = "^plot_.*\\.png$", full.names = TRUE)) {
    system2("cp", args = c(f, file.path(opts$results_dir, basename(f))))
  }
  cat("  Plots saved to results\n")

  # Step 11: Orbitrap Exploris acquisition method exports
  progress_next(prog)
  ms1_file <- file.path(opts$results_dir, paste0(opts$output_prefix, "_MS1_inclusion_list.csv"))
  ms2_file <- file.path(opts$results_dir, paste0(opts$output_prefix, "_MS2_PRM_target_list.csv"))
  frag_file <- file.path(opts$results_dir, paste0(opts$output_prefix, "_MS2_fragment_reference.csv"))
  utils::write.csv(
    thermo_ms1_inclusion_list(mets, ION337_DICT, z_range = z_range,
      h_offset = opts$h_offset, max_oxid = opts$max_oxid),
    ms1_file, row.names = FALSE)
  utils::write.csv(
    thermo_ms2_prm_target_list(mets, ION337_DICT, z_range = z_range,
      h_offset = opts$h_offset, max_oxid = opts$max_oxid),
    ms2_file, row.names = FALSE)
  utils::write.csv(ms2_fragment_reference(mets, ION337_DICT, z_range = 1:3),
                   frag_file, row.names = FALSE)
  cat("  MS1 inclusion list:", ms1_file, "\n")
  cat("  MS2 PRM target list:", ms2_file, "\n")
  cat("  MS2 fragment reference:", frag_file, "\n")

  progress_finish(prog)

  # Summary
  cat("\n=============================================================\n")
  cat("  Pipeline Complete\n")
  cat("=============================================================\n")
  cat("  Metabolites generated:", length(mets), "\n")
  cat("  Fragment ions (parent):", length(parent_frags) + length(parent_ifrags), "\n")
  cat("  PRM inclusion entries:", nrow(prm), "\n")
  if (!is.null(val_results)) {
    cat("  Validation:", val_results$n_match, "/", val_results$n_total,
        "metabolites matched (max", sprintf("%.4f ppm)\n", val_results$max_ppm))
  }
  cat("  Workbook:", results_wb, "\n")
  cat("  Report:", results_report, "\n")
  cat("=============================================================\n")

  invisible(list(
    spec = spec, mets = mets, validation = validation,
    prm = prm, ms_results = ms_results, val_results = val_results,
    workbook = results_wb, report = results_report
  ))
}

## ---- Run ------------------------------------------------------------------
if (!interactive()) {
  main()
}
