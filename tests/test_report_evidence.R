# test_report_evidence.R -- per-metabolite LC-MS evidence for the study report
# (R/report_evidence.R), run end to end on the bundled raw-file batch example:
# deconvolution, matching, then re-reading each originating sample's mzML for
# its apex MS1 spectrum and extracted-ion chromatogram. Skips when Python or
# its packages (inst/python/requirements.txt) are not installed.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "progress_utils.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R", "ms_matching.R", "batch_ms_processing.R",
             "statistics.R", "degradation.R", "export_spectral.R", "mirror_plot.R", "xic_quant.R",
             "report_evidence.R")) {
  source(file.path(.pkg_root, "R", .f))
}
check <- function(label, ok) {
  cat(sprintf("%-78s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Report evidence suite ====\n\n")

py <- find_python()
py_ok <- !is.na(py) && identical(suppressWarnings(system2(py, c("-c", shQuote("import pyteomics, pandas")),
                                                            stdout = FALSE, stderr = FALSE)), 0L)
if (!py_ok) {
  cat("(skipped: Python with pyteomics/pandas not available)\n")
} else {
  ex <- file.path(.pkg_root, "inst", "extdata", "batch_example")
  spec <- parse_input(REFERENCE_OLIGOS$inotersen$triplet)
  mets <- generate_metabolites(spec, list(oligo_name = "inotersen"))
  wl <- write_precursor_watchlist(mets, STANDARD_DICT, out_path = tempfile(fileext = ".txt"))
  dec <- run_batch_deconvolution(normalizePath(list.files(ex, "\\.mzML$", full.names = TRUE)),
                                 precursor_watchlist = wl, n_workers = 1)
  feats <- read_batch_features(dec$features_path)
  ms2 <- read_batch_ms2(dec$ms2_path)
  res <- annotate_metabolites_batch(mets, feats, ms2)
  ids <- unique(res$ms1_matches$met_id[res$ms1_matches$kind %in% c("parent", "exo_3p")])
  ev <- collect_metabolite_evidence(ids, res$ms1_matches, feats, mets, STANDARD_DICT,
                                    ms2_confirmations = res$ms2_confirmations,
                                    ms2_spectra = res$ms2_spectra)
  p <- ev[[ids[1]]]
  check("evidence is collected for every requested metabolite", all(ids %in% names(ev)))
  check("originating sample is the strongest deprotonated observation",
        p$row$adduct == "H" && p$row$k_oxid == 0 &&
          p$row$intensity == max(res$ms1_matches$intensity[res$ms1_matches$met_id == ids[1] &
                                                             res$ms1_matches$adduct == "H" &
                                                             res$ms1_matches$k_oxid == 0]))
  check("theoretical isotopes are nominal peaks starting at the precursor m/z",
        abs(p$theo_iso$mz[1] - p$row$theo_mz) < 1e-6 && all(diff(p$theo_iso$iso) == 1))
  check("apex MS1 spectrum is re-read from the raw file and contains the precursor",
        !is.null(p$spectrum) && any(abs(p$spectrum$mz - p$row$obs_mz) / p$row$obs_mz * 1e6 < 10))
  check("extracted-ion chromatogram peaks at the matched retention time",
        !is.null(p$xic) && abs(p$xic$rt[which.max(p$xic$intensity)] - p$row$rt) < 0.2)
  check("parent carries its MS/MS mirror data", !is.null(p$ms2) && nrow(p$ms2$data$acquired) > 0)
  check("evidence plots build", inherits(plot_evidence_xic(p), "ggplot") &&
          inherits(plot_evidence_isotopes(p), "ggplot") && inherits(plot_evidence_envelope(p), "ggplot"))
  ev0 <- collect_metabolite_evidence(ids, res$ms1_matches, feats, mets, STANDARD_DICT, extract = FALSE)
  check("without raw-file extraction the record keeps matched values only",
        is.null(ev0[[ids[1]]]$spectrum) && !is.null(ev0[[ids[1]]]$theo_iso))
}

cat("\nAll report evidence checks passed.\n")
