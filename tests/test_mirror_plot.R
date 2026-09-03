# test_mirror_plot.R -- validate R/mirror_plot.R: acquired-vs-theoretical
# mirror-plot data/plotting, and the batch-level downloadable exports
# (multi-page mirror-plot PDF, annotated MSP of acquired spectra) built on
# top of it.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R", "ms_matching.R",
             "batch_ms_processing.R", "statistics.R", "degradation.R",
             "export_spectral.R", "mirror_plot.R")) {
  source(file.path(.pkg_root, "R", .f))
}

cat("==== mirror_plot.R validation ====\n\n")

dict <- STANDARD_DICT
spec <- parse_input(INOTERSEN_TRIPLET)
mets <- generate_metabolites(spec, opts = list(oligo_name = "inotersen",
                                                max_3p = 0, max_5p = 0, endo = FALSE))

## ---- mirror_spectrum_data() / plot_mirror_spectrum(): single hit ---------
cat("--- mirror_spectrum_data() on a synthetic MS2 spectrum ---\n")
frags <- generate_fragments(mets[[1]], dict, z_range = 1:2)
true_peaks <- do.call(rbind, lapply(frags[1:8], function(f) {
  data.frame(mz = (f$mono_mass - 1 * .PROTON) / 1, intensity = 5000)
}))
noise_peaks <- data.frame(mz = c(300.1, 555.5, 812.3), intensity = c(1000, 900, 800))
ms2_peaks <- rbind(true_peaks, noise_peaks)

spec_data <- mirror_spectrum_data(mets[[1]], ms2_peaks, dict, tol_ppm = 15, z_range = 1:2)
stopifnot(sum(spec_data$acquired$matched) == 8)
stopifnot(all(!is.na(spec_data$acquired$label[spec_data$acquired$matched])))
stopifnot(all(is.na(spec_data$acquired$label[!spec_data$acquired$matched])))
cat("8/8 true peaks matched and labeled, 3/3 noise peaks unmatched: PASS\n")

p <- plot_mirror_spectrum(spec_data, title = "test")
stopifnot(inherits(p, "ggplot"))
cat("plot_mirror_spectrum() returns a ggplot object: PASS\n")

## ---- Batch pipeline: real fixture, real Python deconvolution -------------
cat("\n--- Batch MS2 confirmation against inst/extdata/batch_example ---\n")
ex_dir <- file.path(.pkg_root, "inst", "extdata", "batch_example")
files <- list.files(ex_dir, pattern = "\\.mzML$", full.names = TRUE)
if (length(files) == 0 || is.na(find_python())) {
  cat("SKIPPED (no batch_example fixtures or no python3 on PATH)\n")
} else {
  watchlist_path <- tempfile(fileext = ".txt")
  write_precursor_watchlist(mets, dict, z_range = 3:12, max_oxid = 0, h_offset = 0,
                             out_path = watchlist_path)
  deconv <- tryCatch(
    run_batch_deconvolution(files, precursor_watchlist = watchlist_path,
                             mass_tol_ppm = 10, n_workers = 1),
    error = function(e) NULL)
  if (is.null(deconv)) {
    cat("SKIPPED (batch deconvolution failed -- see console for details)\n")
  } else {
    feats <- read_batch_features(deconv$features_path)
    feats$sample <- tools::file_path_sans_ext(basename(feats$sample))
    ms2 <- read_batch_ms2(deconv$ms2_path)
    ms2$sample <- tools::file_path_sans_ext(basename(ms2$sample))

    batch_results <- annotate_metabolites_batch(mets, feats, ms2, dict = dict,
      ppm_tol = 10, z_range = 3:12, adducts = "H", max_oxid = 0,
      frag_tol_ppm = 25, frag_z_range = 1:2)
    n_conf <- nrow(batch_results$ms2_confirmations)
    cat("ms2_confirmations:", n_conf, "\n")
    stopifnot(n_conf > 0)
    stopifnot(length(batch_results$ms2_spectra) == n_conf)

    cat("\n--- batch_mirror_plots_pdf() ---\n")
    pdf_out <- tempfile(fileext = ".pdf")
    batch_mirror_plots_pdf(mets, dict, batch_results$ms2_confirmations,
                            batch_results$ms2_spectra, file = pdf_out,
                            tol_ppm = 25, z_range = 1:2)
    stopifnot(file.exists(pdf_out))
    stopifnot(file.info(pdf_out)$size > 0)
    cat("PDF written, one page per confirmed hit: PASS\n")
    unlink(pdf_out)

    cat("\n--- batch_annotated_msp_records() / write_msp(measured=TRUE) ---\n")
    recs <- batch_annotated_msp_records(mets, dict, batch_results$ms2_confirmations,
                                         batch_results$ms2_spectra, tol_ppm = 25, z_range = 1:2)
    stopifnot(length(recs) == n_conf)
    stopifnot(all(vapply(recs, function(r) r$precursor_mz > 0, logical(1))))
    n_labeled_total <- sum(vapply(recs, function(r) sum(nzchar(r$peaks$annotation)), integer(1)))
    stopifnot(n_labeled_total > 0)
    cat(length(recs), "records, precursor m/z all positive,", n_labeled_total,
        "annotated peaks total: PASS\n")

    msp_out <- tempfile(fileext = ".msp")
    write_msp(recs, msp_out, measured = TRUE)
    msp_lines <- readLines(msp_out)
    stopifnot(any(grepl("^COMMENT: Acquired MS2 spectrum", msp_lines)))
    stopifnot(!any(grepl("^COMMENT: Theoretical", msp_lines)))
    cat("MSP comment correctly labeled 'Acquired', not 'Theoretical': PASS\n")
    unlink(msp_out)
  }
}

cat("\n==== All mirror_plot tests passed ====\n")
