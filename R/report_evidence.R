# =============================================================================
# report_evidence.R
# Per-metabolite LC-MS evidence for the regulatory report: for each reported
# metabolite, the originating sample (the study sample with the strongest
# deprotonated, non-desulfurized signal), its precursor ion and charge-state
# envelope, the theoretical isotope pattern at the dominant charge state, the
# acquired MS1 spectrum and extracted-ion chromatogram re-read from that
# sample's raw file, and the MS/MS mirror plot when a confirmation exists.
#
# The raw-file read is a single streaming pass per file through the bundled
# Python reader (inst/python/oligomet_deconv/evidence_cli.py); without Python
# or the raw files the evidence falls back to the matched values and the
# theoretical pattern.
# =============================================================================

.ev_ion_label <- function(z, adduct = "H") {
  if (is.null(adduct) || is.na(adduct) || identical(adduct, "H"))
    return(sprintf("[M-%dH]%d-", z, z))
  sprintf("[M-%dH+%s]%d-", z + 1L, adduct, z)
}

# Strongest deprotonated, intact-backbone observation of a metabolite in the
# study samples; falls back to any adduct / desulfurization level.
.ev_best_row <- function(ms1, met_id, samples) {
  d <- ms1[ms1$met_id == met_id & ms1$sample %in% samples, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  sig <- if ("intensity" %in% names(d)) d$intensity else d$area
  pref <- d$adduct == "H" & d$k_oxid == 0
  if (any(pref, na.rm = TRUE)) { d <- d[which(pref), , drop = FALSE]; sig <- sig[which(pref)] }
  d[which.max(sig), , drop = FALSE]
}

#' Collect per-metabolite LC-MS evidence for the study report
#'
#' @param met_ids Metabolite ids to characterize (in report order).
#' @param ms1_matches Batch MS1 match table (see [annotate_metabolites_batch()]).
#' @param features Batch feature table (see [read_batch_features()]); its
#'   `source_file` column locates each sample's raw file.
#' @param mets Metabolite library.
#' @param dict Chemistry dictionary.
#' @param study_samples Samples eligible as the originating sample.
#' @param ms2_confirmations,ms2_spectra MS/MS confirmation table and its
#'   spectra (see [annotate_metabolites_batch()]), or `NULL`.
#' @param h_offset,use_envipat Mass-calculation settings used for the run.
#' @param ppm Extraction tolerance for the chromatogram (ppm).
#' @param rt_half_window Chromatogram half-width around the apex (min).
#' @param extract Re-read the raw files for spectra and chromatograms.
#' @param python_bin,module_dir Python interpreter and module path.
#' @return A named list (by `met_id`) of evidence records: `row` (the
#'   best match), `envelope` (charge states in the originating sample),
#'   `theo_iso` (theoretical isotope cluster at the precursor charge),
#'   `spectrum` and `xic` (acquired data, or `NULL`), `ms2` (mirror-plot
#'   data and confirmation row, or `NULL`), `file`, `n_detected`.
#' @export
collect_metabolite_evidence <- function(met_ids, ms1_matches, features = NULL, mets = NULL,
                                        dict = STANDARD_DICT, study_samples = NULL,
                                        ms2_confirmations = NULL, ms2_spectra = NULL,
                                        h_offset = 0, use_envipat = TRUE, ppm = 10,
                                        rt_half_window = 1.5, extract = TRUE,
                                        python_bin = find_python(),
                                        module_dir = .find_deconv_module_dir()) {
  if (is.null(ms1_matches) || nrow(ms1_matches) == 0 || !"obs_mz" %in% names(ms1_matches)) return(list())
  if (is.null(study_samples)) study_samples <- unique(ms1_matches$sample)
  ev <- list()
  for (id in met_ids) {
    row <- .ev_best_row(ms1_matches, id, study_samples)
    if (is.null(row)) next
    env <- ms1_matches[ms1_matches$met_id == id & ms1_matches$sample == row$sample &
                         ms1_matches$adduct == row$adduct & ms1_matches$k_oxid == row$k_oxid, , drop = FALSE]
    env <- env[order(env$z), c("z", "theo_mz", "obs_mz", "ppm_error", "intensity"), drop = FALSE]
    iso <- .ev_nominal_isotopes(row$formula, row$z, row$theo_mz, h_offset, use_envipat)
    file <- NULL
    if (!is.null(features) && "source_file" %in% names(features)) {
      f <- unique(features$source_file[features$sample == row$sample])
      f <- f[!is.na(f) & file.exists(f)]
      if (length(f) > 0) file <- f[1]
    }
    ms2 <- NULL
    if (!is.null(ms2_confirmations) && nrow(ms2_confirmations) > 0 && !is.null(ms2_spectra)) {
      cand <- ms2_confirmations[ms2_confirmations$met_id == id, , drop = FALSE]
      if (nrow(cand) > 0) {
        cand <- cand[order(cand$sample != row$sample, -cand$confirmation_score), , drop = FALSE]
        for (j in seq_len(nrow(cand))) {
          key <- paste(cand$sample[j], id, cand$k_oxid[j], cand$z[j], cand$adduct[j], sep = "|")
          sp <- ms2_spectra[[key]]
          met <- if (!is.null(mets)) .find_met(mets, id) else NULL
          if (is.null(sp) || is.null(met)) next
          dat <- tryCatch(mirror_spectrum_data(met, sp, dict, h_offset = h_offset), error = function(e) NULL)
          if (!is.null(dat)) { ms2 <- list(data = dat, conf = cand[j, , drop = FALSE]); break }
        }
      }
    }
    ev[[id]] <- list(row = row, envelope = env, theo_iso = iso, file = file,
                     spectrum = NULL, xic = NULL, ms2 = ms2,
                     n_detected = length(unique(ms1_matches$sample[ms1_matches$met_id == id &
                                                                    ms1_matches$sample %in% study_samples])))
  }
  if (extract && length(ev) > 0) ev <- .ev_extract_raw(ev, ppm, rt_half_window, python_bin, module_dir)
  ev
}

# Theoretical isotope cluster at charge z, with fine-structure isotopologues
# summed into nominal peaks (as a mass spectrometer at typical resolution
# shows them), positioned so the monoisotopic peak sits at theo_mz (which
# carries any adduct and desulfurization shift of the matched ion).
.ev_nominal_isotopes <- function(formula, z, theo_mz, h_offset = 0, use_envipat = TRUE, n_iso = 10) {
  pat <- tryCatch(isotope_pattern(formula, threshold = 1e-6, n_top = 200, use_envipat = use_envipat),
                  error = function(e) NULL)
  if (is.null(pat) || nrow(pat) == 0) return(NULL)
  pat <- pat[order(pat$mass), ]
  k <- round(pat$mass - pat$mass[1])
  agg <- do.call(rbind, lapply(split(seq_len(nrow(pat)), k), function(ix)
    data.frame(iso = k[ix[1]], mass = sum(pat$mass[ix] * pat$abundance[ix]) / sum(pat$abundance[ix]),
               abundance = sum(pat$abundance[ix]))))
  agg <- agg[order(agg$iso), , drop = FALSE]
  agg <- head(agg[agg$abundance >= max(agg$abundance) * 0.01, , drop = FALSE], n_iso)
  agg$mz <- theo_mz + (agg$mass - pat$mass[1]) / z
  agg[, c("iso", "mz", "abundance", "mass")]
}

# One Python pass per raw file: apex MS1 spectrum around each precursor
# cluster, and the summed chromatogram of its monoisotopic ion plus its three
# most abundant isotopes.
.ev_extract_raw <- function(ev, ppm, rt_half_window, python_bin, module_dir) {
  req <- do.call(rbind, lapply(names(ev), function(id) {
    e <- ev[[id]]
    if (is.null(e$file) || is.na(e$row$rt)) return(NULL)
    z <- e$row$z
    iso <- e$theo_iso
    xmz <- if (!is.null(iso) && nrow(iso) > 0) unique(c(e$row$theo_mz, head(iso$mz[order(-iso$abundance)], 3))) else e$row$theo_mz
    data.frame(req_id = id, file = e$file, rt = e$row$rt,
               mz_lo = e$row$theo_mz - 2.5 / z, mz_hi = e$row$theo_mz + 10.5 / z,
               rt_lo = e$row$rt - rt_half_window, rt_hi = e$row$rt + rt_half_window,
               xic_mz = paste(sprintf("%.6f", xmz), collapse = ";"), ppm = ppm,
               stringsAsFactors = FALSE)
  }))
  if (is.null(req) || nrow(req) == 0) return(ev)
  out <- tryCatch({
    od <- tempfile("evidence_"); dir.create(od)
    rp <- file.path(od, "requests.tsv")
    utils::write.table(req, rp, sep = "\t", row.names = FALSE, quote = FALSE)
    .run_python_module(c("-m", "oligomet_deconv.evidence_cli", "--requests", rp, "--output-dir", od),
                       python_bin, module_dir, "Report evidence extraction")
    list(spec = utils::read.delim(file.path(od, "evidence_spectra.tsv"), stringsAsFactors = FALSE,
                                  colClasses = c(req_id = "character")),
         xic = utils::read.delim(file.path(od, "evidence_xic.tsv"), stringsAsFactors = FALSE,
                                 colClasses = c(req_id = "character")))
  }, error = function(e) { message("Evidence extraction skipped: ", conditionMessage(e)); NULL })
  if (is.null(out)) return(ev)
  for (id in names(ev)) {
    s <- out$spec[out$spec$req_id == id, , drop = FALSE]
    x <- out$xic[out$xic$req_id == id, , drop = FALSE]
    if (nrow(s) > 0) ev[[id]]$spectrum <- s[, c("scan_rt", "mz", "intensity")]
    if (nrow(x) > 0) ev[[id]]$xic <- x[order(x$rt), c("rt", "intensity")]
  }
  ev
}

## ---- Plots -------------------------------------------------------------------
.ev_theme <- function() {
  ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 9.5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 8, colour = "#44546a"))
}

#' Plot the extracted-ion chromatogram of a metabolite's precursor
#' @param e One evidence record (see [collect_metabolite_evidence()]).
#' @return A ggplot, or `NULL` without chromatogram data.
#' @export
plot_evidence_xic <- function(e) {
  if (is.null(e$xic) || nrow(e$xic) < 2) return(NULL)
  ggplot2::ggplot(e$xic, ggplot2::aes(x = .data$rt, y = .data$intensity)) +
    ggplot2::geom_area(fill = "#1f3864", alpha = 0.15) +
    ggplot2::geom_line(colour = "#1f3864", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = e$row$rt, linetype = "dashed", colour = "#c0392b", linewidth = 0.3) +
    ggplot2::labs(x = "Retention time (min)", y = "Intensity",
                  title = sprintf("Extracted-ion chromatogram, %s", .ev_ion_label(e$row$z, e$row$adduct)),
                  subtitle = sprintf("Sample %s; monoisotopic m/z %.4f plus the 3 most abundant isotopes, summed; apex %.2f min",
                                     e$row$sample, e$row$theo_mz, e$row$rt)) +
    .ev_theme()
}

#' Plot the acquired MS1 isotope pattern against the theoretical one
#' @param e One evidence record (see [collect_metabolite_evidence()]).
#' @param ppm Tolerance for labelling observed peaks as isotope matches.
#' @return A ggplot, or `NULL` without a theoretical pattern.
#' @export
plot_evidence_isotopes <- function(e, ppm = 10) {
  iso <- e$theo_iso
  if (is.null(iso) || nrow(iso) == 0) return(NULL)
  th <- data.frame(mz = iso$mz, rel = 100 * iso$abundance / max(iso$abundance))
  p <- ggplot2::ggplot()
  sub <- "Theoretical pattern only (raw spectrum not available)"
  if (!is.null(e$spectrum) && nrow(e$spectrum) > 0) {
    ob <- e$spectrum
    ob$rel <- 100 * ob$intensity / max(ob$intensity)
    near <- vapply(ob$mz, function(m) min(abs(th$mz - m)) / m * 1e6, numeric(1))
    ob$matched <- near <= ppm
    lab <- ob[ob$matched & ob$rel >= 5, , drop = FALSE]
    p <- p + ggplot2::geom_segment(data = ob, ggplot2::aes(x = .data$mz, xend = .data$mz, y = 0, yend = .data$rel,
                                                            colour = .data$matched), linewidth = 0.6) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "#1f3864", `FALSE` = "#a5b1c2"),
                                   labels = c(`TRUE` = "observed, isotope match", `FALSE` = "observed, other"),
                                   name = NULL) +
      ggplot2::geom_text(data = lab, ggplot2::aes(x = .data$mz, y = .data$rel, label = sprintf("%.4f", .data$mz)),
                         size = 2.2, vjust = -0.4, colour = "#1f3864")
    sub <- sprintf("Apex MS1 scan, %.2f min, sample %s; diamonds: theoretical isotope distribution",
                   ob$scan_rt[1], e$row$sample)
  }
  p + ggplot2::geom_point(data = th, ggplot2::aes(x = .data$mz, y = .data$rel), shape = 23,
                          fill = "#c0392b", colour = "#c0392b", size = 1.6) +
    ggplot2::scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
    ggplot2::labs(x = "m/z", y = "Relative intensity (%)",
                  title = sprintf("MS1 isotope pattern, %s  (theoretical m/z %.4f)",
                                  .ev_ion_label(e$row$z, e$row$adduct), e$row$theo_mz),
                  subtitle = sub) +
    .ev_theme() + ggplot2::theme(legend.position = "bottom")
}

#' Plot the charge-state envelope of a metabolite in its originating sample
#' @param e One evidence record (see [collect_metabolite_evidence()]).
#' @return A ggplot, or `NULL` with no envelope rows.
#' @export
plot_evidence_envelope <- function(e) {
  env <- e$envelope
  if (is.null(env) || nrow(env) == 0) return(NULL)
  env$zl <- factor(paste0(env$z, "-"), levels = paste0(sort(unique(env$z)), "-"))
  ggplot2::ggplot(env, ggplot2::aes(x = .data$zl, y = .data$intensity)) +
    ggplot2::geom_col(fill = "#1f3864", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("m/z %.4f\n%+.1f ppm", .data$obs_mz, .data$ppm_error)),
                       vjust = -0.25, size = 2.2, lineheight = 0.9) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.3))) +
    ggplot2::labs(x = "Charge state", y = "Intensity",
                  title = "Charge-state envelope",
                  subtitle = sprintf("Sample %s; %s", e$row$sample,
                                     if (identical(e$row$adduct, "H")) "deprotonated ions" else paste0(e$row$adduct, " adduct"))) +
    .ev_theme()
}
