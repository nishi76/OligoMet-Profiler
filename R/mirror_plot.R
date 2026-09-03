# =============================================================================
# mirror_plot.R
# Mirror plot for MS2 confirmation: acquired spectrum vs. theoretical
# fragment library, in the style of the acquired-vs-library mirror plots
# used to report spectral matches (e.g. Fig. 2, Chambers et al. 2020, Nat
# Commun, https://doi.org/10.1038/s41467-020-14665-7) -- acquired spectrum
# on top (positive intensity), theoretical library spectrum on the bottom
# (negative intensity), matched peaks colored, unmatched peaks gray,
# fragment-ion labels on matched peaks only.
# =============================================================================

## ---- Helper: locate a metabolite by id -------------------------------------
.find_met <- function(mets, met_id) {
  idx <- which(vapply(mets, function(m) identical(m$id, met_id), logical(1)))
  if (length(idx) == 0) return(NULL)
  mets[[idx[1]]]
}

## ---- Build mirror-plot data -------------------------------------------------
# Confirms `met` against one acquired MS2 spectrum (`ms2_peaks`: data.frame
# with mz, intensity columns) and returns both sides of the mirror plot,
# each peak flagged matched/unmatched and (for matched peaks) labeled with
# its fragment-ion annotation (e.g. "w5^1-", matching the MS2 library's
# convention -- see .fragment_label() in export_spectral.R).
#
# Both sides are independently rescaled to a 0-100 relative-intensity scale
# so they read as comparable shapes -- acquired intensities are real
# instrument counts, theoretical ones are the rule-based heuristic from
# fragment_intensity_weight() (see that function's header for exactly what
# it does and does not encode), so only the *relative* peak pattern is
# meaningful, never a shared absolute scale.
mirror_spectrum_data <- function(met, ms2_peaks, dict = STANDARD_DICT,
                                  tol_ppm = 25, z_range = 1:2,
                                  ion_types = c("a", "aB", "b", "bB", "w", "y"),
                                  include_internal = FALSE, h_offset = 0) {
  conf <- confirm_metabolite(met, ms2_peaks, dict, tol_ppm = tol_ppm,
                              z_range = z_range, ion_types = ion_types,
                              include_internal = include_internal,
                              h_offset = h_offset)
  frags <- conf$fragments
  matched <- conf$matched

  # Theoretical side: one row per (fragment, charge). `key` disambiguates
  # fragments that happen to share an m/z, so matching against `matched`
  # below never relies on floating-point mz equality.
  theo_mz <- numeric(0); theo_w <- numeric(0)
  theo_label <- character(0); theo_key <- character(0)
  for (f in frags) {
    fw <- fragment_intensity_weight(f, met, dict)
    for (i in seq_len(nrow(f$mz_table))) {
      z <- f$mz_table$z[i]
      theo_mz    <- c(theo_mz, f$mz_table$mz[i])
      theo_w     <- c(theo_w, fw)
      theo_label <- c(theo_label, .fragment_label(f, z))
      theo_key   <- c(theo_key,
                       paste(f$ion_type, f$direction, f$cleavage_site, f$frag_length, z))
    }
  }
  if (length(theo_w) > 0 && max(theo_w) > 0) theo_w <- 100 * theo_w / max(theo_w)

  matched_key <- if (nrow(matched) > 0) {
    paste(matched$ion_type, matched$direction, matched$cleavage_site,
          matched$frag_length, matched$z)
  } else character(0)
  is_matched_theo <- theo_key %in% matched_key

  theoretical <- data.frame(
    mz = theo_mz, intensity = theo_w, matched = is_matched_theo,
    label = ifelse(is_matched_theo, theo_label, NA_character_),
    stringsAsFactors = FALSE
  )

  # Acquired side: every observed peak; labeled (and colored) where it was
  # claimed by a matched theoretical fragment. matched$obs_mz is copied
  # verbatim from ms2_peaks$mz by match_fragments() (not recomputed), so
  # exact equality reliably finds the originating row.
  acquired <- ms2_peaks[, c("mz", "intensity")]
  max_int <- suppressWarnings(max(acquired$intensity, na.rm = TRUE))
  if (is.finite(max_int) && max_int > 0) acquired$intensity <- 100 * acquired$intensity / max_int
  acquired$matched <- FALSE
  acquired$label <- NA_character_
  if (nrow(matched) > 0) {
    for (i in seq_len(nrow(matched))) {
      idx <- which(ms2_peaks$mz == matched$obs_mz[i])
      if (length(idx) == 0) next
      idx <- idx[1]
      lbl <- theo_label[theo_key == matched_key[i]][1]
      acquired$matched[idx] <- TRUE
      acquired$label[idx] <- if (is.na(acquired$label[idx])) lbl
                              else paste(acquired$label[idx], lbl, sep = " / ")
    }
  }

  list(acquired = acquired, theoretical = theoretical,
       score = conf$score, met_name = met$name)
}

## ---- Mirror plot -------------------------------------------------------------
# `spec` is the output of mirror_spectrum_data(). Acquired peaks are drawn
# upward (positive intensity), theoretical library peaks downward (negative
# intensity); matched peaks are colored, unmatched peaks gray; only matched
# peaks are labeled, since an unfiltered theoretical library can carry
# hundreds of unmatched candidate fragments.
plot_mirror_spectrum <- function(spec, title = NULL,
                                  matched_color = "#0279EE",
                                  unmatched_color = "#B7B2A7") {
  acq <- spec$acquired
  theo <- spec$theoretical
  acq$side <- "Acquired"
  theo$side <- "Theoretical"
  theo$intensity <- -theo$intensity

  d <- rbind(acq[, c("mz", "intensity", "matched", "label", "side")],
             theo[, c("mz", "intensity", "matched", "label", "side")])
  d$label_y <- d$intensity + ifelse(d$intensity >= 0, 4, -4)
  d$label_hjust <- ifelse(d$intensity >= 0, 0, 1)
  labels_df <- d[!is.na(d$label), , drop = FALSE]

  subtitle <- if (!is.null(spec$score)) {
    sprintf("Acquired (top) vs. theoretical library (bottom) -- %d matched ions, %.0f%% sequence coverage",
            spec$score$n_matches, 100 * spec$score$coverage)
  } else {
    "Acquired (top) vs. theoretical library (bottom)"
  }

  ggplot2::ggplot(d, ggplot2::aes(x = mz, y = intensity, color = matched)) +
    ggplot2::geom_hline(yintercept = 0, color = "#333333", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0), linewidth = 0.4) +
    ggplot2::geom_text(data = labels_df,
                        ggplot2::aes(y = label_y, label = label, hjust = label_hjust),
                        angle = 90, size = 2.6, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c(`TRUE` = matched_color, `FALSE` = unmatched_color),
                                 labels = c(`TRUE` = "Matched", `FALSE` = "Unmatched"),
                                 name = NULL) +
    ggplot2::scale_y_continuous(labels = function(x) abs(x)) +
    ggplot2::expand_limits(y = c(-135, 135)) +
    ggplot2::labs(x = "m/z", y = "Relative intensity (%)",
                  title = title, subtitle = subtitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
}

## ---- Batch MS2 confirmations: shared per-hit resolution ---------------------
# One row per confirmed batch hit, each carrying its acquired spectrum, the
# metabolite it was confirmed against, and a computed theoretical precursor
# m/z (ms2_confirmations itself doesn't carry precursor_mz -- back-computed
# here the same way match_ms1() does: oxidized mass + adduct shift, charge-
# reduced). Shared by both the mirror-plot PDF and the annotated MSP export
# below so they stay in lockstep with exactly the same hit list.
#
# ms2_confirmations: data.frame from confirm_ms2_batch()/annotate_metabolites_batch().
# ms2_spectra:        named list from the same call, keyed "sample|met_id|k_oxid|z|adduct".
.resolve_batch_ms2_hits <- function(mets, dict, ms2_confirmations, ms2_spectra,
                                     h_offset = 0) {
  if (is.null(ms2_confirmations) || nrow(ms2_confirmations) == 0) return(list())
  out <- list()
  for (i in seq_len(nrow(ms2_confirmations))) {
    row <- ms2_confirmations[i, ]
    key <- paste(row$sample, row$met_id, row$k_oxid, row$z, row$adduct, sep = "|")
    best_spec <- ms2_spectra[[key]]
    met <- .find_met(mets, row$met_id)
    if (is.null(best_spec) || is.null(met)) next

    info <- metabolite_mass_info(met, dict)
    mass <- ps_oxid_mass(info$mono_mass, row$k_oxid)
    ad_shift <- if (identical(row$adduct, "H")) 0 else adduct_shift(row$adduct)
    theo_mass <- mass + ad_shift
    precursor_mz <- (theo_mass + h_offset - row$z * .PROTON) / row$z

    out[[key]] <- list(row = row, met = met, best_spec = best_spec,
                        info = info, precursor_mz = precursor_mz)
  }
  out
}

## ---- Downloadable mirror plots (multi-page PDF) ------------------------------
# One page per confirmed batch MS2 hit -- reuses mirror_spectrum_data()/
# plot_mirror_spectrum() exactly as the interactive single-hit plot does, so
# the PDF matches what's shown in the app.
batch_mirror_plots_pdf <- function(mets, dict, ms2_confirmations, ms2_spectra,
                                    file, tol_ppm = 25, z_range = 1:2,
                                    h_offset = 0, width = 9, height = 5) {
  hits <- .resolve_batch_ms2_hits(mets, dict, ms2_confirmations, ms2_spectra, h_offset)
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (length(hits) == 0) {
    plot.new()
    text(0.5, 0.5, "No confirmed batch MS2 hits to plot.")
    return(invisible(file))
  }
  for (h in hits) {
    spec <- mirror_spectrum_data(h$met, h$best_spec, dict, tol_ppm = tol_ppm,
                                  z_range = z_range, h_offset = h_offset)
    p <- plot_mirror_spectrum(spec, title = paste0(
      h$met$name, " (", h$row$sample, ") -- z=", h$row$z, " -- MS2 mirror plot"))
    print(p)
  }
  invisible(file)
}

## ---- Downloadable annotated MSP (acquired batch MS2 spectra) ----------------
# Exports the ACQUIRED spectrum for every confirmed batch hit as an MSP
# record -- peaks are real instrument data, and each peak's annotation is
# the fragment ion the app's own in-memory confirmation matched it to
# (blank where unmatched). Unlike re-matching against an exported MSP
# library file (precursor-mz + charge lookup against the closest library
# entry), this reuses the exact per-hit confirmation already computed for
# that (sample, metabolite, charge, adduct) -- more precise, since it knows
# exactly which fragment matched rather than the nearest library candidate.
batch_annotated_msp_records <- function(mets, dict, ms2_confirmations, ms2_spectra,
                                         tol_ppm = 25, z_range = 1:2, h_offset = 0) {
  hits <- .resolve_batch_ms2_hits(mets, dict, ms2_confirmations, ms2_spectra, h_offset)
  recs <- list()
  for (h in hits) {
    spec <- mirror_spectrum_data(h$met, h$best_spec, dict, tol_ppm = tol_ppm,
                                  z_range = z_range, h_offset = h_offset)
    acq <- spec$acquired
    acq <- acq[order(acq$mz), ]
    recs[[length(recs) + 1]] <- .spectrum_record(
      name = paste0(h$met$name, " (", h$row$sample, ") [z=", h$row$z, "] MS2 -- acquired"),
      precursor_mz = h$precursor_mz, z = h$row$z,
      formula = h$info$formula_str, mono_mass = h$info$mono_mass,
      level = "MS2",
      fields = c(sample = h$row$sample, met_id = h$row$met_id,
                 confirmation_score = as.character(h$row$confirmation_score),
                 coverage = as.character(round(h$row$coverage, 3))),
      peaks = data.frame(mz = acq$mz, intensity = acq$intensity,
                         annotation = ifelse(is.na(acq$label), "", acq$label),
                         stringsAsFactors = FALSE))
  }
  recs
}
