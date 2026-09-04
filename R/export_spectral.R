# =============================================================================
# export_spectral.R
# Spectral-library export: theoretical MS1 and MS2 libraries as MGF and MSP.
#
# Both formats are plain text and are read by the common spectral-library
# tools (MS-DIAL, mzVault/Compound Discoverer, MZmine, matchms, SIRIUS,
# GNPS). Each "spectrum" here is theoretical, not measured:
#
#   MS1 library -- one spectrum per (metabolite, PS-oxidation level, charge
#     state). Peaks are the theoretical isotope cluster at that charge, so
#     the peak intensities are real relative abundances from the isotope
#     calculation and can be matched against an acquired isotope pattern.
#
#   MS2 library -- one spectrum per (metabolite, precursor charge state).
#     Peaks are the McLuckey fragment ions (terminal a/a-B/b/b-B/c/w/x/y and,
#     optionally, w-a/w-b/w-d internal ions) at the requested fragment
#     charges. Intensities come from a rule-based heuristic
#     (fragment_intensity_weight() in fragments.R) that encodes known
#     oligo fragmentation propensities -- PS-linkage lability, MOE-PS vs
#     DNA-PS dominant ion series, purine vs pyrimidine base-loss lability --
#     it is not fit to any measured spectra, so treat it as a coarse
#     relative ranking, not a predicted abundance. Use these libraries for
#     m/z matching and annotation; do not run quantitative dot-product /
#     cosine scoring against them as if the intensities were measured.
#
# Format notes:
#   MGF  -- BEGIN IONS / END IONS blocks, TITLE / PEPMASS / CHARGE headers.
#           Peak lines are "m/z intensity" only; per-peak annotations do not
#           have a portable place in MGF, so they are carried in the TITLE
#           and in the MSP twin of the same library.
#   MSP  -- NIST-style "Name:" ... "Num Peaks:" blocks. Peak lines carry the
#           annotation as a quoted third column, which MS-DIAL and matchms
#           both read.
#
# Negative ESI throughout ([M - zH]^z-), matching the rest of the pipeline.
# =============================================================================

## ---- Spectrum record --------------------------------------------------------
# Internal representation shared by both writers.
#   name          spectrum name (unique within the library)
#   precursor_mz  precursor m/z (negative mode)
#   z             charge state as a positive integer (polarity is implied)
#   formula       neutral molecular formula string
#   mono_mass     neutral monoisotopic mass
#   level         "MS1" or "MS2"
#   fields        named character vector of extra metadata (written to the
#                 MGF TITLE and as MSP header lines / Comment)
#   peaks         data.frame(mz, intensity, annotation)
.spectrum_record <- function(name, precursor_mz, z, formula, mono_mass,
                             level, fields, peaks) {
  list(name = name, precursor_mz = precursor_mz, z = z, formula = formula,
       mono_mass = mono_mass, level = level, fields = fields, peaks = peaks)
}

# "[M-5H]5-" for z = 5.
.precursor_type <- function(z) {
  if (z == 1) "[M-H]-" else sprintf("[M-%dH]%d-", z, z)
}

# Flatten a named vector to "k=v; k=v" for TITLE / Comment lines.
.fields_to_string <- function(fields) {
  if (length(fields) == 0) return("")
  paste(sprintf("%s=%s", names(fields), as.character(fields)), collapse = "; ")
}

# Truncate a record list to a cap, warning once with the reason.
.cap_records <- function(recs, max_spectra, what, knobs) {
  if (length(recs) <= max_spectra) return(recs)
  warning(what, " has ", length(recs), " spectra; truncating to ",
          max_spectra, ". Narrow ", knobs, ", or raise max_spectra, to keep ",
          "more.", call. = FALSE)
  recs[seq_len(max_spectra)]
}

## ---- MS1 library ------------------------------------------------------------
# One spectrum per (metabolite, oxidation level k, charge state z); peaks are
# the theoretical isotope cluster. Intensities are relative isotope
# abundances rescaled so the largest peak in each spectrum is 100.
build_ms1_library <- function(mets, dict = STANDARD_DICT, z_range = 3:12,
                              n_iso = 8, max_oxid = 0, h_offset = 0,
                              use_envipat = TRUE, oligo_name = NULL,
                              max_spectra = 5000) {
  recs <- list()
  for (met in mets) {
    env <- compute_envelope(met, z_range = z_range, n_iso = n_iso,
                            max_oxid = max_oxid, h_offset = h_offset,
                            use_envipat = use_envipat, dict = dict)
    if (is.null(env) || nrow(env) == 0) next
    grp <- split(env, list(env$k_oxid, env$z), drop = TRUE)
    for (g in grp) {
      g <- g[order(g$mz), ]
      k <- g$k_oxid[1]
      z <- g$z[1]
      nm <- paste0(met$name, if (k > 0) paste0(" +", k, "Ox") else "",
                   " [z=", z, "]")
      inten <- g$abundance
      if (max(inten) > 0) inten <- 100 * inten / max(inten)
      fields <- c(
        met_id     = met$id,
        met_kind   = met$kind,
        length     = as.character(met$n),
        k_oxid     = as.character(k),
        n_ps       = as.character(met$n_ps),
        if (!is.null(oligo_name)) c(oligo = oligo_name) else character(0)
      )
      # compute_envelope()'s `iso` column is the index of the peak within the
      # selected cluster, not the isotopologue's nominal offset -- the
      # selection keeps the monoisotopic peak plus the top-(n_iso-1) by
      # abundance, so for a 7 kDa oligo the second peak kept is typically
      # M+7, not M+1. Recover the real offset from the neutral mass.
      neutral <- g$mz * z + z * .PROTON - h_offset
      offset <- round(neutral - g$mono_mass)

      # The isotope engine bins masses to 0.0001 Da at every convolution
      # step, and for a 7 kDa formula that rounding compounds to ~1 mDa by
      # the end -- small (0.2 ppm) but enough to put the library's peaks off
      # the exact masses the rest of the pipeline reports. Peak *spacing* is
      # unaffected by the accumulated offset, so anchor the cluster by
      # shifting it rigidly onto the m/z computed directly from the formula.
      prec_mz <- (g$mono_mass[1] + h_offset - z * .PROTON) / z
      mz <- g$mz
      if (any(offset == 0)) mz <- mz + (prec_mz - mz[offset == 0][1])

      recs[[length(recs) + 1]] <- .spectrum_record(
        name = nm, precursor_mz = prec_mz, z = z,
        formula = g$formula[1], mono_mass = g$mono_mass[1],
        level = "MS1", fields = fields,
        peaks = data.frame(mz = mz, intensity = inten,
                           annotation = paste0("M+", offset),
                           stringsAsFactors = FALSE))
    }
  }
  .cap_records(recs, max_spectra, "MS1 spectral library",
               "z_range/max_oxid/metabolite selection")
}

## ---- MS2 library ------------------------------------------------------------
# Compact ion label: "w4^2-", "a-B5^1-", "w-a(3,8)^1-".
.fragment_label <- function(f, z, include_formula = FALSE) {
  base <- switch(f$ion_type,
                 aB = "a-B", bB = "b-B", f$ion_type)
  site <- if (!is.na(f$internal_5) && !is.na(f$internal_3)) {
    sprintf("(%d,%d)", f$internal_5, f$internal_3)
  } else {
    as.character(f$frag_length)
  }
  lbl <- sprintf("%s%s^%d-", base, site, z)
  if (include_formula) lbl <- paste0(lbl, " [", f$formula, "]")
  lbl
}

# Drop low-relative-intensity peaks from an already mz-sorted, 0-100-rescaled
# peak table before export -- a long tail of near-zero-weight ions (internal
# fragments, weak base losses) clutters exported libraries with entries no
# acquisition method would realistically detect at typical noise floors.
# Note: under the current fragment_intensity_weight() constants (fragments.R:
# PS_LINKAGE_BOOST=1.5, CHEMISTRY_ION_BOOST=1.6, PURINE_LOSS_BOOST=1.3,
# INTERNAL_ION_WEIGHT=0.35), the weakest possible fragment is never below
# ~11% of the strongest (0.35 / (1.6*1.5*1.3)), so the default 5% threshold
# is a forward-looking safety net -- it won't visibly thin today's libraries,
# but guards against a future weight-model change (or a custom-chemistry
# combination) producing a wider intensity spread than exists today.
.filter_min_rel_intensity <- function(peaks, min_rel_intensity = 0.05) {
  if (min_rel_intensity <= 0 || nrow(peaks) == 0) return(peaks)
  peaks[peaks$intensity >= min_rel_intensity * max(peaks$intensity), , drop = FALSE]
}

# One spectrum per (metabolite, precursor charge state). Peaks are the
# theoretical fragment ions at the requested fragment charges; intensities
# are a rule-based relative-abundance heuristic (fragment_intensity_weight()
# in fragments.R), not measured or calibrated data -- see that function's
# header for exactly what it does and does not encode.
build_ms2_library <- function(mets, dict = STANDARD_DICT,
                              precursor_z_range = 4:7, frag_z_range = 1:2,
                              ion_types = c("a", "aB", "b", "bB", "w", "y"),
                              include_internal = FALSE, h_offset = 0,
                              mz_min = 100, mz_max = 6000,
                              min_rel_intensity = 0.05,
                              oligo_name = NULL, max_spectra = 2000) {
  recs <- list()
  for (met in mets) {
    if (met$n < 3) next
    info <- metabolite_mass_info(met, dict)
    frags <- generate_fragments(met, dict, ion_types = ion_types,
                                z_range = frag_z_range)
    if (include_internal) {
      frags <- c(frags, generate_internal_fragments(met, dict,
                                                    z_range = frag_z_range))
    }
    if (length(frags) == 0) next

    # Flatten fragments to a peak table once; it does not depend on the
    # precursor charge state. Each fragment's heuristic weight is the same
    # across its own charge states (the heuristic doesn't model fragment
    # charge), so it's computed once per fragment, not per peak.
    mz <- numeric(0)
    annot <- character(0)
    weight <- numeric(0)
    for (f in frags) {
      fw <- fragment_intensity_weight(f, met, dict)
      for (i in seq_len(nrow(f$mz_table))) {
        m <- f$mz_table$mz[i]
        if (m < mz_min || m > mz_max) next
        mz <- c(mz, m)
        annot <- c(annot, .fragment_label(f, f$mz_table$z[i], include_formula = TRUE))
        weight <- c(weight, fw)
      }
    }
    if (length(mz) == 0) next
    ord <- order(mz)
    inten <- weight[ord]
    if (max(inten) > 0) inten <- 100 * inten / max(inten)
    peaks <- data.frame(mz = mz[ord], intensity = inten,
                        annotation = annot[ord], stringsAsFactors = FALSE)
    peaks <- .filter_min_rel_intensity(peaks, min_rel_intensity)
    if (nrow(peaks) == 0) next

    for (z in precursor_z_range) {
      prec_mz <- (info$mono_mass + h_offset - z * .PROTON) / z
      fields <- c(
        met_id      = met$id,
        met_kind    = met$kind,
        length      = as.character(met$n),
        n_fragments = as.character(nrow(peaks)),
        if (!is.null(oligo_name)) c(oligo = oligo_name) else character(0)
      )
      recs[[length(recs) + 1]] <- .spectrum_record(
        name = paste0(met$name, " [z=", z, "] MS2"),
        precursor_mz = prec_mz, z = z, formula = info$formula_str,
        mono_mass = info$mono_mass, level = "MS2", fields = fields,
        peaks = peaks)
    }
  }
  .cap_records(recs, max_spectra, "MS2 spectral library",
               "precursor_z_range/metabolite selection")
}

## ---- MGF writer -------------------------------------------------------------
write_mgf <- function(records, file) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "# Theoretical spectral library generated by OligoMetProfiler",
    "# Negative ESI; precursor type [M-zH]z-.",
    "# MS2 peak intensities are a rule-based heuristic, not measured data.",
    .disclaimer_comment_block(),
    ""), con)
  for (r in records) {
    title <- paste0(r$name,
                    " | level=", r$level,
                    " | formula=", r$formula,
                    " | mono_mass=", sprintf("%.4f", r$mono_mass),
                    if (nzchar(.fields_to_string(r$fields)))
                      paste0(" | ", .fields_to_string(r$fields)) else "")
    writeLines(c(
      "BEGIN IONS",
      paste0("TITLE=", title),
      paste0("PEPMASS=", sprintf("%.4f", r$precursor_mz)),
      paste0("CHARGE=", r$z, "-"),
      sprintf("%.4f %.2f", r$peaks$mz, r$peaks$intensity),
      "END IONS",
      ""), con)
  }
  invisible(normalizePath(file))
}

## ---- MSP writer -------------------------------------------------------------
# measured = TRUE labels records as acquired/annotated data (peaks and
# precursor mz/intensity are real instrument values; only the per-peak
# fragment-ion ANNOTATION is a prediction) rather than a fully theoretical
# spectrum -- see batch_annotated_msp_records() in R/mirror_plot.R for the
# one caller that uses this.
write_msp <- function(records, file, measured = FALSE) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  for (r in records) {
    # MSP has no portable comment-line syntax, so the footer goes in each
    # record's Comment field rather than in a file header -- a stray "#" block
    # ahead of the first Name: is not something every MSP reader tolerates.
    # Kept to the footer line; the full statement is in DISCLAIMER.md.
    comment <- if (measured) {
      paste0("Acquired ", r$level, " spectrum -- ", OLIGOMET_FOOTER,
             " Peaks are measured instrument data; per-peak fragment-ion ",
             "annotations (where present) are a computed match, not a ",
             "measurement",
             if (nzchar(.fields_to_string(r$fields)))
               paste0("; ", .fields_to_string(r$fields)) else "")
    } else {
      paste0("Theoretical ", r$level, " spectrum -- ",
                      OLIGOMET_FOOTER,
                      " Computed prediction, not a measurement",
                      if (r$level == "MS2")
                        "; intensities are a rule-based heuristic, not measured" else "",
                      if (nzchar(.fields_to_string(r$fields)))
                        paste0("; ", .fields_to_string(r$fields)) else "")
    }
    writeLines(c(
      paste0("NAME: ", r$name),
      paste0("PRECURSORMZ: ", sprintf("%.4f", r$precursor_mz)),
      paste0("PRECURSORTYPE: ", .precursor_type(r$z)),
      paste0("FORMULA: ", r$formula),
      paste0("EXACTMASS: ", sprintf("%.4f", r$mono_mass)),
      "IONMODE: Negative",
      # No COLLISIONENERGY line: build_ms2_library()'s peak intensities are
      # not NCE-dependent (fragment_intensity_weight() doesn't model NCE),
      # so there is no principled per-record NCE value to report yet -- see
      # this file's header / DEFAULT_PIPELINE_PARAMS$hcd_nce, which is an
      # acquisition-method starting-value suggestion, not a spectrum property.
      paste0("SPECTRUMTYPE: ", if (r$level == "MS1") "MS1" else "MS2"),
      paste0("COMMENT: ", comment),
      paste0("Num Peaks: ", nrow(r$peaks)),
      sprintf("%.4f\t%.2f\t\"%s\"", r$peaks$mz, r$peaks$intensity,
              r$peaks$annotation),
      ""), con)
  }
  invisible(normalizePath(file))
}

## ---- Convenience: write all four files --------------------------------------
# Writes <prefix>_MS1_library.{mgf,msp} and <prefix>_MS2_library.{mgf,msp}
# into out_dir. Returns the four paths, named ms1_mgf/ms1_msp/ms2_mgf/ms2_msp.
export_spectral_libraries <- function(mets, dict = STANDARD_DICT,
                                      out_dir = ".", prefix = "oligo",
                                      z_range = 3:12, n_iso = 8, max_oxid = 0,
                                      precursor_z_range = 4:7,
                                      frag_z_range = 1:2,
                                      include_internal = FALSE, h_offset = 0,
                                      use_envipat = TRUE, oligo_name = NULL,
                                      max_ms1_spectra = 5000,
                                      max_ms2_spectra = 2000) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ms1 <- build_ms1_library(mets, dict, z_range = z_range, n_iso = n_iso,
                           max_oxid = max_oxid, h_offset = h_offset,
                           use_envipat = use_envipat, oligo_name = oligo_name,
                           max_spectra = max_ms1_spectra)
  ms2 <- build_ms2_library(mets, dict, precursor_z_range = precursor_z_range,
                           frag_z_range = frag_z_range,
                           include_internal = include_internal,
                           h_offset = h_offset, oligo_name = oligo_name,
                           max_spectra = max_ms2_spectra)
  paths <- c(
    ms1_mgf = file.path(out_dir, paste0(prefix, "_MS1_library.mgf")),
    ms1_msp = file.path(out_dir, paste0(prefix, "_MS1_library.msp")),
    ms2_mgf = file.path(out_dir, paste0(prefix, "_MS2_library.mgf")),
    ms2_msp = file.path(out_dir, paste0(prefix, "_MS2_library.msp"))
  )
  write_mgf(ms1, paths[["ms1_mgf"]])
  write_msp(ms1, paths[["ms1_msp"]])
  write_mgf(ms2, paths[["ms2_mgf"]])
  write_msp(ms2, paths[["ms2_msp"]])
  paths
}
