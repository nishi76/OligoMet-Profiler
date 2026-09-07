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
.fragment_label <- function(f, z) {
  base <- switch(f$ion_type,
                 aB = "a-B", bB = "b-B", f$ion_type)
  site <- if (!is.na(f$internal_5) && !is.na(f$internal_3)) {
    sprintf("(%d,%d)", f$internal_5, f$internal_3)
  } else {
    as.character(f$frag_length)
  }
  sprintf("%s%s^%d-", base, site, z)
}

# One spectrum per (metabolite, precursor charge state). Peaks are the
# theoretical fragment ions at the requested fragment charges; intensities
# are a rule-based relative-abundance heuristic (fragment_intensity_weight()
# in fragments.R), not measured or calibrated data -- see that function's
# header for exactly what it does and does not encode.
build_ms2_library <- function(mets, dict = STANDARD_DICT,
                              precursor_z_range = 4:7, frag_z_range = 1:2,
                              ion_types = c("aB", "w", "y", "b"),
                              include_internal = FALSE, h_offset = 0,
                              mz_min = 100, mz_max = 6000,
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
        annot <- c(annot, .fragment_label(f, f$mz_table$z[i]))
        weight <- c(weight, fw)
      }
    }
    if (length(mz) == 0) next
    ord <- order(mz)
    inten <- weight[ord]
    if (max(inten) > 0) inten <- 100 * inten / max(inten)
    peaks <- data.frame(mz = mz[ord], intensity = inten,
                        annotation = annot[ord], stringsAsFactors = FALSE)

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

## ---- Empirical (measured) MS2 library from real DDA data --------------------
# build_ms2_library() above is a rule-based heuristic: it encodes known
# fragmentation propensities, but its peak intensities are not fit to any
# measured spectrum. This builds the real alternative -- one consensus MS2
# spectrum per (metabolite, PS-oxidation level, charge, adduct) from the
# ACTUAL acquired fragment spectra already confirmed by the batch MS2
# pipeline (confirm_ms2_batch()/annotate_metabolites_batch() in
# R/batch_ms_processing.R), pooled across every sample/replicate that
# confirmed that species.
#
# Why pooling helps: any single acquired MS2 spectrum carries chemical/
# electronic noise particular to that one scan. Requiring a peak to recur,
# within fragment_ppm, across a large-enough fraction of the INDEPENDENT
# spectra pooled for the same species is what turns a pile of individually
# noisy real spectra into one clean, repeatable fragment-ion set -- the
# same logic a NIST-style consensus spectrum uses. A species confirmed in
# only one sample/replicate has nothing to vote against yet, so all of its
# peaks are kept as-is.
#
# ms2_confirmations / ms2_spectra: the two objects returned together by
# confirm_ms2_batch() (also inside annotate_metabolites_batch()'s result) --
# ms2_confirmations is one row per (sample, met_id, k_oxid, z, adduct) hit
# that had a real MS2 spectrum nearby, and ms2_spectra is that spectrum
# itself (data.frame(mz, intensity)), keyed "sample|met_id|k_oxid|z|adduct".
# Pooling here is NOT gated on confirmation_score/coverage -- doing so would
# bias the empirical library toward spectra that already agree with the
# very heuristic model this function exists to replace. Any real spectrum
# whose precursor was MS1-confirmed contributes; peak recurrence across
# replicates is the only denoising step.
#
# label_peaks optionally annotates each consensus peak with the matching
# theoretical fragment ion (a-B9^2- etc., generate_fragments()/
# match_fragments() in R/fragments.R) purely for readability -- exactly
# like build_ms2_library()'s labeling, this never removes or adds a peak,
# it only names one that's already there when a theoretical fragment lines
# up within fragment_ppm.
#
# Returns list(records = <spectrum records for write_msp()/write_mgf()>,
#              summary = one row per group, whether or not it made the cut).
build_empirical_ms2_library <- function(mets, dict = STANDARD_DICT,
                                         ms2_confirmations, ms2_spectra,
                                         ion_types = c("aB", "w", "y", "b"),
                                         frag_z_range = 1:2, h_offset = 0,
                                         fragment_ppm = 10,
                                         min_consensus_fraction = 0.5,
                                         min_consensus_peak_intensity = 0.01,
                                         min_peaks_for_library = 2,
                                         label_peaks = TRUE,
                                         max_spectra = 2000) {
  empty <- list(records = list(), summary = data.frame())
  if (is.null(ms2_confirmations) || nrow(ms2_confirmations) == 0) return(empty)

  grp_key <- paste(ms2_confirmations$met_id, ms2_confirmations$k_oxid,
                   ms2_confirmations$z, ms2_confirmations$adduct, sep = "|")
  summary_rows <- list()
  recs <- list()

  for (key in unique(grp_key)) {
    sub <- ms2_confirmations[grp_key == key, , drop = FALSE]
    met <- mets[[which(vapply(mets, function(m) identical(m$id, sub$met_id[1]), logical(1)))]]

    # Keep each spectrum paired with its originating sample through the
    # filter below -- contributors can drop rows at any position (a
    # confirmed hit whose spectrum key doesn't resolve), so slicing
    # sub$sample by position after filtering would misattribute samples.
    raw <- lapply(seq_len(nrow(sub)), function(i) {
      spec_key <- paste(sub$sample[i], sub$met_id[i], sub$k_oxid[i], sub$z[i],
                        sub$adduct[i], sep = "|")
      list(sample = sub$sample[i], spectrum = ms2_spectra[[spec_key]])
    })
    raw <- Filter(function(x) !is.null(x$spectrum) && nrow(x$spectrum) > 0, raw)
    contributors <- lapply(raw, function(x) x$spectrum)
    contrib_samples <- vapply(raw, function(x) x$sample, character(1))

    row <- data.frame(
      met_id = met$id, met_name = met$name, kind = met$kind, n = met$n,
      k_oxid = sub$k_oxid[1], z = sub$z[1], adduct = sub$adduct[1],
      n_source_spectra = length(contributors),
      source_samples = paste(sort(unique(contrib_samples)), collapse = ";"),
      stringsAsFactors = FALSE
    )
    if (length(contributors) == 0) {
      row$n_consensus_peaks <- 0L; row$written_to_library <- FALSE
      summary_rows[[length(summary_rows) + 1]] <- row
      next
    }

    consensus <- .consensus_ms2_peaks(contributors, fragment_ppm,
                                      min_consensus_fraction,
                                      min_consensus_peak_intensity)
    row$n_consensus_peaks <- nrow(consensus)
    if (nrow(consensus) < min_peaks_for_library) {
      row$written_to_library <- FALSE
      summary_rows[[length(summary_rows) + 1]] <- row
      next
    }
    row$written_to_library <- TRUE
    summary_rows[[length(summary_rows) + 1]] <- row

    consensus$annotation <- NA_character_
    if (label_peaks) {
      frags <- generate_fragments(met, dict, ion_types = ion_types, z_range = frag_z_range)
      hits <- match_fragments(frags, consensus[, c("mz", "intensity")],
                              tol_ppm = fragment_ppm, z_range = frag_z_range)
      if (nrow(hits) > 0) {
        for (h in seq_len(nrow(hits))) {
          idx <- which(consensus$mz == hits$obs_mz[h])
          if (length(idx) == 0) next
          base <- switch(hits$ion_type[h], aB = "a-B", bB = "b-B", hits$ion_type[h])
          lbl <- sprintf("%s%d^%d-", base, hits$frag_length[h], hits$z[h])
          consensus$annotation[idx[1]] <- if (is.na(consensus$annotation[idx[1]])) lbl
                                          else paste(consensus$annotation[idx[1]], lbl, sep = " / ")
        }
      }
    }
    consensus$annotation <- ifelse(is.na(consensus$annotation), "", consensus$annotation)

    info <- metabolite_mass_info(met, dict)
    mass <- ps_oxid_mass(info$mono_mass, sub$k_oxid[1])
    fv <- ps_oxid_formula(info$formula_vec, sub$k_oxid[1])
    ad_shift <- if (identical(sub$adduct[1], "H")) 0 else adduct_shift(sub$adduct[1])
    prec_mz <- (mass + ad_shift + h_offset - sub$z[1] * .PROTON) / sub$z[1]

    recs[[length(recs) + 1]] <- .spectrum_record(
      name = paste0(met$name, if (sub$k_oxid[1] > 0) paste0(" +", sub$k_oxid[1], "Ox") else "",
                   " [z=", sub$z[1], if (!identical(sub$adduct[1], "H")) paste0(" ", sub$adduct[1]) else "",
                   "] MS2 -- empirical"),
      precursor_mz = prec_mz, z = sub$z[1], formula = format_formula(fv),
      mono_mass = mass, level = "MS2",
      fields = c(met_id = met$id, met_kind = met$kind, k_oxid = as.character(sub$k_oxid[1]),
                adduct = sub$adduct[1], n_source_spectra = as.character(row$n_source_spectra),
                source_samples = row$source_samples),
      peaks = data.frame(mz = consensus$mz, intensity = consensus$intensity,
                        annotation = consensus$annotation, stringsAsFactors = FALSE))
  }

  recs <- .cap_records(recs, max_spectra, "Empirical MS2 spectral library",
                       "the number of confirmed batch MS2 hits")
  list(records = recs, summary = do.call(rbind, summary_rows))
}

# Pool real MS2 spectra (one data.frame(mz, intensity) per contributing
# scan) into one consensus spectrum: each spectrum is independently
# rescaled to its own 0-100 base peak, all peaks are pooled and sorted by
# mz, then walked once left-to-right, growing a cluster while the next
# peak sits within fragment_ppm of the cluster's running mean AND isn't a
# second peak from a spectrum already represented in it (two neighbouring
# peaks in the same acquired spectrum must never count as each other's
# "recurrence" vote). This is a simple greedy pass, not a full clustering
# solve, but it is adequate at fragment_ppm-scale tolerances.
.consensus_ms2_peaks <- function(spectra_list, fragment_ppm,
                                 min_consensus_fraction,
                                 min_consensus_peak_intensity) {
  n <- length(spectra_list)
  tagged <- lapply(seq_len(n), function(i) {
    sp <- spectra_list[[i]]
    sp$intensity <- sp$intensity / max(sp$intensity) * 100
    sp$spectrum_id <- i
    sp[, c("mz", "intensity", "spectrum_id")]
  })
  allpk <- do.call(rbind, tagged)
  allpk <- allpk[order(allpk$mz), ]

  clusters <- list()
  cur <- allpk[1, , drop = FALSE]
  cur_ids <- cur$spectrum_id
  if (nrow(allpk) > 1) {
    for (r in 2:nrow(allpk)) {
      row <- allpk[r, , drop = FALSE]
      cur_mean_mz <- mean(cur$mz)
      tol <- cur_mean_mz * fragment_ppm / 1e6
      if (abs(row$mz - cur_mean_mz) <= tol && !(row$spectrum_id %in% cur_ids)) {
        cur <- rbind(cur, row)
        cur_ids <- c(cur_ids, row$spectrum_id)
      } else {
        clusters[[length(clusters) + 1]] <- cur
        cur <- row
        cur_ids <- row$spectrum_id
      }
    }
  }
  clusters[[length(clusters) + 1]] <- cur

  out <- do.call(rbind, lapply(clusters, function(cl) {
    data.frame(mz = mean(cl$mz), intensity = mean(cl$intensity),
              n_spectra_with_peak = nrow(cl), stringsAsFactors = FALSE)
  }))
  out$n_spectra_total <- n
  out$consensus_fraction <- out$n_spectra_with_peak / n
  out <- out[out$n_spectra_total == 1 | out$consensus_fraction >= min_consensus_fraction, , drop = FALSE]
  if (nrow(out) == 0) return(out[, c("mz", "intensity")])

  out$intensity <- out$intensity / max(out$intensity) * 100
  out <- out[out$intensity >= min_consensus_peak_intensity * 100, , drop = FALSE]
  out <- out[order(out$mz), , drop = FALSE]
  rownames(out) <- NULL
  out[, c("mz", "intensity")]
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
# measured = TRUE labels records as acquired/annotated data (peaks are real
# instrument values -- pooled consensus peaks for build_empirical_ms2_library(),
# or a single acquired spectrum for batch_annotated_msp_records() in
# R/mirror_plot.R -- only the per-peak fragment-ion ANNOTATION is a
# prediction) rather than a fully theoretical spectrum.
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
