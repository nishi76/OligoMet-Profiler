# =============================================================================
# build_workbook.R
# Excel workbook generation for oligonucleotide metabolite identification.
#
# Produces a multi-sheet Excel workbook (.xlsx) with:
#   Sheet 1: Metabolite Library     — one row per metabolite
#   Sheet 2: Charge Envelopes       — one row per (metabolite, z, isotope)
#   Sheet 3: PS Oxidation Series    — one row per (metabolite, k_oxid, z)
#   Sheet 4: Fragment Ions          — one row per (metabolite, ion, z)
#   Sheet 5: PRM Inclusion List     — for MS acquisition
#   Sheet 6: MS Matching Results    — if MS data provided
#   Sheet 7: Summary                — overview and validation
#
# Uses openxlsx for formatting (headers, column widths, freeze panes).
#
# Performance note: every sheet is assembled as a single data.frame and
# written with ONE openxlsx::writeData() call, and cell-level styling is
# limited to the header row plus a handful of accent cells. The base font
# for all data cells comes from modifyBaseFont() in build_workbook(), so no
# per-cell "normal" style is needed. An earlier version wrote row by row
# and attached a style object to every data cell, which made saveWorkbook()
# serialize hundreds of thousands of styled cells one at a time -- that is
# what made the final Excel export take longer than the entire computation.
# =============================================================================

## ---- Style definitions ----------------------------------------------------
.wb_styles <- function() {
  list(
    header = openxlsx::createStyle(
      fontName = "Arial", fontSize = 11, fontColour = "#FFFFFF",
      fgFill = "#0279EE", halign = "center", textDecoration = "bold",
      border = "TopBottomLeftRight", borderColour = "#FFFFFF", borderStyle = "thin"),
    subheader = openxlsx::createStyle(
      fontName = "Arial", fontSize = 10, fontColour = "#333333",
      fgFill = "#ECE9E2", halign = "center", textDecoration = "bold",
      border = "TopBottomLeftRight", borderColour = "#CCCCCC", borderStyle = "thin"),
    parent = openxlsx::createStyle(
      fontName = "Arial", fontSize = 10, fgFill = "#FAF9F3",
      textDecoration = "bold"),
    normal = openxlsx::createStyle(
      fontName = "Arial", fontSize = 10),
    mono = openxlsx::createStyle(
      fontName = "Consolas", fontSize = 9),
    good = openxlsx::createStyle(
      fontName = "Arial", fontSize = 10, fontColour = "#75A025"),
    warn = openxlsx::createStyle(
      fontName = "Arial", fontSize = 10, fontColour = "#FF9400")
  )
}

## ---- Shared sheet scaffolding ----------------------------------------------
# Write custom header text on row 1, the data block below it (one call), and
# apply the header style + column widths + freeze pane.
.write_sheet <- function(wb, sheet, headers, df, styles, widths) {
  openxlsx::writeData(wb, sheet, as.data.frame(t(headers)),
                       startRow = 1, colNames = FALSE)
  openxlsx::addStyle(wb, sheet, styles$header,
                      rows = 1, cols = seq_along(headers), gridExpand = TRUE)
  if (!is.null(df) && nrow(df) > 0) {
    openxlsx::writeData(wb, sheet, df, startRow = 2, colNames = FALSE)
  }
  openxlsx::setColWidths(wb, sheet, cols = seq_along(widths), widths = widths)
  openxlsx::freezePane(wb, sheet, firstActiveRow = 2)
}

## ---- Sheet 1: Metabolite Library ------------------------------------------
.build_metabolite_sheet <- function(wb, mets, dict, styles) {
  openxlsx::addWorksheet(wb, "Metabolite Library")
  headers <- c("ID", "Name", "Kind", "Length", "Bases (5'→3')", "Sugars",
               "Linkages", "n_PS", "n_PO", "Conj_5", "Conj_3",
               "Formula", "Monoisotopic Mass (Da)", "Average Mass (Da)",
               "Modification", "Site")

  rows <- lapply(mets, function(met) {
    info <- metabolite_mass_info(met, dict)
    data.frame(
      id = met$id, name = met$name, kind = met$kind, n = met$n,
      bases = paste(met$bases, collapse = ""),
      sugars = paste(met$sugars, collapse = ""),
      linkages = paste(ifelse(is.na(met$linkages), ".", met$linkages),
                       collapse = ""),
      n_ps = met$n_ps, n_po = met$n_po,
      conj5 = met$conj5, conj3 = met$conj3,
      formula = info$formula_str,
      mono_mass = round(info$mono_mass, 6),
      avg_mass = round(info$avg_mass, 4),
      modification = met$modification,
      site = ifelse(is.null(met$site) || is.na(met$site), "", met$site),
      stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)

  widths <- c(6, 22, 12, 6, 18, 16, 16, 6, 6, 10, 10, 28, 18, 16, 28, 6)
  .write_sheet(wb, "Metabolite Library", headers, df, styles, widths)

  # Accent styling: highlight parent rows, monospace the formula column.
  parent_rows <- which(vapply(mets, function(m) m$kind == "parent",
                              logical(1))) + 1
  if (length(parent_rows) > 0) {
    openxlsx::addStyle(wb, "Metabolite Library", styles$parent,
                        rows = parent_rows, cols = seq_along(headers),
                        gridExpand = TRUE)
  }
  openxlsx::addStyle(wb, "Metabolite Library", styles$mono,
                      rows = seq_len(nrow(df)) + 1, cols = 12,
                      gridExpand = TRUE)
}

## ---- Sheet 2: Charge Envelopes --------------------------------------------
# This is the most expensive sheet to compute: one isotope-pattern
# calculation per (metabolite, PS-oxidation level) -- see compute_envelope().
# If a progress_utils.R tracker is supplied (console_tracker), progress is
# shown as a live in-place updating bar with a step-local ETA
# (progress_tick()); otherwise falls back to a plain text line every 10
# metabolites so this still works if progress_utils.R hasn't been sourced.
# `progress()` is a separate, optional callback used by the Shiny app to
# nudge its own browser-side progress bar; both can be supplied at once.
.build_envelope_sheet <- function(wb, mets, dict, styles,
                                   z_range, n_iso, max_oxid, h_offset,
                                   use_envipat, progress = NULL,
                                   console_tracker = NULL) {
  openxlsx::addWorksheet(wb, "Charge Envelopes")
  headers <- c("Met ID", "Met Name", "k_Oxid", "Formula", "z", "iso",
               "m/z", "Abundance", "Monoisotopic Mass (Da)")

  n_mets <- length(mets)
  t0 <- Sys.time()
  use_tick <- !is.null(console_tracker) && exists("progress_tick")
  chunks <- vector("list", n_mets)
  for (i in seq_along(mets)) {
    met <- mets[[i]]
    env <- compute_envelope(met, z_range = z_range, n_iso = n_iso,
                             max_oxid = max_oxid, h_offset = h_offset,
                             use_envipat = use_envipat, dict = dict)
    if (!is.null(env) && nrow(env) > 0) chunks[[i]] <- env
    if (use_tick) {
      progress_tick(console_tracker, i, n_mets, "Charge envelopes")
      if (!is.null(progress) && (i %% 10 == 0 || i == n_mets)) {
        progress(sprintf("Charge envelopes: %d/%d metabolites", i, n_mets))
      }
    } else if (i %% 10 == 0 || i == n_mets) {
      elapsed <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)
      msg <- sprintf("Charge envelopes: %d/%d metabolites (%.1fs elapsed)",
                      i, n_mets, elapsed)
      cat(" ", msg, "\n")
      if (!is.null(progress)) progress(msg)
    }
  }
  if (use_tick && exists("progress_tick_end")) progress_tick_end()

  df <- do.call(rbind, chunks)
  if (!is.null(df) && nrow(df) > 0) {
    # Column order must match the header row above (an earlier version
    # wrote compute_envelope()'s native order, which put z under the
    # "Formula" header and shifted every column after it).
    df <- df[, c("met_id", "met_name", "k_oxid", "formula", "z", "iso",
                 "mz", "abundance", "mono_mass")]
  }

  widths <- c(6, 22, 8, 28, 5, 5, 14, 12, 18)
  .write_sheet(wb, "Charge Envelopes", headers, df, styles, widths)
}

## ---- Sheet 3: PS Oxidation Series -----------------------------------------
.build_oxidation_sheet <- function(wb, mets, dict, styles, max_oxid, z_range,
                                    h_offset) {
  openxlsx::addWorksheet(wb, "PS Oxidation Series")
  headers <- c("Met ID", "Met Name", "k_Oxid", "Formula",
               "Monoisotopic Mass (Da)", "Average Mass (Da)",
               paste0("z=", z_range))

  chunks <- lapply(mets, function(met) {
    series <- ps_oxidation_series(met, max_oxid = max_oxid, dict = dict)
    series <- series[series$k <= min(met$n_ps, max_oxid), , drop = FALSE]
    mz <- vapply(z_range,
                 function(z) round((series$mono_mass + h_offset - z * .PROTON) / z, 4),
                 numeric(nrow(series)))
    mz <- matrix(mz, nrow = nrow(series))
    colnames(mz) <- paste0("z", z_range)
    cbind(
      data.frame(met_id = met$id, met_name = met$name, k = series$k,
                 formula = series$formula_str,
                 mono_mass = round(series$mono_mass, 6),
                 avg_mass = round(series$avg_mass, 4),
                 stringsAsFactors = FALSE),
      as.data.frame(mz))
  })
  df <- do.call(rbind, chunks)

  widths <- c(6, 22, 8, 28, 18, 16, rep(12, length(z_range)))
  .write_sheet(wb, "PS Oxidation Series", headers, df, styles, widths)
  if (!is.null(df) && nrow(df) > 0) {
    openxlsx::addStyle(wb, "PS Oxidation Series", styles$mono,
                        rows = seq_len(nrow(df)) + 1, cols = 4,
                        gridExpand = TRUE)
  }
}

## ---- Sheet 4: Fragment Ions -----------------------------------------------
.build_fragment_sheet <- function(wb, mets, dict, styles, z_range,
                                   include_internal = TRUE) {
  openxlsx::addWorksheet(wb, "Fragment Ions")
  headers <- c("Met ID", "Met Name", "Ion Type", "Direction", "Cleavage Site",
               "Fragment Length", "Formula", "Monoisotopic Mass (Da)",
               "Base Loss", "z", "m/z")

  chunks <- lapply(mets, function(met) {
    if (met$n < 3) return(NULL)
    frags <- generate_fragments(met, dict, z_range = z_range)
    if (include_internal) {
      frags <- c(frags, generate_internal_fragments(met, dict, z_range = z_range))
    }
    if (length(frags) == 0) return(NULL)
    base <- data.frame(
      met_id = vapply(frags, function(f) as.character(f$met_id), character(1)),
      met_name = met$name,
      ion_type = vapply(frags, function(f) f$ion_type, character(1)),
      direction = vapply(frags, function(f) f$direction, character(1)),
      cleavage_site = vapply(frags, function(f) as.numeric(f$cleavage_site), numeric(1)),
      frag_length = vapply(frags, function(f) as.numeric(f$frag_length), numeric(1)),
      formula = vapply(frags, function(f) f$formula, character(1)),
      mono_mass = round(vapply(frags, function(f) f$mono_mass, numeric(1)), 4),
      base_loss = vapply(frags, function(f)
        ifelse(is.na(f$base_loss), "", as.character(f$base_loss)), character(1)),
      stringsAsFactors = FALSE)
    # Expand: one row per (fragment, charge state).
    idx <- rep(seq_len(nrow(base)), each = length(z_range))
    out <- base[idx, , drop = FALSE]
    out$z <- rep(z_range, times = nrow(base))
    out$mz <- round((out$mono_mass - out$z * .PROTON) / out$z, 4)
    out
  })
  df <- do.call(rbind, chunks)

  widths <- c(6, 22, 10, 10, 10, 10, 28, 18, 10, 5, 14)
  .write_sheet(wb, "Fragment Ions", headers, df, styles, widths)
}

## ---- Sheet 5: PRM Inclusion List ------------------------------------------
.build_prm_sheet <- function(wb, mets, dict, styles, z_range, max_oxid, h_offset) {
  openxlsx::addWorksheet(wb, "PRM Inclusion List")
  prm <- prm_inclusion_list(mets, dict, z_range = z_range,
                             max_oxid = max_oxid, h_offset = h_offset)
  if (nrow(prm) == 0) {
    openxlsx::writeData(wb, "PRM Inclusion List", "No entries")
    return()
  }
  openxlsx::writeData(wb, "PRM Inclusion List", prm,
                       startRow = 1, colNames = TRUE)
  openxlsx::addStyle(wb, "PRM Inclusion List", styles$header,
                      rows = 1, cols = 1:ncol(prm), gridExpand = TRUE)
  widths <- c(6, 22, 12, 5, 8, 5, 14, 12)
  openxlsx::setColWidths(wb, "PRM Inclusion List", cols = 1:length(widths),
                          widths = widths)
  openxlsx::freezePane(wb, "PRM Inclusion List", firstActiveRow = 2)
}

## ---- Sheet 6: MS Matching Results -----------------------------------------
.build_matching_sheet <- function(wb, results, styles) {
  if (is.null(results) || nrow(results$ms1_matches) == 0) {
    openxlsx::addWorksheet(wb, "MS Matching")
    openxlsx::writeData(wb, "MS Matching", "No MS data provided or no matches found")
    return()
  }
  openxlsx::addWorksheet(wb, "MS Matching")
  # MS1 matches
  openxlsx::writeData(wb, "MS Matching", "MS1 Matches",
                       startRow = 1, colNames = FALSE)
  openxlsx::addStyle(wb, "MS Matching", styles$subheader,
                      rows = 1, cols = 1:12, gridExpand = TRUE)
  openxlsx::writeData(wb, "MS Matching", results$ms1_matches,
                       startRow = 2, colNames = TRUE)
  openxlsx::addStyle(wb, "MS Matching", styles$header,
                      rows = 2, cols = 1:ncol(results$ms1_matches), gridExpand = TRUE)
  n_ms1 <- nrow(results$ms1_matches)

  # Summary
  if (!is.null(results$summary) && nrow(results$summary) > 0) {
    srow <- n_ms1 + 4
    openxlsx::writeData(wb, "MS Matching", "Annotation Summary",
                         startRow = srow, colNames = FALSE)
    openxlsx::addStyle(wb, "MS Matching", styles$subheader,
                        rows = srow, cols = 1:12, gridExpand = TRUE)
    openxlsx::writeData(wb, "MS Matching", results$summary,
                         startRow = srow + 1, colNames = TRUE)
    openxlsx::addStyle(wb, "MS Matching", styles$header,
                        rows = srow + 1, cols = 1:ncol(results$summary), gridExpand = TRUE)
  }

  openxlsx::setColWidths(wb, "MS Matching", cols = 1:15, widths = 14)
}

## ---- Sheet 7: Summary -----------------------------------------------------
.build_summary_sheet <- function(wb, spec, mets, dict, validation, styles,
                                  opts, ms_info) {
  openxlsx::addWorksheet(wb, "Summary")
  st <- styles$subheader
  nm <- styles$normal

  openxlsx::writeData(wb, "Summary", "Oligonucleotide Metabolite Identification",
                       startRow = 1, colNames = FALSE)
  openxlsx::addStyle(wb, "Summary",
    openxlsx::createStyle(fontSize = 14, textDecoration = "bold", fontColour = "#0279EE"),
    rows = 1, cols = 1)

  row <- 3
  write_pair <- function(label, value, row) {
    openxlsx::writeData(wb, "Summary", label, startRow = row, colNames = FALSE)
    openxlsx::writeData(wb, "Summary", as.character(value), startRow = row,
                         startCol = 2, colNames = FALSE)
    openxlsx::addStyle(wb, "Summary", st, rows = row, cols = 1)
    openxlsx::addStyle(wb, "Summary", nm, rows = row, cols = 2)
    row + 1
  }

  row <- write_pair("Oligonucleotide", spec$raw %||% "N/A", row)
  row <- write_pair("Notation", spec$notation, row)
  row <- write_pair("Length (nt)", spec$n, row)
  row <- write_pair("Bases (5'→3')", paste(spec$bases, collapse = ""), row)
  row <- write_pair("Sugars", paste(spec$sugars, collapse = ""), row)
  row <- write_pair("Linkages", paste(ifelse(is.na(spec$linkages), ".", spec$linkages), collapse = ""), row)
  row <- write_pair("5' Conjugate", spec$conj5, row)
  row <- write_pair("3' Conjugate", spec$conj3, row)

  parent_info <- metabolite_mass_info(mets[[1]], dict)
  row <- row + 1
  row <- write_pair("Parent Formula", parent_info$formula_str, row)
  row <- write_pair("Parent Monoisotopic Mass (Da)", round(parent_info$mono_mass, 6), row)
  row <- write_pair("Parent Average Mass (Da)", round(parent_info$avg_mass, 4), row)

  row <- row + 1
  row <- write_pair("Total Metabolites Generated", length(mets), row)
  row <- write_pair("  Parent", sum(sapply(mets, function(m) m$kind == "parent")), row)
  row <- write_pair("  3' Exonuclease", sum(sapply(mets, function(m) m$kind == "exo_3p")), row)
  row <- write_pair("  5' Exonuclease", sum(sapply(mets, function(m) m$kind == "exo_5p")), row)
  row <- write_pair("  Endonuclease", sum(sapply(mets, function(m) grepl("endo", m$kind))), row)

  row <- row + 1
  row <- write_pair("Max 3' Truncation", opts$max_3p, row)
  row <- write_pair("Max 5' Truncation", opts$max_5p, row)
  row <- write_pair("Endonuclease", opts$endo, row)
  row <- write_pair("Max PS Oxidation", opts$max_oxid, row)
  row <- write_pair("Charge State Range", paste(range(opts$z_range), collapse = " - "), row)
  row <- write_pair("Adducts", paste(opts$adducts, collapse = ", "), row)
  row <- write_pair("Mass Tolerance (ppm)", opts$ppm_tol, row)

  if (!is.null(validation)) {
    row <- row + 1
    row <- write_pair("Validation: Formula Match", validation$ok, row)
    row <- write_pair("Validation: Mass Delta (ppm)", round(validation$ppm, 4), row)
  }

  if (!is.null(ms_info) && !is.na(ms_info$file)) {
    row <- row + 1
    row <- write_pair("MS Data File", ms_info$file, row)
    row <- write_pair("MS1 Spectra", ms_info$n_ms1, row)
    row <- write_pair("MS2 Spectra", ms_info$n_ms2, row)
  }

  openxlsx::setColWidths(wb, "Summary", cols = 1:2, widths = c(30, 50))
}

## ---- Main workbook builder -------------------------------------------------
# Build the complete Excel workbook and save to file. Writes to output_dir
# (a scratch/temp location by default) -- callers that want the file
# somewhere specific should copy it from the returned path, since binary
# formats like xlsx write more reliably to a local scratch dir first.
build_workbook <- function(spec, mets, dict = STANDARD_DICT,
                           validation = NULL, ms_results = NULL,
                           ms_info = NULL, opts = list(),
                           output_file = "oligo_metabolite_library.xlsx",
                           output_dir = tempdir(), progress = NULL,
                           console_tracker = NULL) {
  z_range <- opts$z_range %||% 3:12
  n_iso <- opts$n_iso %||% 5
  max_oxid <- opts$max_oxid %||% 6
  h_offset <- opts$h_offset %||% 0
  use_envipat <- opts$use_envipat %||% TRUE
  include_internal <- opts$include_internal %||% TRUE

  styles <- .wb_styles()
  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 10)

  # Build sheets
  .build_summary_sheet(wb, spec, mets, dict, validation, styles, opts, ms_info)
  .build_metabolite_sheet(wb, mets, dict, styles)
  .build_envelope_sheet(wb, mets, dict, styles, z_range, n_iso, max_oxid,
                         h_offset, use_envipat, progress = progress,
                         console_tracker = console_tracker)
  .build_oxidation_sheet(wb, mets, dict, styles, max_oxid, z_range, h_offset)
  .build_fragment_sheet(wb, mets, dict, styles, 1:2, include_internal)
  .build_prm_sheet(wb, mets, dict, styles, z_range, max_oxid, h_offset)
  .build_matching_sheet(wb, ms_results, styles)

  # Save (write to a scratch dir first for binary format, then copy)
  out_path <- file.path(output_dir, output_file)
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  cat("Workbook saved to", out_path, "\n")
  out_path
}
