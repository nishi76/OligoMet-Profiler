# =============================================================================
# blank_correction.R
# Reagent/matrix-blank background subtraction. A blank sample (identified
# via `sample_type` OR an RB/MB token in the sample name -- see
# .is_blank_sample()/.detect_blank_type_from_name() in R/chemistry_dict.R)
# carries no biological signal by definition, so its per-metabolite signal
# is instrument/reagent/matrix background. Subtracting the pooled mean
# blank signal for each metabolite from every other sample's own signal
# gives a background-corrected estimate.
#
#   corrected = max(0, signal - mean(blank signal for that metabolite))
#
# Reagent and matrix blanks are pooled into a single per-metabolite value
# rather than kept separate -- there is no per-sample-family assignment of
# which blank "belongs" to which study sample in this codebase's sample
# metadata, so pooling is the only scope that doesn't invent one.
#
# Raw intensity/area are never modified. This module only ever ADDS new
# `intensity_bcorr`/`area_bcorr` columns alongside them (see
# apply_blank_correction()), so every existing consumer of `ms1_matches`
# is unaffected unless it explicitly opts into the corrected column.
# =============================================================================

## ---- Per-metabolite blank value ---------------------------------------------
#' Compute the pooled per-metabolite blank signal
#'
#' Averages a metabolite's raw signal across every row whose sample is a
#' reagent or matrix blank (see [.is_blank_sample()]), pooling both blank
#' types into a single value per metabolite. A metabolite with zero blank
#' observations gets `blank_mean = 0` (so [apply_blank_correction()]
#' leaves it unchanged), rather than being dropped.
#'
#' @param batch_matches The `ms1_matches` data.frame from
#'   [match_ms1_batch()]/[annotate_metabolites_batch()].
#' @param sample_meta A data.frame with `sample` and `sample_type`
#'   columns (used together with the sample name itself to identify
#'   blanks -- see [.is_blank_sample()]).
#' @param signal_col Which column to summarize. `NULL` (the default)
#'   auto-detects: `"area"` if present and not all-`NA`, else
#'   `"intensity"`.
#' @return A list: `table` (a data.frame with `met_id`, `met_name`,
#'   `blank_mean`, `blank_sd`, `n_blank` -- one row per metabolite present
#'   in `batch_matches`), and `note` (empty string on success, else why no
#'   blank samples were found).
#' @seealso [apply_blank_correction()]
#' @export
compute_blank_signal <- function(batch_matches, sample_meta, signal_col = NULL) {
  empty <- list(table = data.frame(), note = "")
  if (is.null(batch_matches) || nrow(batch_matches) == 0) {
    empty$note <- "no MS1 matches to compute a blank value from"
    return(empty)
  }
  if (is.null(signal_col)) signal_col <- .auto_signal_col(batch_matches)
  if (!signal_col %in% names(batch_matches) || all(is.na(batch_matches[[signal_col]]))) {
    empty$note <- paste0("no '", signal_col, "' signal to compute a blank value from")
    return(empty)
  }
  if (is.null(sample_meta) || !all(c("sample", "sample_type") %in% names(sample_meta))) {
    empty$note <- "sample_meta has no sample/sample_type column"
    return(empty)
  }

  blank_samples <- sample_meta$sample[.is_blank_sample(sample_meta$sample_type, sample_meta$sample)]
  m <- batch_matches[batch_matches$sample %in% blank_samples & !is.na(batch_matches[[signal_col]]), ]

  met_info <- unique(batch_matches[, c("met_id", "met_name")])
  met_info <- met_info[order(met_info$met_id), ]

  if (nrow(m) == 0) {
    empty$note <- "no reagent_blank/matrix_blank samples found (by Sample Type or RB/MB in name)"
    empty$table <- data.frame(met_id = met_info$met_id, met_name = met_info$met_name,
                               blank_mean = 0, blank_sd = NA_real_, n_blank = 0L,
                               stringsAsFactors = FALSE)
    return(empty)
  }

  best <- .best_signal_per_met(m, signal_col)  # one row per (blank sample, met_id)
  agg_mean <- stats::aggregate(stats::as.formula(paste(signal_col, "~ met_id")), data = best, FUN = mean)
  agg_sd <- stats::aggregate(stats::as.formula(paste(signal_col, "~ met_id")), data = best, FUN = stats::sd)
  agg_n <- stats::aggregate(stats::as.formula(paste(signal_col, "~ met_id")), data = best, FUN = length)

  tbl <- data.frame(met_id = met_info$met_id, met_name = met_info$met_name, stringsAsFactors = FALSE)
  tbl$blank_mean <- agg_mean[[signal_col]][match(tbl$met_id, agg_mean$met_id)]
  tbl$blank_mean[is.na(tbl$blank_mean)] <- 0
  tbl$blank_sd <- agg_sd[[signal_col]][match(tbl$met_id, agg_sd$met_id)]
  tbl$n_blank <- agg_n[[signal_col]][match(tbl$met_id, agg_n$met_id)]
  tbl$n_blank[is.na(tbl$n_blank)] <- 0L

  list(table = tbl, note = "")
}

## ---- Apply the correction ----------------------------------------------------
#' Add blank-corrected signal columns to `ms1_matches`
#'
#' Computes the pooled per-metabolite blank value (see
#' [compute_blank_signal()]) for `intensity` and/or `area` (whichever are
#' present in `batch_matches`) and adds `intensity_bcorr`/`area_bcorr`
#' columns: `max(0, signal - blank_mean)` for that metabolite, floored at
#' zero. If no reagent/matrix blank samples are found at all, the input
#' is returned unchanged (no `_bcorr` columns added), so callers can check
#' `"area_bcorr" %in% names(m)` to know whether blank correction is
#' available.
#'
#' @param batch_matches The `ms1_matches` data.frame from
#'   [match_ms1_batch()]/[annotate_metabolites_batch()].
#' @param sample_meta A data.frame with `sample` and `sample_type` columns.
#' @param signal_col Which column(s) to correct. `NULL` (the default)
#'   corrects every one of `intensity`/`area` present in `batch_matches`.
#'   A single column name corrects only that one.
#' @return `batch_matches` with `intensity_bcorr`/`area_bcorr` columns
#'   added when blanks were found; unchanged otherwise.
#' @seealso [compute_blank_signal()]
#' @export
apply_blank_correction <- function(batch_matches, sample_meta, signal_col = NULL) {
  if (is.null(batch_matches) || nrow(batch_matches) == 0) return(batch_matches)

  cols <- if (!is.null(signal_col)) signal_col else intersect(c("intensity", "area"), names(batch_matches))
  cols <- cols[vapply(cols, function(cc) cc %in% names(batch_matches) && any(!is.na(batch_matches[[cc]])), logical(1))]
  if (length(cols) == 0) return(batch_matches)

  out <- batch_matches
  for (cc in cols) {
    blank <- compute_blank_signal(batch_matches, sample_meta, signal_col = cc)
    if (nzchar(blank$note) && all(blank$table$n_blank == 0)) next  # no blanks at all -- don't add the column
    blank_mean <- blank$table$blank_mean[match(out$met_id, blank$table$met_id)]
    blank_mean[is.na(blank_mean)] <- 0
    out[[paste0(cc, "_bcorr")]] <- pmax(0, out[[cc]] - blank_mean)
  }
  out
}
