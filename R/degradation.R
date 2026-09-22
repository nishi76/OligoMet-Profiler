# =============================================================================
# degradation.R
# Peak-area-based % degradation of the parent oligonucleotide, computed per
# sample from annotate_metabolites_batch()'s ms1_matches (R/batch_ms_processing.R),
# grouped by the existing metabolite `kind` taxonomy from generate_metabolites()
# (R/metabolites.R): parent, exo_3p, exo_5p, endo_5frag, endo_3frag, and the
# mass-balanced endo_5frag_p/endo_3frag_p (the phosphorylated counterpart of
# each endonuclease product -- see endo_cleave()).
#
#   % degradation = 1 - (parent_signal / (parent_signal + degradant_signal))
#
# "signal" is peak area (AUC, threaded through match_ms1()/match_ms1_batch()
# in R/ms_matching.R -- see that module's header) by default, falling back
# to max intensity per-sample if area is entirely unavailable (e.g. the
# single-file R-native reading path, which does not compute AUC).
# =============================================================================

## ---- Collapse to one row per (sample, met_id) --------------------------------
# Same "max across charge/adduct/oxidation combos" rule build_abundance_matrix()
# uses, on whichever of area/intensity is being used as the signal.
.best_signal_per_met <- function(matches, signal_col) {
  stats::aggregate(stats::as.formula(paste(signal_col, "~ sample + met_id")),
                    data = matches, FUN = max)
}

## ---- Degradation summary ------------------------------------------------------
#' Summarize % degradation of the parent oligonucleotide, per sample
#'
#' `% degradation = 1 - (parent_signal / (parent_signal + degradant_signal))`,
#' computed per sample from [annotate_metabolites_batch()]'s
#' `ms1_matches`, grouped by the existing metabolite `kind` taxonomy from
#' [generate_metabolites()] (parent, exo_3p, exo_5p, endo_5frag,
#' endo_3frag, and their phosphorylated counterparts). Uses peak area
#' (AUC) as signal when available, falling back to max intensity
#' per-sample if area is entirely unavailable (e.g. the single-file
#' R-native reading path, which doesn't compute AUC).
#'
#' @param ms1_matches The `ms1_matches` data.frame from
#'   [annotate_metabolites_batch()]/[match_ms1_batch()].
#' @param top_n How many top degradant species to report per sample in
#'   `top_degradants`.
#' @param sample_meta Optional data.frame with `sample` plus `sample_type`
#'   (and, if available, `group`/`timepoint`). When given, calibration
#'   standards, quality-control, and reagent/matrix-blank samples (every
#'   `sample_type` other than `"unknown"`/blank/`NA`) are excluded before
#'   computing anything below -- they're for assessing assay performance,
#'   not the biological experiment, and shouldn't appear as their own
#'   "sample" in the composition-by-class chart or skew a parent's %
#'   degradation. Any `group`/`timepoint` columns present are also joined
#'   onto `per_sample`/`composition` so degradation can be tracked across
#'   the study design directly, without a separate manual merge. `NULL`
#'   (the default) keeps every sample, for callers with no sample-type
#'   information at all.
#' @return A list: `per_sample` (`sample`, `signal_used`, `parent_signal`,
#'   `degradant_signal`, `total_signal`, `pct_degradation`, plus `group`/
#'   `timepoint` if `sample_meta` supplied them); `composition`
#'   (`sample`, `kind`, `class_signal`, `pct_of_total`,
#'   `pct_of_degradants` -- grouped by the RAW `kind` value, e.g.
#'   `exo_3p`/`endo_5frag`, not pre-collapsed into a 5'/3'/endo framing,
#'   so a caller can still render either view; plus `group`/`timepoint` as
#'   above); `top_degradants`
#'   (`sample`, `met_id`, `met_name`, `kind`, `signal`,
#'   `pct_of_degradants`, `rank`); and `signal_used` (`"area"`,
#'   `"intensity"`, or `NA` if there's nothing to summarize).
#' @export
# Returns list(
#   per_sample    = data.frame(sample, signal_used, parent_signal,
#                               degradant_signal, total_signal, pct_degradation),
#   composition   = data.frame(sample, kind, class_signal, pct_of_total,
#                               pct_of_degradants),
#   top_degradants = data.frame(sample, met_id, met_name, kind, signal,
#                                pct_of_degradants, rank),
#   signal_used   = "area" | "intensity" | NA (NA if there's nothing to summarize)
# )
#
# `composition` groups by the raw `kind` value (exo_3p/exo_5p/endo_5frag/
# endo_3frag/endo_5frag_p/endo_3frag_p), not pre-collapsed into a 3-class
# "5' exo / 3' exo / endo"
# framing -- that collapse is a trivial display-layer ifelse() at render
# time (Shiny/Excel), not baked into this function, so callers keep the
# more informative breakdown and can still show either view.
degradation_summary <- function(ms1_matches, top_n = 10, sample_meta = NULL) {
  empty <- list(per_sample = data.frame(), composition = data.frame(),
                top_degradants = data.frame(), signal_used = NA_character_)
  if (is.null(ms1_matches) || nrow(ms1_matches) == 0) return(empty)

  signal_col <- if ("area" %in% names(ms1_matches) && any(!is.na(ms1_matches$area))) {
    "area"
  } else if ("intensity" %in% names(ms1_matches)) {
    "intensity"
  } else {
    return(empty)
  }
  m <- ms1_matches[!is.na(ms1_matches[[signal_col]]), ]

  if (!is.null(sample_meta) && "sample_type" %in% names(sample_meta)) {
    keep_samples <- sample_meta$sample[.is_study_sample(sample_meta$sample_type)]
    m <- m[m$sample %in% keep_samples, ]
  }
  if (nrow(m) == 0) return(empty)

  best <- .best_signal_per_met(m, signal_col)  # one row per (sample, met_id)
  best <- merge(best, unique(m[, c("met_id", "met_name", "kind")]), by = "met_id")

  samples <- sort(unique(best$sample))
  per_sample_rows <- list(); comp_rows <- list(); top_rows <- list()

  for (s in samples) {
    sub <- best[best$sample == s, ]
    parent_sig <- sum(sub[[signal_col]][sub$kind == "parent"], na.rm = TRUE)
    degr_sub <- sub[sub$kind != "parent", ]
    degr_sig <- sum(degr_sub[[signal_col]], na.rm = TRUE)
    total_sig <- parent_sig + degr_sig
    pct_degr <- if (total_sig > 0) 1 - (parent_sig / total_sig) else NA_real_

    per_sample_rows[[length(per_sample_rows) + 1]] <- data.frame(
      sample = s, signal_used = signal_col,
      parent_signal = parent_sig, degradant_signal = degr_sig,
      total_signal = total_sig,
      pct_degradation = if (is.na(pct_degr)) NA_real_ else round(pct_degr * 100, 2),
      stringsAsFactors = FALSE)

    for (k in unique(sub$kind)) {
      class_sig <- sum(sub[[signal_col]][sub$kind == k], na.rm = TRUE)
      comp_rows[[length(comp_rows) + 1]] <- data.frame(
        sample = s, kind = k, class_signal = class_sig,
        pct_of_total = if (total_sig > 0) round(100 * class_sig / total_sig, 2) else NA_real_,
        pct_of_degradants = if (k == "parent" || degr_sig == 0) NA_real_
                             else round(100 * class_sig / degr_sig, 2),
        stringsAsFactors = FALSE)
    }

    if (nrow(degr_sub) > 0 && degr_sig > 0) {
      ord <- degr_sub[order(-degr_sub[[signal_col]]), ]
      ord <- ord[seq_len(min(top_n, nrow(ord))), ]
      top_rows[[length(top_rows) + 1]] <- data.frame(
        sample = s, met_id = ord$met_id, met_name = ord$met_name, kind = ord$kind,
        signal = ord[[signal_col]],
        pct_of_degradants = round(100 * ord[[signal_col]] / degr_sig, 2),
        rank = seq_len(nrow(ord)), stringsAsFactors = FALSE)
    }
  }

  per_sample <- do.call(rbind, per_sample_rows)
  composition <- do.call(rbind, comp_rows)

  # Join group/timepoint (whichever is present) onto both tables so
  # degradation can be read off across the study design directly -- a
  # trend over time or between groups, not just an unordered list of
  # sample names -- without a second manual merge downstream.
  context_cols <- if (!is.null(sample_meta)) intersect(c("sample", "group", "timepoint"), names(sample_meta)) else "sample"
  if (length(context_cols) > 1) {
    ctx <- unique(sample_meta[, context_cols, drop = FALSE])
    per_sample <- merge(per_sample, ctx, by = "sample", all.x = TRUE)
    composition <- merge(composition, ctx, by = "sample", all.x = TRUE)
  }

  list(
    per_sample = per_sample,
    composition = composition,
    top_degradants = if (length(top_rows) > 0) do.call(rbind, top_rows) else data.frame(),
    signal_used = signal_col
  )
}

## ---- Composition plot ---------------------------------------------------------
# Stacked bar of % of total signal by kind, one bar per sample -- matches
# the multi-category palette convention used elsewhere (R/statistics.R's
# plot_group_boxplot()).
plot_degradation_composition <- function(degradation) {
  comp <- degradation$composition
  if (is.null(comp) || nrow(comp) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = "No degradation data to plot"))
  }
  kinds <- sort(unique(comp$kind))
  ggplot2::ggplot(comp, ggplot2::aes(x = .data$sample, y = .data$pct_of_total, fill = .data$kind)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = grDevices::colorRampPalette(
      c("#0279EE", "#FF9400", "#75A025", "#FD9BED", "#E9ED4C"))(length(kinds)),
      name = "Kind") +
    ggplot2::labs(x = NULL, y = "% of total signal",
                  title = "Composition by metabolite class",
                  subtitle = paste0("Signal: ", degradation$signal_used %||% "n/a")) +
    ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
