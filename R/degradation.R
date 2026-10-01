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
#' @param signal_col Which column to use as signal (e.g. a blank-corrected
#'   column such as `"intensity_bcorr"` -- see [apply_blank_correction()]).
#'   `NULL` (the default) auto-detects: `"area"` if present and not
#'   all-`NA`, else `"intensity"`.
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
degradation_summary <- function(ms1_matches, top_n = 10, sample_meta = NULL, signal_col = NULL) {
  empty <- list(per_sample = data.frame(), composition = data.frame(),
                top_degradants = data.frame(), signal_used = NA_character_)
  if (is.null(ms1_matches) || nrow(ms1_matches) == 0) return(empty)

  if (is.null(signal_col)) {
    signal_col <- if ("area" %in% names(ms1_matches) && any(!is.na(ms1_matches$area))) {
      "area"
    } else if ("intensity" %in% names(ms1_matches)) {
      "intensity"
    } else {
      return(empty)
    }
  } else if (!signal_col %in% names(ms1_matches) || all(is.na(ms1_matches[[signal_col]]))) {
    return(empty)
  }
  m <- ms1_matches[!is.na(ms1_matches[[signal_col]]), ]
  # An internal standard is neither parent nor degradant.
  if ("kind" %in% names(m)) m <- m[is.na(m$kind) | m$kind != "internal_standard", ]

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
plot_degradation_composition <- function(degradation, by = c("sample", "condition")) {
  by <- match.arg(by)
  comp <- degradation$composition
  if (is.null(comp) || nrow(comp) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = "No degradation data to plot"))
  }
  pal <- function(k) grDevices::colorRampPalette(c("#0279EE", "#FF9400", "#75A025", "#FD9BED", "#E9ED4C"))(length(k))
  theme <- ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans")
  sub <- paste0("Signal: ", degradation$signal_used %||% "n/a")

  if (by == "condition") {
    bc <- degradation_by_condition(degradation)
    cc <- bc$composition
    if (is.null(cc) || nrow(cc) == 0) {
      return(ggplot2::ggplot() + ggplot2::theme_void() +
               ggplot2::labs(title = "No Group/Timepoint in the sample table to group by"))
    }
    kinds <- sort(unique(cc$kind))
    x_var <- if (bc$design == "group") "group" else "timepoint_lab"
    p <- ggplot2::ggplot(cc, ggplot2::aes(x = .data[[x_var]], y = .data$mean_pct_of_total, fill = .data$kind)) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::scale_fill_manual(values = pal(kinds), name = "Kind") +
      ggplot2::labs(x = if (x_var == "group") "Group" else "Timepoint", y = "Mean % of total signal",
                    title = "Composition by metabolite class, per group / timepoint",
                    subtitle = paste0(sub, "; mean of replicates (n per bar in the table)")) + theme
    if (bc$design == "group_timepoint") p <- p + ggplot2::facet_wrap(~ group, nrow = 1)
    return(p)
  }

  kinds <- sort(unique(comp$kind))
  has_g <- "group" %in% names(comp) && any(!is.na(comp$group) & nzchar(comp$group))
  has_t <- "timepoint" %in% names(comp) && any(!is.na(suppressWarnings(as.numeric(comp$timepoint))))
  if (has_g || has_t) {
    tn <- if ("timepoint" %in% names(comp)) suppressWarnings(as.numeric(comp$timepoint)) else NA_real_
    gg <- if ("group" %in% names(comp)) ifelse(is.na(comp$group), "", comp$group) else ""
    lab <- .condition_label(if (has_g) comp$group else NULL, if (has_t) comp$timepoint else NULL)
    ord <- order(gg, tn)
    comp$.cond <- factor(lab, levels = unique(lab[ord]))
  }
  p <- ggplot2::ggplot(comp, ggplot2::aes(x = .data$sample, y = .data$pct_of_total, fill = .data$kind)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = pal(kinds), name = "Kind") +
    ggplot2::labs(x = NULL, y = "% of total signal",
                  title = "Composition by metabolite class, per sample", subtitle = sub) + theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  # Samples sit inside their own group/timepoint panel instead of one
  # alphabetical row, so replicates of a condition stay together.
  if (has_g || has_t) {
    p <- p + ggplot2::facet_grid(~ .cond, scales = "free_x", space = "free_x") +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 9))
  }
  p
}

## ---- Grouping by experimental design --------------------------------------------

.condition_label <- function(group, timepoint) {
  g <- if (is.null(group)) rep("", length(timepoint)) else ifelse(is.na(group), "", as.character(group))
  t <- if (is.null(timepoint)) rep("", length(g)) else ifelse(is.na(timepoint) | !nzchar(as.character(timepoint)), "",
                                                             paste0("t=", timepoint))
  lab <- trimws(paste(g, t))
  ifelse(nzchar(lab), lab, "(no group/timepoint)")
}

#' Degradation summarised by experimental group and timepoint
#'
#' Collapses [degradation_summary()]'s per-sample tables to one row per
#' condition (group, timepoint, or group x timepoint -- whichever the
#' sample table carries), so replicates are reported as mean, SD and n.
#' A class absent from a replicate counts as 0% for that replicate.
#'
#' @param degradation A [degradation_summary()] result computed with
#'   `sample_meta`, so `per_sample`/`composition` carry `group`/`timepoint`.
#' @return A list: `design` (`"group"`, `"timepoint"`, `"group_timepoint"`,
#'   or `"none"`); `per_condition` (`group`, `timepoint`, `n`,
#'   `mean_pct_degradation`, `sd_pct_degradation`, `sem_pct_degradation`,
#'   `mean_parent_signal`, `mean_degradant_signal`); `composition`
#'   (`group`, `timepoint`, `kind`, `n`, `mean_pct_of_total`,
#'   `sd_pct_of_total`, `mean_pct_of_degradants`).
#' @export
degradation_by_condition <- function(degradation) {
  empty <- list(design = "none", per_condition = data.frame(), composition = data.frame())
  ps <- degradation$per_sample
  if (is.null(ps) || nrow(ps) == 0) return(empty)
  has_g <- "group" %in% names(ps) && any(!is.na(ps$group) & nzchar(ps$group))
  tnum <- if ("timepoint" %in% names(ps)) suppressWarnings(as.numeric(ps$timepoint)) else rep(NA_real_, nrow(ps))
  has_t <- any(!is.na(tnum))
  if (!has_g && !has_t) return(empty)
  design <- if (has_g && has_t) "group_timepoint" else if (has_g) "group" else "timepoint"

  ps$group <- if (has_g) ifelse(is.na(ps$group), "", as.character(ps$group)) else ""
  ps$timepoint <- if (has_t) tnum else NA_real_
  keep <- (if (has_g) nzchar(ps$group) else TRUE) & (if (has_t) !is.na(ps$timepoint) else TRUE)
  ps <- ps[keep, , drop = FALSE]
  if (nrow(ps) == 0) return(empty)
  key <- paste(ps$group, ps$timepoint, sep = "\r")
  sd0 <- function(x) if (sum(!is.na(x)) > 1) stats::sd(x, na.rm = TRUE) else NA_real_

  per <- do.call(rbind, lapply(split(ps, key), function(d) {
    n <- sum(!is.na(d$pct_degradation))
    data.frame(group = d$group[1], timepoint = d$timepoint[1], n = nrow(d),
               mean_pct_degradation = mean(d$pct_degradation, na.rm = TRUE),
               sd_pct_degradation = sd0(d$pct_degradation),
               sem_pct_degradation = if (n > 1) sd0(d$pct_degradation) / sqrt(n) else NA_real_,
               mean_parent_signal = mean(d$parent_signal, na.rm = TRUE),
               mean_degradant_signal = mean(d$degradant_signal, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  per <- per[order(per$group, per$timepoint), , drop = FALSE]
  rownames(per) <- NULL

  comp <- degradation$composition
  comp <- comp[comp$sample %in% ps$sample, , drop = FALSE]
  kinds <- sort(unique(comp$kind))
  # Complete every sample x kind with 0% so a class missing from one
  # replicate pulls that condition's mean down instead of being ignored.
  grid <- expand.grid(sample = unique(ps$sample), kind = kinds, stringsAsFactors = FALSE)
  full <- merge(grid, comp[, c("sample", "kind", "pct_of_total", "pct_of_degradants")],
                by = c("sample", "kind"), all.x = TRUE)
  full$pct_of_total[is.na(full$pct_of_total)] <- 0
  full$pct_of_degradants[is.na(full$pct_of_degradants) & full$kind != "parent"] <- 0
  full$group <- ps$group[match(full$sample, ps$sample)]
  full$timepoint <- ps$timepoint[match(full$sample, ps$sample)]
  ckey <- paste(full$group, full$timepoint, full$kind, sep = "\r")
  cc <- do.call(rbind, lapply(split(full, ckey), function(d) data.frame(
    group = d$group[1], timepoint = d$timepoint[1], kind = d$kind[1], n = nrow(d),
    mean_pct_of_total = mean(d$pct_of_total), sd_pct_of_total = sd0(d$pct_of_total),
    mean_pct_of_degradants = if (d$kind[1] == "parent") NA_real_ else mean(d$pct_of_degradants, na.rm = TRUE),
    stringsAsFactors = FALSE)))
  cc <- cc[order(cc$group, cc$timepoint, cc$kind), , drop = FALSE]
  tl <- sort(unique(cc$timepoint))
  cc$timepoint_lab <- factor(ifelse(is.na(cc$timepoint), "", format(cc$timepoint, trim = TRUE)),
                             levels = format(tl, trim = TRUE))
  rownames(cc) <- NULL
  if (!has_t) { per$timepoint <- NULL; cc$timepoint <- NULL }
  if (!has_g) { per$group <- NULL; cc$group <- NULL }
  list(design = design, per_condition = per, composition = cc)
}

#' Plot % degradation by group and timepoint
#'
#' @param degradation A [degradation_summary()] result with `sample_meta`.
#' @return A ggplot: mean % degradation (+/- SD) per condition with the
#'   individual replicates as points -- lines over time per group for a
#'   time course, bars per group otherwise.
#' @export
plot_degradation_by_condition <- function(degradation) {
  bc <- degradation_by_condition(degradation)
  if (bc$design == "none") {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = "No Group/Timepoint in the sample table to group by"))
  }
  per <- bc$per_condition
  ps <- degradation$per_sample
  if (!"group" %in% names(per)) per$group <- ""
  if (!"group" %in% names(ps)) ps$group <- ""
  ps$group <- ifelse(is.na(ps$group), "", ps$group)
  per$lo <- per$mean_pct_degradation - ifelse(is.na(per$sd_pct_degradation), 0, per$sd_pct_degradation)
  per$hi <- per$mean_pct_degradation + ifelse(is.na(per$sd_pct_degradation), 0, per$sd_pct_degradation)
  groups <- sort(unique(per$group))
  cols <- stats::setNames(grDevices::colorRampPalette(c("#0279EE", "#FF9400", "#75A025", "#FD9BED", "#E9ED4C"))(length(groups)), groups)
  theme <- ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans")
  if (bc$design == "group") {
    ps <- ps[ps$group %in% groups, , drop = FALSE]
    return(ggplot2::ggplot(per, ggplot2::aes(x = .data$group, y = .data$mean_pct_degradation, fill = .data$group)) +
      ggplot2::geom_col(width = 0.6, alpha = 0.85) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), width = 0.2) +
      ggplot2::geom_jitter(data = ps, ggplot2::aes(x = .data$group, y = .data$pct_degradation),
                           inherit.aes = FALSE, width = 0.08, height = 0, size = 1.8, color = "#1f2430") +
      ggplot2::scale_fill_manual(values = cols, guide = "none") +
      ggplot2::labs(x = "Group", y = "% degradation", title = "% degradation by group",
                    subtitle = "Bar = mean, error bar = SD, points = replicates") + theme)
  }
  ps$timepoint <- suppressWarnings(as.numeric(ps$timepoint))
  ps <- ps[!is.na(ps$timepoint) & ps$group %in% groups, , drop = FALSE]
  ggplot2::ggplot(per, ggplot2::aes(x = .data$timepoint, y = .data$mean_pct_degradation,
                                    color = .data$group, group = .data$group)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), width = 0) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_point(data = ps, ggplot2::aes(x = .data$timepoint, y = .data$pct_degradation,
                                                color = .data$group), inherit.aes = FALSE,
                        size = 1.3, alpha = 0.5) +
    ggplot2::scale_color_manual(values = cols, name = NULL,
                                labels = function(x) ifelse(nzchar(x), x, "all samples")) +
    ggplot2::labs(x = "Timepoint", y = "% degradation", title = "% degradation over time",
                  subtitle = "Line = mean per group, error bar = SD, faint points = replicates") + theme
}

## ---- Degradation relative to an earlier timepoint ----------------------------
#' Parent loss and % degradation relative to a reference timepoint
#'
#' Summarizes [degradation_summary()]'s `per_sample` table by timepoint
#' (and by group/arm, when one is present) and expresses each timepoint
#' against the reference timepoint of the SAME arm and against the
#' immediately preceding timepoint:
#'
#' * `pct_parent_remaining` = 100 x mean parent signal / mean parent signal
#'   at the reference timepoint; `pct_parent_loss` = 100 - that.
#' * `delta_pct_degradation` = mean % degradation minus the reference
#'   timepoint's mean % degradation (percentage points). This one is
#'   ratio-based, so it is robust to injection-to-injection signal drift
#'   that a raw parent-signal ratio is not.
#' * `pct_parent_change_vs_previous` = 100 x (mean parent signal / mean
#'   parent signal at the previous timepoint - 1).
#'
#' @param per_sample `degradation_summary(...)$per_sample`, which carries
#'   `timepoint` (and `group`) when `sample_meta` was supplied.
#' @param reference_timepoint Numeric timepoint to compare against. `NULL`
#'   uses the earliest timepoint of each arm.
#' @return A data.frame, one row per (group, timepoint), or an empty
#'   data.frame when there is no numeric timepoint information.
#' @export
degradation_vs_reference <- function(per_sample, reference_timepoint = NULL) {
  if (is.null(per_sample) || nrow(per_sample) == 0 || !"timepoint" %in% names(per_sample)) {
    return(data.frame())
  }
  d <- per_sample
  d$timepoint <- suppressWarnings(as.numeric(d$timepoint))
  d <- d[!is.na(d$timepoint), , drop = FALSE]
  if (nrow(d) == 0) return(data.frame())
  d$group <- if ("group" %in% names(d)) ifelse(is.na(d$group), "", as.character(d$group)) else ""

  rows <- lapply(split(d, d$group), function(g) {
    tps <- sort(unique(g$timepoint))
    agg <- do.call(rbind, lapply(tps, function(t) {
      s <- g[g$timepoint == t, , drop = FALSE]
      data.frame(group = s$group[1], timepoint = t, n = nrow(s),
                 mean_parent_signal = mean(s$parent_signal, na.rm = TRUE),
                 sd_parent_signal = if (nrow(s) > 1) stats::sd(s$parent_signal, na.rm = TRUE) else NA_real_,
                 mean_pct_degradation = mean(s$pct_degradation, na.rm = TRUE),
                 sd_pct_degradation = if (nrow(s) > 1) stats::sd(s$pct_degradation, na.rm = TRUE) else NA_real_,
                 stringsAsFactors = FALSE)
    }))
    ref_t <- if (is.null(reference_timepoint) || is.na(reference_timepoint)) min(tps) else reference_timepoint
    ref <- agg[agg$timepoint == ref_t, , drop = FALSE]
    ref_parent <- if (nrow(ref) > 0) ref$mean_parent_signal[1] else NA_real_
    ref_degr <- if (nrow(ref) > 0) ref$mean_pct_degradation[1] else NA_real_
    agg$reference_timepoint <- ref_t
    agg$pct_parent_remaining <- if (is.finite(ref_parent) && ref_parent > 0)
      round(100 * agg$mean_parent_signal / ref_parent, 2) else NA_real_
    agg$pct_parent_loss <- round(100 - agg$pct_parent_remaining, 2)
    agg$delta_pct_degradation <- round(agg$mean_pct_degradation - ref_degr, 2)
    prev <- c(NA_real_, utils::head(agg$mean_parent_signal, -1))
    agg$pct_parent_change_vs_previous <- ifelse(is.finite(prev) & prev > 0,
                                                round(100 * (agg$mean_parent_signal / prev - 1), 2),
                                                NA_real_)
    agg
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Plot parent remaining vs. time
#'
#' @param deg_ref Output of [degradation_vs_reference()].
#' @return A ggplot object: % parent remaining (relative to the reference
#'   timepoint) against timepoint, one line per group.
#' @export
plot_degradation_vs_reference <- function(deg_ref) {
  if (is.null(deg_ref) || nrow(deg_ref) == 0 || all(is.na(deg_ref$pct_parent_remaining))) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = "No timepoint information to plot degradation over time"))
  }
  d <- deg_ref
  d$group <- ifelse(nzchar(d$group), d$group, "all samples")
  n_grp <- length(unique(d$group))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$timepoint, y = .data$pct_parent_remaining,
                                  color = .data$group, group = .data$group)) +
    ggplot2::geom_hline(yintercept = 100, linetype = "dotted", color = "#B7B2A7") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_manual(values = grDevices::colorRampPalette(
      c("#0279EE", "#FF9400", "#75A025", "#FD9BED", "#E9ED4C"))(n_grp), name = NULL) +
    ggplot2::labs(x = "Timepoint", y = "% parent remaining",
                  title = "Parent remaining relative to reference timepoint",
                  subtitle = paste0("Reference timepoint: ",
                                    paste(unique(d$reference_timepoint), collapse = ", "))) +
    ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans")
}
