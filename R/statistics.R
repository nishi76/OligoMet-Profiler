# =============================================================================
# statistics.R
# Statistical comparison of confirmed metabolite abundances across samples --
# 2-group, 3+-group, or time-series designs. Independent-samples only (no
# repeated-measures/mixed-effects modeling), keeping this dependency-light:
# only stats:: (t.test, aov, TukeyHSD, lm, p.adjust) and ggplot2, both
# already Imports.
# =============================================================================

## ---- Abundance matrix / long format ----------------------------------------
# Collapses multi-sample MS1 matches (one row per match, from
# match_ms1_batch()/annotate_metabolites_batch()) into a metabolite x sample
# abundance table, taking the most-intense match per (met_id, sample) across
# charge states/adducts/oxidation levels as that metabolite's abundance in
# that sample. Missing (met_id, sample) combinations are NA, not zero.
#
#' Build a metabolite x sample abundance matrix
#'
#' Collapses multi-sample MS1 matches into one row per metabolite, one
#' column per sample, taking the most-intense match per (metabolite,
#' sample) across charge states/adducts/oxidation levels as that
#' metabolite's abundance in that sample.
#'
#' @param batch_matches The `ms1_matches` data.frame from
#'   [match_ms1_batch()]/[annotate_metabolites_batch()]: one row per
#'   (sample, metabolite, charge, adduct) hit, with `met_id`, `met_name`,
#'   `kind`, `sample`, and a signal column (`intensity` by default).
#' @param signal_col Which column to use as signal (e.g. `"intensity"`,
#'   `"area"`, or a blank-corrected column such as `"intensity_bcorr"` --
#'   see [apply_blank_correction()]). Defaults to `"intensity"`, preserving
#'   this function's original behavior.
#' @return A data.frame with one row per metabolite (`met_id`, `met_name`,
#'   `kind`) and one column per sample (named after that sample) holding
#'   its max signal for that metabolite. A (metabolite, sample)
#'   combination with no match is `NA`, not zero. Empty data.frame if
#'   `batch_matches` is `NULL`/empty, or if `signal_col` isn't present at
#'   all (e.g. `"area"` on single-file-mode features, which the R-native
#'   reader doesn't compute AUC for).
#' @seealso [abundance_long()] to reshape this into long format for
#'   [compare_two_groups()]/[compare_multi_groups()]/[compare_time_series()].
#' @export
build_abundance_matrix <- function(batch_matches, signal_col = "intensity") {
  if (is.null(batch_matches) || nrow(batch_matches) == 0) return(data.frame())
  if (!signal_col %in% names(batch_matches) || all(is.na(batch_matches[[signal_col]]))) return(data.frame())
  met_info <- unique(batch_matches[, c("met_id", "met_name", "kind")])
  met_info <- met_info[order(met_info$met_id), ]
  samples <- unique(batch_matches$sample)
  out <- met_info
  for (s in samples) {
    sub <- batch_matches[batch_matches$sample == s & !is.na(batch_matches[[signal_col]]), ]
    if (nrow(sub) == 0) { out[[s]] <- NA_real_; next }
    best <- stats::aggregate(stats::as.formula(paste(signal_col, "~ met_id")), data = sub, FUN = max)
    out[[s]] <- best[[signal_col]][match(out$met_id, best$met_id)]
  }
  out
}

# Same as build_abundance_matrix(), but on peak area (trapezoidal AUC) --
# the batch/Python ROI pipeline's alternative to max intensity (see the
# `area` column threaded through match_ms1()/match_ms1_batch() in
# R/ms_matching.R and R/batch_ms_processing.R). Kept as a thin wrapper so
# existing callers referencing it by name keep working unchanged; callers
# (degradation_summary() in R/degradation.R) fall back to
# build_abundance_matrix()'s intensity when this returns empty.
build_abundance_matrix_area <- function(batch_matches) {
  build_abundance_matrix(batch_matches, signal_col = "area")
}

# Long-format companion: sample_meta is data.frame(sample, group) and/or
# data.frame(sample, timepoint) -- whichever columns are present besides
# `sample` are carried through untouched.
#
#' Reshape an abundance matrix to long format
#'
#' Companion to [build_abundance_matrix()]/[build_kind_abundance_matrix()]:
#' pivots the wide metabolite x sample matrix to one row per (metabolite,
#' sample) with sample metadata attached, the shape
#' [compare_two_groups()]/[compare_multi_groups()]/[compare_time_series()]
#' and the quantification functions expect.
#'
#' @param abundance_matrix Output of [build_abundance_matrix()] or
#'   [build_kind_abundance_matrix()].
#' @param sample_meta A data.frame with a `sample` column plus whichever
#'   of `group`/`timepoint`/`sample_type`/`concentration` are relevant --
#'   every non-`sample` column is carried through untouched.
#' @return A data.frame: one row per (metabolite, sample), with `met_id`,
#'   `met_name`, `kind`, `sample`, `intensity`, and every `sample_meta`
#'   column besides `sample`. Empty data.frame if `abundance_matrix` is
#'   empty, or if a sample has no matching row in `sample_meta`.
#' @export
abundance_long <- function(abundance_matrix, sample_meta) {
  if (is.null(abundance_matrix) || nrow(abundance_matrix) == 0) return(data.frame())
  meta_cols <- setdiff(names(sample_meta), "sample")
  sample_cols <- setdiff(names(abundance_matrix), c("met_id", "met_name", "kind"))
  rows <- lapply(sample_cols, function(s) {
    meta <- sample_meta[sample_meta$sample == s, , drop = FALSE]
    if (nrow(meta) == 0) return(NULL)
    out <- abundance_matrix[, c("met_id", "met_name", "kind")]
    out$sample <- s
    out$intensity <- abundance_matrix[[s]]
    for (col in meta_cols) out[[col]] <- meta[[col]][1]
    out
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

## ---- Kind-level (composition-class) abundance --------------------------------
# Sums signal (intensity, or area if you have it -- see build_abundance_matrix_area())
# across all metabolites sharing a `kind` value (parent/exo_3p/exo_5p/
# endo_5frag/endo_3frag/endo_5frag_p/endo_3frag_p -- the same rollup
# degradation_summary() does
# internally, see R/degradation.R), per sample. Exposed here as a standalone
# long table so it can be run through the EXISTING compare_two_groups()/
# compare_multi_groups()/compare_time_series() completely unmodified, by
# presenting `kind` values as if they were `met_id` -- those three
# functions are structurally generic over "whatever the grouping-ID column
# is called," so this costs zero new statistical code. The resulting
# met_id/met_name columns will literally hold "exo_3p"/"endo_5frag" etc.;
# rename for display at the call site (Shiny/Excel), not here.
build_kind_abundance_matrix <- function(batch_matches, signal_col = "intensity") {
  if (is.null(batch_matches) || nrow(batch_matches) == 0) return(data.frame())
  if (!signal_col %in% names(batch_matches)) return(data.frame())
  m <- batch_matches[!is.na(batch_matches[[signal_col]]), ]
  if (nrow(m) == 0) return(data.frame())
  # Collapse to one row per (sample, met_id) first -- max across charge
  # states/adducts/oxidation levels of the SAME metabolite, the same rule
  # build_abundance_matrix()/degradation_summary() use -- so a metabolite
  # matched at more than one charge state doesn't get double-counted when
  # summed into its kind's total below.
  best <- stats::aggregate(stats::as.formula(paste(signal_col, "~ sample + met_id")),
                            data = m, FUN = max)
  best <- merge(best, unique(m[, c("met_id", "kind")]), by = "met_id")
  kinds <- sort(unique(best$kind))
  samples <- unique(best$sample)
  out <- data.frame(met_id = kinds, met_name = kinds, kind = kinds, stringsAsFactors = FALSE)
  for (s in samples) {
    sub <- best[best$sample == s, ]
    agg <- stats::aggregate(stats::as.formula(paste(signal_col, "~ kind")), data = sub, FUN = sum)
    out[[s]] <- agg[[signal_col]][match(out$kind, agg$kind)]
  }
  out
}

# Long-format companion, identical shape/contract to abundance_long() -- a
# verbatim pass-through, since abundance_long() only ever references
# met_id/met_name/kind/sample columns generically.
kind_abundance_long <- function(kind_abundance_matrix, sample_meta) {
  abundance_long(kind_abundance_matrix, sample_meta)
}

## ---- Two-group comparison ---------------------------------------------------
#' Compare metabolite abundance between two groups
#'
#' Independent-samples Welch's t-test per metabolite, with
#' Benjamini-Hochberg FDR correction across metabolites.
#'
#' @param abundance_long Long-format abundance data.frame from
#'   [abundance_long()], with a `group` column.
#' @param group_a,group_b The two `group` values to compare. `log2fc` is
#'   `log2(mean(group_b) / mean(group_a))`, so `group_a` is the reference.
#' @param min_n Minimum replicates required per group to run the test
#'   (clamped to at least 2, since a t-test needs to estimate variance).
#'   A metabolite with fewer replicates gets `note = "insufficient
#'   replicates"` and `NA` statistics instead of being dropped.
#' @param p_adjust_method Method passed to [stats::p.adjust()] for the
#'   across-metabolite correction (`"BH"` by default).
#' @return A data.frame with one row per metabolite: `met_id`, `met_name`,
#'   `mean_a`, `mean_b`, `log2fc`, `t_stat`, `p_value`, `p_adj` (adjusted per
#'   `p_adjust_method`), `n_a`, `n_b`, `note`.
#' @seealso [plot_volcano()] to visualize this result.
#' @export
compare_two_groups <- function(abundance_long, group_a, group_b, min_n = 2,
                                p_adjust_method = "BH") {
  if (is.null(abundance_long) || nrow(abundance_long) == 0) return(data.frame())
  min_n <- max(min_n, 2)  # stats::t.test() cannot estimate variance from a single observation
  mets <- unique(abundance_long[, c("met_id", "met_name")])
  rows <- lapply(seq_len(nrow(mets)), function(i) {
    mid <- mets$met_id[i]
    sub <- abundance_long[abundance_long$met_id == mid, ]
    a <- sub$intensity[sub$group == group_a]; a <- a[!is.na(a)]
    b <- sub$intensity[sub$group == group_b]; b <- b[!is.na(b)]
    if (length(a) < min_n || length(b) < min_n) {
      return(data.frame(
        met_id = mid, met_name = mets$met_name[i],
        mean_a = if (length(a) > 0) mean(a) else NA_real_,
        mean_b = if (length(b) > 0) mean(b) else NA_real_,
        log2fc = NA_real_, t_stat = NA_real_, p_value = NA_real_,
        n_a = length(a), n_b = length(b), note = "insufficient replicates",
        stringsAsFactors = FALSE))
    }
    tt <- stats::t.test(a, b, var.equal = FALSE)
    data.frame(met_id = mid, met_name = mets$met_name[i],
               mean_a = mean(a), mean_b = mean(b), log2fc = log2(mean(b) / mean(a)),
               t_stat = unname(tt$statistic), p_value = tt$p.value,
               n_a = length(a), n_b = length(b), note = "", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust_method)
  out
}

## ---- Three-or-more-group comparison -----------------------------------------
#' Compare metabolite abundance across 3+ groups
#'
#' One-way ANOVA per metabolite (omnibus test across all groups), with
#' Tukey HSD post-hoc pairwise contrasts, and Benjamini-Hochberg FDR
#' correction on the omnibus p-values across metabolites.
#'
#' @param abundance_long Long-format abundance data.frame from
#'   [abundance_long()], with a `group` column.
#' @param groups Which `group` values to include, and in what order for
#'   the ANOVA factor levels. Defaults to every distinct value in
#'   `abundance_long$group`.
#' @param min_n Minimum replicates required in EVERY group to run the
#'   test for a metabolite (clamped to at least 2). Below that, the
#'   metabolite gets `note = "insufficient replicates"` and `NA`
#'   statistics instead of being dropped.
#' @param p_adjust_method Method passed to [stats::p.adjust()] for the
#'   across-metabolite correction on the omnibus p-values (`"BH"` by
#'   default). The Tukey HSD post-hoc `p_adj` is unaffected -- that's
#'   already family-wise adjusted by [stats::TukeyHSD()] itself.
#' @return `list(omnibus, posthoc)`: `omnibus` has one row per metabolite
#'   (`met_id`, `met_name`, `f_stat`, `p_value`, `p_adj`, `note`);
#'   `posthoc` has one row per (metabolite, pairwise contrast) with the
#'   Tukey HSD `diff`/`lwr`/`upr`/`p_adj`.
#' @export
compare_multi_groups <- function(abundance_long, groups = NULL, min_n = 2,
                                  p_adjust_method = "BH") {
  empty <- list(omnibus = data.frame(), posthoc = data.frame())
  if (is.null(abundance_long) || nrow(abundance_long) == 0) return(empty)
  if (is.null(groups)) groups <- unique(abundance_long$group)
  mets <- unique(abundance_long[, c("met_id", "met_name")])

  omnibus_rows <- list(); posthoc_rows <- list()
  for (i in seq_len(nrow(mets))) {
    mid <- mets$met_id[i]
    sub <- abundance_long[abundance_long$met_id == mid & abundance_long$group %in% groups, ]
    sub <- sub[!is.na(sub$intensity), ]
    counts <- table(factor(sub$group, levels = groups))

    if (length(groups) < 2 || any(counts < min_n)) {
      omnibus_rows[[length(omnibus_rows) + 1]] <- data.frame(
        met_id = mid, met_name = mets$met_name[i], f_stat = NA_real_,
        p_value = NA_real_, note = "insufficient replicates", stringsAsFactors = FALSE)
      next
    }

    sub$group <- factor(sub$group, levels = groups)
    fit <- stats::aov(intensity ~ group, data = sub)
    at <- summary(fit)[[1]]
    omnibus_rows[[length(omnibus_rows) + 1]] <- data.frame(
      met_id = mid, met_name = mets$met_name[i],
      f_stat = at["group", "F value"], p_value = at["group", "Pr(>F)"],
      note = "", stringsAsFactors = FALSE)

    tk <- stats::TukeyHSD(fit)
    tk_df <- as.data.frame(tk$group)
    posthoc_rows[[length(posthoc_rows) + 1]] <- data.frame(
      met_id = mid, met_name = mets$met_name[i], contrast = rownames(tk_df),
      diff = tk_df$diff, lwr = tk_df$lwr, upr = tk_df$upr,
      p_adj = tk_df[["p adj"]], stringsAsFactors = FALSE)
  }

  omnibus <- do.call(rbind, omnibus_rows)
  omnibus$p_adj <- stats::p.adjust(omnibus$p_value, method = p_adjust_method)
  posthoc <- if (length(posthoc_rows) > 0) do.call(rbind, posthoc_rows) else data.frame()
  list(omnibus = omnibus, posthoc = posthoc)
}

## ---- Time series (linear trend) ---------------------------------------------
# Independent-samples trend test: per metabolite, a linear regression of
# intensity on timepoint. If replicates exist per timepoint but represent
# the SAME tracked biological replicate across time (a paired/repeated-
# measures design), this OLS slope test is still a valid trend test but
# does not correct for within-subject correlation -- a documented
# simplification; mixed-effects modeling would need a new R dependency
# (nlme/lme4) not otherwise required by this package.
#' Test for a linear trend across timepoints
#'
#' Per metabolite, an independent-samples linear regression of abundance
#' on timepoint. If replicates exist per timepoint but represent the SAME
#' tracked biological replicate across time (a paired/repeated-measures
#' design), this OLS slope test is still valid as a trend test but does
#' NOT correct for within-subject correlation -- a documented
#' simplification, since mixed-effects modeling (nlme/lme4) isn't
#' otherwise a dependency of this package.
#'
#' @param abundance_long Long-format abundance data.frame from
#'   [abundance_long()], with a `timepoint` (or `time_var`) column.
#' @param time_var Name of the timepoint column; coerced to numeric via
#'   `as.numeric()` (non-numeric/`NA` values are dropped).
#' @param min_n Minimum total observations, across at least 2 distinct
#'   timepoints, required to fit the regression for a metabolite. Below
#'   that, the metabolite gets `note = "insufficient data"` and `NA`
#'   statistics instead of being dropped.
#' @param p_adjust_method Method passed to [stats::p.adjust()] for the
#'   across-metabolite correction (`"BH"` by default).
#' @return A data.frame with one row per metabolite: `met_id`, `met_name`,
#'   `slope`, `slope_se`, `p_value`, `p_adj` (adjusted per
#'   `p_adjust_method`), `r_squared`, `n_timepoints`, `n_total`, `note`.
#' @seealso [plot_trend()] to visualize this result for one metabolite.
#' @export
compare_time_series <- function(abundance_long, time_var = "timepoint", min_n = 2,
                                 p_adjust_method = "BH") {
  if (is.null(abundance_long) || nrow(abundance_long) == 0) return(data.frame())
  if (!time_var %in% names(abundance_long)) {
    stop("abundance_long has no '", time_var, "' column")
  }
  mets <- unique(abundance_long[, c("met_id", "met_name")])
  rows <- lapply(seq_len(nrow(mets)), function(i) {
    mid <- mets$met_id[i]
    sub <- abundance_long[abundance_long$met_id == mid, ]
    sub$.time <- suppressWarnings(as.numeric(sub[[time_var]]))
    sub <- sub[!is.na(sub$intensity) & !is.na(sub$.time), ]
    n_tp <- length(unique(sub$.time))

    if (nrow(sub) < min_n || n_tp < 2) {
      return(data.frame(met_id = mid, met_name = mets$met_name[i],
                         slope = NA_real_, slope_se = NA_real_, p_value = NA_real_,
                         r_squared = NA_real_, n_timepoints = n_tp, n_total = nrow(sub),
                         note = "insufficient data", stringsAsFactors = FALSE))
    }

    fit <- stats::lm(intensity ~ .time, data = sub)
    s <- summary(fit)
    coefs <- s$coefficients
    data.frame(met_id = mid, met_name = mets$met_name[i],
               slope = coefs[2, 1], slope_se = coefs[2, 2], p_value = coefs[2, 4],
               r_squared = s$r.squared, n_timepoints = n_tp, n_total = nrow(sub),
               note = "", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust_method)
  out
}

## ---- Plots (ggplot2, matching R/build_report.R conventions) ----------------
# All return a ggplot object; callers decide whether to render (Shiny) or
# save (ggplot2::ggsave(..., width=7, height=4, dpi=150, bg="white"), the
# same call build_report.R already uses for its other plots).

plot_volcano <- function(two_group_result, fc_thresh = 1, p_thresh = 0.05, label_top_n = 10) {
  df <- two_group_result[!is.na(two_group_result$p_adj), ]
  df$sig <- abs(df$log2fc) >= fc_thresh & df$p_adj <= p_thresh
  df$neglog10p <- -log10(pmax(df$p_adj, .Machine$double.xmin))
  top <- df[order(df$p_adj), ][seq_len(min(label_top_n, nrow(df))), ]

  ggplot2::ggplot(df, ggplot2::aes(x = .data$log2fc, y = .data$neglog10p, color = .data$sig)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#0279EE", `FALSE` = "#B7B2A7"), guide = "none") +
    ggplot2::geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "#75A025") +
    ggplot2::geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "#75A025") +
    ggplot2::geom_text(data = top, ggplot2::aes(label = .data$met_name),
                        size = 3, color = "#333333", vjust = -0.6, check_overlap = TRUE) +
    ggplot2::labs(x = "log2 fold change", y = "-log10(adjusted p)", title = "Group comparison") +
    ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans")
}

plot_trend <- function(time_series_long, met_id) {
  sub <- time_series_long[time_series_long$met_id == met_id, ]
  ggplot2::ggplot(sub, ggplot2::aes(x = .data$timepoint, y = .data$intensity)) +
    ggplot2::geom_jitter(width = 0, height = 0, size = 2, color = "#0279EE", alpha = 0.8) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#FF9400", fill = "#FF9400", alpha = 0.15) +
    ggplot2::labs(x = "Timepoint", y = "Intensity",
                  title = if (nrow(sub) > 0) sub$met_name[1] else met_id) +
    ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans")
}

plot_group_boxplot <- function(abundance_long, met_id) {
  sub <- abundance_long[abundance_long$met_id == met_id, ]
  ggplot2::ggplot(sub, ggplot2::aes(x = .data$group, y = .data$intensity, fill = .data$group)) +
    ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, size = 2, color = "#333333") +
    ggplot2::scale_fill_manual(values = grDevices::colorRampPalette(
      c("#0279EE", "#FF9400", "#75A025", "#FD9BED", "#E9ED4C"))(length(unique(sub$group))), guide = "none") +
    ggplot2::labs(x = NULL, y = "Intensity",
                  title = if (nrow(sub) > 0) sub$met_name[1] else met_id) +
    ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans")
}

## ---- Multi-metabolite time-course trend (Time Course tab) ------------------
# Several metabolites overlaid on raw intensity would be dominated by
# whichever one happens to have the largest signal -- unlike plot_trend()'s
# single metabolite (its own y-axis, an lm trend line), each metabolite here
# is normalized to its own mean signal at `reference_timepoint` by default,
# the same baseline quantify_relative(mode = "time_series") already uses,
# so metabolites spanning orders of magnitude in raw signal are comparable
# on one shared scale.
#' Summarize a multi-metabolite time-course trend
#'
#' One row per (metabolite, timepoint): mean, SD, SEM, and replicate count,
#' optionally normalized to each metabolite's own mean signal at a
#' reference timepoint.
#'
#' @param time_series_long Long-format abundance data.frame (see
#'   [abundance_long()]) with a `timepoint` column, already numeric (the
#'   caller is expected to have converted/filtered it, the same way
#'   [compare_time_series()]'s own callers do).
#' @param met_ids Which metabolites to summarize; `NULL` (the default)
#'   summarizes every metabolite present.
#' @param normalize `TRUE` (the default): divide each metabolite's signal
#'   by its own mean at `reference_timepoint`. `FALSE` keeps raw
#'   intensity.
#' @param reference_timepoint The numeric `timepoint` value to normalize
#'   against when `normalize = TRUE`. `NULL` (the default) uses the
#'   earliest timepoint present -- see [quantify_relative()]'s own
#'   `reference_timepoint` for why that isn't always "pre-dose".
#' @return A data.frame: `met_id`, `met_name`, `timepoint`, `n`,
#'   `mean_value`, `sd`, `sem`. Empty if there's nothing to summarize.
#' @seealso [plot_multi_trend()], [quantify_relative()] for the same
#'   baseline definition used elsewhere.
#' @export
multi_trend_summary <- function(time_series_long, met_ids = NULL, normalize = TRUE,
                                 reference_timepoint = NULL) {
  df <- time_series_long
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  if (!is.null(met_ids)) df <- df[df$met_id %in% met_ids, ]
  df <- df[!is.na(df$intensity) & !is.na(df$timepoint), ]
  if (nrow(df) == 0) return(data.frame())

  if (normalize) {
    ref_time <- if (is.null(reference_timepoint)) min(df$timepoint) else reference_timepoint
    baseline <- stats::aggregate(intensity ~ met_id,
                                  data = df[df$timepoint == ref_time, , drop = FALSE], FUN = mean)
    df$.baseline <- baseline$intensity[match(df$met_id, baseline$met_id)]
    df <- df[!is.na(df$.baseline) & df$.baseline > 0, ]
    if (nrow(df) == 0) return(data.frame())
    df$.value <- df$intensity / df$.baseline
  } else {
    df$.value <- df$intensity
  }

  rows <- lapply(split(df, list(df$met_id, df$timepoint), drop = TRUE), function(g) {
    v <- g$.value[!is.na(g$.value)]
    data.frame(met_id = g$met_id[1], met_name = g$met_name[1], timepoint = g$timepoint[1],
               n = length(v), mean_value = mean(v),
               sd = if (length(v) > 1) stats::sd(v) else NA_real_,
               sem = if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out[order(out$met_id, out$timepoint), ]
}

#' Plot a multi-metabolite time-course trend
#'
#' One line per metabolite (mean +/- SEM per timepoint), from a
#' [multi_trend_summary()] result -- for comparing several metabolites'
#' kinetics on one shared axis, which [plot_trend()]'s single-metabolite,
#' own-scale view can't do.
#'
#' @param trend_summary A [multi_trend_summary()] result.
#' @param normalize Whether `trend_summary` was built with
#'   `normalize = TRUE` -- only affects the y-axis label.
#' @param reference_timepoint The `reference_timepoint` (if any) that
#'   `trend_summary` was normalized against -- only affects the y-axis
#'   label (`NULL` labels it "earliest timepoint", matching
#'   [multi_trend_summary()]'s own default).
#' @return A ggplot object.
#' @export
plot_multi_trend <- function(trend_summary, normalize = TRUE, reference_timepoint = NULL) {
  if (is.null(trend_summary) || nrow(trend_summary) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = "No data for the selected metabolite(s)"))
  }
  n_met <- length(unique(trend_summary$met_name))
  y_lab <- if (!normalize) "Intensity" else if (is.null(reference_timepoint)) {
    "Relative signal (vs. earliest timepoint)"
  } else {
    paste0("Relative signal (vs. timepoint ", reference_timepoint, ")")
  }
  ggplot2::ggplot(trend_summary, ggplot2::aes(x = .data$timepoint, y = .data$mean_value,
                                               color = .data$met_name, group = .data$met_name)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$mean_value - .data$sem, ymax = .data$mean_value + .data$sem),
                            width = 0, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = grDevices::colorRampPalette(
      c("#0279EE", "#FF9400", "#75A025", "#FD9BED", "#E9ED4C"))(n_met), name = NULL) +
    ggplot2::labs(x = "Timepoint", y = y_lab, title = "Multi-metabolite time course") +
    ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans")
}

## =============================================================================
## Quantification: absolute (calibration curve) + relative (fold-change over
## a baseline). A user-selected subset of metabolites gets absolute
## quantification; everything else gets relative quantification -- see
## quantify_metabolites() at the bottom, which splits met_ids and dispatches
## to the two halves below. Reuses .best_signal_per_met() from
## R/degradation.R (same package, no export needed) for the "max across
## charge states/adducts/oxidation levels" collapse every other
## per-(sample, met_id) signal calculation in this codebase already uses.
## =============================================================================

.auto_signal_col <- function(batch_matches) {
  if ("area" %in% names(batch_matches) && any(!is.na(batch_matches$area))) "area" else "intensity"
}

.calibration_weights <- function(conc, scheme) {
  switch(scheme,
    "1/x2" = 1 / conc^2,
    "1/x"  = 1 / conc,
    "none" = rep(1, length(conc)),
    stop("Unknown weighting scheme: '", scheme, "' (expected one of 1/x2, 1/x, none)"))
}

## ---- Calibration curve for one metabolite -----------------------------------
# ONLY sample_meta rows with sample_type == "standard" and a positive numeric
# concentration build the curve. Every other row (unknown study samples,
# quality_control, reagent_blank, matrix_blank) is deliberately excluded
# here -- quantify_absolute() below back-calculates concentration FOR them
# FROM this fitted curve, they never feed into fitting it. That split (curve
# built from standards only, validated against QC) is the actual point of
# having a `sample_type` vocabulary at all, not just a labeling convenience.
#
# The `weighting` scheme (1/x^2, 1/x, or unweighted OLS) is a real
# methodology choice with no universally-correct default for LC-MS
# bioanalytical quantification -- left as an explicit, required argument
# rather than silently defaulting, so a caller (the Shiny UI) has to pick.
#' Fit a calibration curve for one metabolite
#'
#' Builds a (weighted) linear regression of signal vs. concentration from
#' calibration standards, for back-calculating unknown concentrations
#' (see [quantify_absolute()]). ONLY `sample_meta` rows with
#' `sample_type == "standard"` and a positive numeric `concentration`
#' contribute to the fit -- quality control and unknown study samples are
#' back-calculated FROM this curve elsewhere, never used to build it.
#'
#' @param batch_matches The `ms1_matches` data.frame (see
#'   [match_ms1_batch()]/[annotate_metabolites_batch()]).
#' @param met_id The single metabolite ID to fit a curve for.
#' @param sample_meta A data.frame with `sample`, `sample_type`, and
#'   `concentration` columns (as produced by the Shiny app's batch sample
#'   metadata table, or hand-built to match).
#' @param weighting Regression weighting scheme: `"1/x2"` (weight =
#'   `1/concentration^2`, the common default for LC-MS/MS bioanalytical
#'   quantification), `"1/x"`, or `"none"` (unweighted OLS). No universally
#'   correct default exists across labs/instruments/SOPs, so this has no
#'   default of its own -- callers must pick.
#' @param signal_col Which column to use as signal. `NULL` (the default)
#'   auto-detects: `"area"` if present and not all-`NA`, else
#'   `"intensity"`.
#' @param min_points Minimum number of standard points (after dropping
#'   zero/non-positive concentrations, which can't be weighted by 1/x or
#'   1/x^2 and carry no line-fitting information anyway) required to fit
#'   the curve. Below this, an empty curve is returned with a `note`
#'   explaining why, rather than erroring.
#' @return A list: `met_id`, `signal_col`, `weighting`, `model` (the
#'   fitted `lm` object, or `NULL` if the curve couldn't be built),
#'   `intercept`, `slope`, `r_squared`, `n_points`, `conc_range` (the
#'   observed standard concentration range -- see the `extrapolated` flag
#'   in [quantify_absolute()]), `points` (the standards used, with their
#'   own back-calculated concentration and `percent_re`), and `note`
#'   (empty string on success, else why the curve is empty).
#' @seealso [quantify_absolute()], [quantify_metabolites()]
#' @export
fit_calibration_curve <- function(batch_matches, met_id, sample_meta,
                                   weighting = c("1/x2", "1/x", "none"),
                                   signal_col = NULL, min_points = 3) {
  weighting <- match.arg(weighting)
  if (is.null(signal_col)) signal_col <- .auto_signal_col(batch_matches)

  empty <- list(met_id = met_id, signal_col = signal_col, weighting = weighting,
                model = NULL, intercept = NA_real_, slope = NA_real_, r_squared = NA_real_,
                n_points = 0, conc_range = c(NA_real_, NA_real_),
                points = data.frame(), note = "")

  if (!all(c("sample_type", "concentration") %in% names(sample_meta))) {
    empty$note <- "sample_meta has no sample_type/concentration column"
    return(empty)
  }
  std_meta <- sample_meta[!is.na(sample_meta$sample_type) & sample_meta$sample_type == "standard" &
                             !is.na(sample_meta$concentration) & nzchar(sample_meta$concentration), ,
                           drop = FALSE]
  if (nrow(std_meta) == 0) {
    empty$note <- "no standard-type samples with a concentration value"
    return(empty)
  }
  std_meta$concentration <- suppressWarnings(as.numeric(std_meta$concentration))
  # A zero-concentration calibrator (double blank) is legitimate for
  # establishing the LLOQ but breaks 1/x and 1/x^2 weighting (division by
  # zero) and contributes no information to a line's slope anyway -- drop
  # it from the fit itself rather than let it silently produce Inf weights.
  std_meta <- std_meta[!is.na(std_meta$concentration) & std_meta$concentration > 0, ]
  if (nrow(std_meta) == 0) {
    empty$note <- "no standards with a positive numeric concentration"
    return(empty)
  }

  met_matches <- batch_matches[batch_matches$met_id == met_id, ]
  if (nrow(met_matches) == 0 || !signal_col %in% names(met_matches)) {
    empty$note <- paste0("no matches for ", met_id)
    return(empty)
  }
  best <- .best_signal_per_met(met_matches[!is.na(met_matches[[signal_col]]), ], signal_col)
  pts <- merge(std_meta[, c("sample", "concentration")], best[, c("sample", signal_col)], by = "sample")
  names(pts)[names(pts) == signal_col] <- "signal"
  pts <- pts[!is.na(pts$signal), ]

  if (nrow(pts) < max(min_points, 2)) {
    empty$note <- sprintf("only %d standard point(s) with signal for %s (need >= %d)",
                           nrow(pts), met_id, min_points)
    empty$points <- pts
    return(empty)
  }

  w <- .calibration_weights(pts$concentration, weighting)
  fit <- stats::lm(signal ~ concentration, data = pts, weights = w)
  s <- summary(fit)
  intercept <- unname(stats::coef(fit)[1]); slope <- unname(stats::coef(fit)[2])
  pts$back_calc_concentration <- (pts$signal - intercept) / slope
  pts$percent_re <- 100 * (pts$back_calc_concentration - pts$concentration) / pts$concentration

  list(met_id = met_id, signal_col = signal_col, weighting = weighting, model = fit,
       intercept = intercept, slope = slope, r_squared = s$r.squared,
       n_points = nrow(pts), conc_range = range(pts$concentration),
       points = pts, note = "")
}

## ---- Calibration curve plot -------------------------------------------------
# Standards (fit_calibration_curve()'s own `points`) as the scatter the
# line is actually fit to, plus dotted vlines at the fitted conc_range --
# the same boundary quantify_absolute()'s `extrapolated` flag is computed
# from, so a point past the line visually IS an extrapolated one. QC/
# unknown samples are optional: quantify_absolute()'s `quant` table
# back-calculates concentration for every sample with signal (standards
# included), not just the ones used to build the curve, so passing that
# subset in lets a QC or unknown sample's own (back-calculated
# concentration, signal) show where it actually falls against the line --
# still informative relative to the extrapolation boundary even though,
# by construction, back-calculated points always sit exactly on the line.
#' Plot a calibration curve
#'
#' Scatter of standard signal vs. concentration with the fitted regression
#' line, dotted vertical lines marking the calibrated (non-extrapolated)
#' concentration range, and an r-squared/weighting/n annotation.
#' Optionally overlays QC/unknown samples' own back-calculated
#' concentration in a second color.
#'
#' @param curve_result A [fit_calibration_curve()] result.
#' @param quant_points Optional: the subset of [quantify_absolute()]'s
#'   `quant` data.frame for this SAME `met_id` (`sample_type`,
#'   `concentration_calc`, `signal` columns) -- standard rows are dropped
#'   automatically (already shown from `curve_result$points`), so this can
#'   just be every row for the metabolite.
#' @return A ggplot object. If `curve_result$model` is `NULL` (curve
#'   couldn't be fit), returns an empty plot titled with
#'   `curve_result$note` instead of erroring.
#' @seealso [fit_calibration_curve()], [quantify_absolute()].
#' @export
plot_calibration_curve <- function(curve_result, quant_points = NULL) {
  # Checks note/intercept/slope, NOT $model -- a curve restored from a
  # saved Analysis State has $model set to NULL (stripped for
  # serialization, see .strip_for_analysis_state() in app.R) even though
  # everything actually needed to plot (intercept, slope, r_squared,
  # points) survives the round-trip untouched.
  if (is.null(curve_result) || nzchar(curve_result$note %||% "") ||
      is.null(curve_result$intercept) || is.na(curve_result$intercept)) {
    msg <- if (!is.null(curve_result) && nzchar(curve_result$note %||% "")) curve_result$note else "No calibration curve to plot"
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = msg))
  }
  pts <- curve_result$points
  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = curve_result$conc_range, linetype = "dotted", color = "#B7B2A7") +
    ggplot2::geom_abline(intercept = curve_result$intercept, slope = curve_result$slope,
                          color = "#75A025", linewidth = 0.8) +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = .data$concentration, y = .data$signal),
                         color = "#0279EE", size = 2.5, alpha = 0.85)

  if (!is.null(quant_points) && nrow(quant_points) > 0) {
    qp <- quant_points[!is.na(quant_points$concentration_calc) & !is.na(quant_points$signal), ]
    is_std <- !is.na(qp$sample_type) & qp$sample_type == "standard"
    qp <- qp[!is_std, ]
    if (nrow(qp) > 0) {
      qp$.type <- ifelse(!is.na(qp$sample_type) & qp$sample_type == "quality_control", "QC", "unknown")
      p <- p + ggplot2::geom_point(
        data = qp, ggplot2::aes(x = .data$concentration_calc, y = .data$signal, shape = .data$.type),
        color = "#FD9BED", size = 2.5, alpha = 0.85) +
        ggplot2::scale_shape_manual(values = c(QC = 17, unknown = 15), name = NULL)
    }
  }

  p + ggplot2::labs(
    x = "Concentration", y = "Signal",
    title = paste0("Calibration curve -- ", curve_result$met_id),
    subtitle = sprintf("R^2 = %.4f   n = %d standards   weighting = %s   signal = %s",
                        curve_result$r_squared, curve_result$n_points,
                        curve_result$weighting, curve_result$signal_col)) +
    ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans")
}

## ---- Absolute quantification (a set of metabolites, all study samples) -----
# For each metabolite in met_ids, fits its own calibration curve (see
# fit_calibration_curve() above) and back-calculates concentration for
# EVERY sample that has signal for it -- standards and QC included, not
# just unknowns, so a standard's or QC's back-calculated value can be
# checked against its own nominal concentration (percent_re) as an
# accuracy/data-quality signal, the same way a bioanalytical run's
# calibration report would. `extrapolated` flags a back-calculated
# concentration outside the curve's own observed standard range --
# extrapolating past the calibrated range is unreliable, not just
# borderline, so this is worth carrying as its own column rather than
# leaving it to look identical to an interpolated result.
#' Absolute-quantify a set of metabolites from their own calibration curves
#'
#' For each metabolite in `met_ids`, fits its own calibration curve (see
#' [fit_calibration_curve()]) and back-calculates concentration for EVERY
#' sample with signal for it -- standards and quality-control samples
#' included, not just unknowns, so a standard's/QC's back-calculated value
#' can be checked against its own nominal concentration as an
#' accuracy/data-quality signal.
#'
#' @param batch_matches The `ms1_matches` data.frame (see
#'   [match_ms1_batch()]/[annotate_metabolites_batch()]).
#' @param met_ids Metabolite IDs to absolute-quantify.
#' @param sample_meta A data.frame with `sample`, `sample_type`, and
#'   `concentration` columns.
#' @param weighting Calibration curve weighting -- see
#'   [fit_calibration_curve()].
#' @param signal_col Which column to use as signal; `NULL` auto-detects
#'   (see [fit_calibration_curve()]).
#' @param min_points Minimum calibration standard points required per
#'   metabolite -- see [fit_calibration_curve()].
#' @return `list(quant, curves)`: `quant` is a data.frame with one row
#'   per (metabolite, sample) -- `met_id`, `met_name`, `sample`,
#'   `sample_type`, `signal`, `concentration_calc` (back-calculated from
#'   the curve), `extrapolated` (`TRUE` if `concentration_calc` falls
#'   outside the curve's own observed standard range -- unreliable, not
#'   just borderline), `nominal_concentration`, `percent_re` (only
#'   populated for standard/QC rows with a nominal concentration), and
#'   the curve's own `curve_r_squared`/`curve_n_points`/`curve_weighting`.
#'   `curves` is a named list (by `met_id`) of each metabolite's
#'   [fit_calibration_curve()] result, e.g. for plotting or QC review.
#'   A metabolite whose curve couldn't be fit contributes no rows to
#'   `quant` but its (empty) curve is still in `curves`.
#' @seealso [quantify_relative()] for metabolites NOT selected for
#'   absolute quantification, [quantify_metabolites()] to dispatch both
#'   in one call.
#' @export
quantify_absolute <- function(batch_matches, met_ids, sample_meta,
                               weighting = c("1/x2", "1/x", "none"),
                               signal_col = NULL, min_points = 3) {
  weighting <- match.arg(weighting)
  curves <- stats::setNames(
    lapply(met_ids, function(mid)
      fit_calibration_curve(batch_matches, mid, sample_meta, weighting, signal_col, min_points)),
    met_ids)

  rows <- lapply(met_ids, function(mid) {
    curve <- curves[[mid]]
    if (is.null(curve$model)) return(NULL)
    sig_col <- curve$signal_col
    met_matches <- batch_matches[batch_matches$met_id == mid & !is.na(batch_matches[[sig_col]]), ]
    best <- .best_signal_per_met(met_matches, sig_col)
    met_name <- unique(batch_matches$met_name[batch_matches$met_id == mid])[1]

    meta_cols <- intersect(c("sample", "sample_type", "concentration"), names(sample_meta))
    d <- merge(sample_meta[, meta_cols, drop = FALSE], best[, c("sample", sig_col)], by = "sample")
    if (nrow(d) == 0) return(NULL)
    names(d)[names(d) == sig_col] <- "signal"
    if (!"sample_type" %in% names(d)) d$sample_type <- NA_character_
    if (!"concentration" %in% names(d)) d$concentration <- NA_character_

    d$concentration_calc <- (d$signal - curve$intercept) / curve$slope
    d$extrapolated <- d$concentration_calc < curve$conc_range[1] | d$concentration_calc > curve$conc_range[2]
    d$nominal_concentration <- suppressWarnings(as.numeric(d$concentration))
    is_ref <- !is.na(d$sample_type) & d$sample_type %in% c("standard", "quality_control") &
      !is.na(d$nominal_concentration)
    d$percent_re <- ifelse(is_ref,
                            100 * (d$concentration_calc - d$nominal_concentration) / d$nominal_concentration,
                            NA_real_)

    data.frame(met_id = mid, met_name = met_name, sample = d$sample, sample_type = d$sample_type,
               signal = d$signal, concentration_calc = d$concentration_calc,
               extrapolated = d$extrapolated, nominal_concentration = d$nominal_concentration,
               percent_re = d$percent_re, curve_r_squared = curve$r_squared,
               curve_n_points = curve$n_points, curve_weighting = curve$weighting,
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  list(quant = if (length(rows) > 0) do.call(rbind, rows) else data.frame(), curves = curves)
}

## ---- Relative quantification (fold-change over a baseline) -----------------
# For metabolites NOT selected for absolute quantification. Baseline is a
# COHORT mean, not a paired per-subject baseline -- sample_meta has no
# subject/animal-ID column, consistent with how compare_time_series()/
# compare_multi_groups() already treat replicates as independent samples
# rather than repeated measures on the same individual (see that
# function's own header comment).
#
#   mode = "time_series": baseline = mean signal at `reference_timepoint`
#     (the earliest timepoint present, by default -- overridable, since
#     "earliest" and "pre-dose" are not always the same thing: a real
#     report had pre-dose coded as timepoint 0 with 2h/6h post-dose, where
#     0 IS the earliest value and this default is already correct, but a
#     study that only samples DURING dosing, or codes pre-dose as a
#     negative offset, needs to say explicitly which timepoint is the
#     reference rather than have "smallest number" silently decide it),
#     computed PER group/arm if a group column is also present, so e.g. a
#     "treated" arm is normalized to ITS OWN reference-timepoint mean, not
#     pooled with control's.
#   mode = "group": baseline = mean signal in `control_group`; every
#     sample's relative_signal is its own signal over that one number, per
#     the request that "control vs treatment should be relative to
#     control."
#
# Calibration standards/QC/blanks are excluded up front: fold-change over
# a study baseline is a statement about the biological experiment, and
# those sample_type rows aren't part of it even if a group/timepoint value
# was accidentally left on them.
#' Relative-quantify a set of metabolites (fold-change over a baseline)
#'
#' For metabolites NOT selected for absolute quantification: every
#' sample's signal divided by a baseline. Baseline is a COHORT mean, not
#' a paired per-subject baseline -- `sample_meta` has no subject/animal-ID
#' column, consistent with how [compare_time_series()]/
#' [compare_multi_groups()] already treat replicates as independent
#' samples rather than repeated measures on the same individual.
#' Calibration standards/QC/blanks are excluded up front (identified via
#' `sample_type`), even if a group/timepoint value was left on them.
#'
#' @param batch_matches The `ms1_matches` data.frame (see
#'   [match_ms1_batch()]/[annotate_metabolites_batch()]).
#' @param met_ids Metabolite IDs to relative-quantify.
#' @param sample_meta A data.frame with `sample` plus `group` and/or
#'   `timepoint` (and optionally `sample_type`, used to exclude
#'   standards/QC/blanks).
#' @param mode `"time_series"`: baseline is the mean signal at
#'   `reference_timepoint`, computed per group/arm if a `group` column is
#'   also present, so e.g. a "treated" arm normalizes to ITS OWN
#'   reference-timepoint mean rather than pooling with control's.
#'   `"group"`: baseline is the mean signal in `control_group` -- every
#'   sample's `relative_signal` is its own signal over that one number.
#' @param control_group Required when `mode = "group"`: the `group` value
#'   to treat as baseline.
#' @param reference_timepoint Only used when `mode = "time_series"`: the
#'   numeric `timepoint` value to treat as the pre-dose/reference
#'   baseline. `NULL` (the default) uses the earliest timepoint present,
#'   which is only the same thing as "pre-dose" if pre-dose happens to be
#'   coded as the smallest number -- pass this explicitly whenever that
#'   isn't the case (e.g. a study with no baseline draw, or a non-zero/
#'   negative pre-dose code). A sample whose exact timepoint doesn't
#'   match any observed value contributes no rows to `relative_signal`
#'   for that arm (a `NA` baseline), rather than erroring.
#' @param signal_col Which column to use as signal; `NULL` auto-detects
#'   (`"area"` if present and not all-`NA`, else `"intensity"`).
#' @return A data.frame with one row per (metabolite, sample):
#'   `met_id`, `met_name`, the `sample_meta` columns, `signal`,
#'   `baseline_signal`, and `relative_signal` (`signal / baseline_signal`).
#' @seealso [quantify_absolute()] for the calibration-curve half,
#'   [quantify_metabolites()] to dispatch both in one call.
#' @export
quantify_relative <- function(batch_matches, met_ids, sample_meta,
                               mode = c("time_series", "group"),
                               control_group = NULL, reference_timepoint = NULL,
                               signal_col = NULL) {
  mode <- match.arg(mode)
  if (is.null(signal_col)) signal_col <- .auto_signal_col(batch_matches)
  if (mode == "group" && is.null(control_group)) {
    stop("quantify_relative(mode = 'group') needs control_group")
  }

  study_meta <- sample_meta
  if ("sample_type" %in% names(study_meta)) {
    study_meta <- study_meta[.is_study_sample(study_meta$sample_type), , drop = FALSE]
  }

  rows <- lapply(met_ids, function(mid) {
    met_matches <- batch_matches[batch_matches$met_id == mid & !is.na(batch_matches[[signal_col]]), ]
    best <- .best_signal_per_met(met_matches, signal_col)
    met_name <- unique(batch_matches$met_name[batch_matches$met_id == mid])[1]
    d <- merge(study_meta, best[, c("sample", signal_col)], by = "sample")
    if (nrow(d) == 0) return(NULL)
    names(d)[names(d) == signal_col] <- "signal"

    if (mode == "time_series") {
      if (!"timepoint" %in% names(d)) stop("quantify_relative(mode = 'time_series') needs a timepoint column in sample_meta")
      d$.time <- suppressWarnings(as.numeric(d$timepoint))
      d <- d[!is.na(d$.time), ]
      if (nrow(d) == 0) return(NULL)
      d$.arm <- if ("group" %in% names(d)) ifelse(nzchar(d$group), d$group, "") else ""
      ref_time <- if (is.null(reference_timepoint)) min(d$.time) else reference_timepoint
      baseline <- stats::aggregate(signal ~ .arm, data = d[d$.time == ref_time, , drop = FALSE], FUN = mean)
      d$baseline_signal <- baseline$signal[match(d$.arm, baseline$.arm)]
      d$.time <- NULL; d$.arm <- NULL
    } else {
      if (!"group" %in% names(d)) stop("quantify_relative(mode = 'group') needs a group column in sample_meta")
      ctrl_signal <- d$signal[d$group == control_group]
      if (length(ctrl_signal) == 0) return(NULL)
      d$baseline_signal <- mean(ctrl_signal, na.rm = TRUE)
    }
    d$relative_signal <- d$signal / d$baseline_signal
    cbind(met_id = mid, met_name = met_name, d, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

## ---- Combined orchestrator ---------------------------------------------------
# Splits every metabolite present in batch_matches into the user-selected
# absolute_met_ids (calibration curve back-calculation) and everything
# else (fold-change relative to pre-dose/time-0 or to control, per
# `mode`) -- one call site for the whole "some metabolites absolute, the
# rest relative" workflow, rather than making every caller remember to
# call both halves and merge them itself.
#' Quantify every metabolite: absolute for a selected subset, relative for the rest
#'
#' Splits every metabolite present in `batch_matches` into the
#' user-selected `absolute_met_ids` (calibration-curve back-calculation,
#' via [quantify_absolute()]) and everything else (fold-change relative
#' to pre-dose/time-0 or to `control_group`, via [quantify_relative()]) --
#' one call for the whole "some metabolites absolute, the rest relative"
#' workflow.
#'
#' @param batch_matches The `ms1_matches` data.frame (see
#'   [match_ms1_batch()]/[annotate_metabolites_batch()]).
#' @param sample_meta A data.frame with `sample` plus `group`/`timepoint`
#'   (for relative quantification) and `sample_type`/`concentration` (for
#'   absolute quantification) as relevant.
#' @param absolute_met_ids Metabolite IDs to absolute-quantify from their
#'   own calibration curve. Every other metabolite actually present in
#'   `batch_matches` is relative-quantified instead. Defaults to none
#'   (every metabolite relative-quantified).
#' @param mode `"time_series"` or `"group"` -- see [quantify_relative()].
#' @param control_group Required when `mode = "group"` -- see
#'   [quantify_relative()].
#' @param reference_timepoint Only used when `mode = "time_series"` -- see
#'   [quantify_relative()].
#' @param weighting Calibration curve weighting -- see
#'   [fit_calibration_curve()].
#' @param signal_col Which column to use as signal; `NULL` auto-detects.
#' @param min_points Minimum calibration standard points required per
#'   metabolite -- see [fit_calibration_curve()].
#' @return `list(absolute, calibration_curves, relative)`: `absolute` and
#'   `relative` are the `quant`/result data.frames from
#'   [quantify_absolute()]/[quantify_relative()] respectively;
#'   `calibration_curves` is the named list of fitted curves from
#'   [quantify_absolute()].
#' @examples
#' \dontrun{
#' quant <- quantify_metabolites(
#'   batch_results$ms1_matches, sample_meta,
#'   absolute_met_ids = "PARENT", mode = "group", control_group = "control",
#'   weighting = "1/x2")
#' quant$absolute   # PARENT's back-calculated concentrations
#' quant$relative   # every other metabolite's fold-change vs control
#' }
#' @export
quantify_metabolites <- function(batch_matches, sample_meta, absolute_met_ids = character(0),
                                  mode = c("time_series", "group"), control_group = NULL,
                                  reference_timepoint = NULL,
                                  weighting = c("1/x2", "1/x", "none"),
                                  signal_col = NULL, min_points = 3) {
  mode <- match.arg(mode)
  weighting <- match.arg(weighting)
  if (is.null(batch_matches) || nrow(batch_matches) == 0) {
    return(list(absolute = data.frame(), calibration_curves = list(), relative = data.frame()))
  }

  all_met_ids <- unique(batch_matches$met_id)
  absolute_met_ids <- intersect(absolute_met_ids, all_met_ids)
  relative_met_ids <- setdiff(all_met_ids, absolute_met_ids)

  abs_res <- if (length(absolute_met_ids) > 0) {
    quantify_absolute(batch_matches, absolute_met_ids, sample_meta, weighting = weighting,
                       signal_col = signal_col, min_points = min_points)
  } else list(quant = data.frame(), curves = list())

  rel_res <- if (length(relative_met_ids) > 0) {
    quantify_relative(batch_matches, relative_met_ids, sample_meta, mode = mode,
                       control_group = control_group, reference_timepoint = reference_timepoint,
                       signal_col = signal_col)
  } else data.frame()

  list(absolute = abs_res$quant, calibration_curves = abs_res$curves, relative = rel_res)
}
