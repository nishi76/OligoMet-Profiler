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
build_abundance_matrix <- function(batch_matches) {
  if (is.null(batch_matches) || nrow(batch_matches) == 0) return(data.frame())
  met_info <- unique(batch_matches[, c("met_id", "met_name", "kind")])
  met_info <- met_info[order(met_info$met_id), ]
  samples <- unique(batch_matches$sample)
  out <- met_info
  for (s in samples) {
    sub <- batch_matches[batch_matches$sample == s, ]
    best <- stats::aggregate(intensity ~ met_id, data = sub, FUN = max)
    out[[s]] <- best$intensity[match(out$met_id, best$met_id)]
  }
  out
}

# Long-format companion: sample_meta is data.frame(sample, group) and/or
# data.frame(sample, timepoint) -- whichever columns are present besides
# `sample` are carried through untouched.
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

## ---- Two-group comparison ---------------------------------------------------
compare_two_groups <- function(abundance_long, group_a, group_b, min_n = 2) {
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
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out
}

## ---- Three-or-more-group comparison -----------------------------------------
compare_multi_groups <- function(abundance_long, groups = NULL, min_n = 2) {
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
  omnibus$p_adj <- stats::p.adjust(omnibus$p_value, method = "BH")
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
compare_time_series <- function(abundance_long, time_var = "timepoint", min_n = 2) {
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
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
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
