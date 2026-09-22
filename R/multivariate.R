# =============================================================================
# multivariate.R
# PCA and hierarchical clustering across metabolites/samples, on top of
# annotate_metabolites_batch()'s ms1_matches -- the "Multivariate" tab's
# backing functions. Dependency-light on purpose, same as the rest of the
# statistics suite: stats::prcomp/hclust/cutree/dist (base R) plus ggplot2,
# both already Imports. No new dependency.
#
# Both PCA and clustering work on the SAME samples x metabolites matrix --
# build_abundance_matrix() (R/statistics.R) transposed, missing (metabolite,
# sample) combinations treated as zero signal (a real absence, not a random
# missing value -- a different, deliberate choice from build_abundance_matrix()'s
# own NA-for-no-match convention), log2-transformed by default (MS intensity
# is heavy right-skewed), with calibration standards/QC/blanks excluded the
# same way degradation_summary()/quantify_relative() already exclude them
# (.is_study_sample(), R/chemistry_dict.R) -- assay-performance samples have
# no place clustering alongside the biological experiment.
# =============================================================================

## ---- Shared matrix builder ---------------------------------------------------
# Not exported -- run_pca()/run_hclust() are the public entry points, both
# built on the exact same (filtered, zeroed, log-transformed) matrix so a
# PCA plot and a dendrogram of the same run are never looking at subtly
# different data. Scaling (unit variance per metabolite) is deliberately
# NOT done here -- prcomp() has its own scale. argument and doing it twice
# would double-apply, so each caller below applies it at its own call site.
.build_multivariate_matrix <- function(batch_matches, sample_meta, met_ids,
                                        signal_col, log_transform, min_samples) {
  empty <- list(mat = matrix(numeric(0), nrow = 0, ncol = 0), met_info = data.frame(),
                dropped_zero_variance = character(0), note = "")
  if (is.null(batch_matches) || nrow(batch_matches) == 0) {
    empty$note <- "no MS1 matches to run on"
    return(empty)
  }
  bm <- batch_matches
  if (!is.null(met_ids)) bm <- bm[bm$met_id %in% met_ids, ]
  if (nrow(bm) == 0) {
    empty$note <- "none of the selected metabolites have any matches"
    return(empty)
  }

  sig_col <- signal_col %||% .auto_signal_col(bm)
  abund <- if (sig_col == "area") build_abundance_matrix_area(bm) else build_abundance_matrix(bm)
  if (nrow(abund) == 0) {
    empty$note <- "no metabolites matched"
    return(empty)
  }

  samples <- setdiff(names(abund), c("met_id", "met_name", "kind"))
  if (!is.null(sample_meta) && "sample_type" %in% names(sample_meta)) {
    keep <- sample_meta$sample[.is_study_sample(sample_meta$sample_type)]
    samples <- intersect(samples, keep)
  }
  if (length(samples) < min_samples) {
    empty$note <- sprintf(
      "only %d study sample(s) with matched signal (calibration standards/QC/blanks excluded) -- need >= %d",
      length(samples), min_samples)
    return(empty)
  }

  met_info <- abund[, c("met_id", "met_name", "kind")]
  mat <- as.matrix(abund[, samples, drop = FALSE])
  mat[is.na(mat)] <- 0
  mat <- t(mat)
  rownames(mat) <- samples
  colnames(mat) <- met_info$met_id

  if (log_transform) mat <- log2(mat + 1)

  sds <- apply(mat, 2, stats::sd)
  dropped <- met_info$met_id[is.na(sds) | sds == 0]
  if (length(dropped) > 0) {
    keep_cols <- !(met_info$met_id %in% dropped)
    mat <- mat[, keep_cols, drop = FALSE]
    met_info <- met_info[keep_cols, ]
  }
  if (ncol(mat) < 2) {
    empty$note <- "fewer than 2 metabolites vary across the study samples (constant/zero everywhere else)"
    empty$dropped_zero_variance <- dropped
    return(empty)
  }

  list(mat = mat, met_info = met_info, dropped_zero_variance = dropped, note = "")
}

## ---- PCA -----------------------------------------------------------------------
#' Principal component analysis of metabolite signal across samples
#'
#' Runs PCA on the samples x metabolites abundance matrix (samples as
#' observations, metabolites as variables) -- do study samples separate by
#' group/timepoint in metabolite-profile space? Calibration standards,
#' quality-control, and reagent/matrix-blank samples are excluded via
#' `sample_meta$sample_type` (see [degradation_summary()] for the same
#' rule); a missing (metabolite, sample) match is treated as zero signal
#' (a real absence), not dropped. Metabolites with zero variance across
#' the kept samples (after any log transform) are dropped before fitting,
#' since [stats::prcomp()] with `scale. = TRUE` divides by standard
#' deviation.
#'
#' @param batch_matches The `ms1_matches` data.frame from
#'   [annotate_metabolites_batch()]/[match_ms1_batch()].
#' @param sample_meta Optional data.frame with `sample` plus `sample_type`
#'   (to exclude standards/QC/blanks) and/or `group`/`timepoint` (joined
#'   onto `scores` for coloring a plot). `NULL` keeps every sample.
#' @param met_ids Optional subset of metabolite IDs to include; `NULL`
#'   uses every metabolite present in `batch_matches`.
#' @param signal_col Which column to use as signal; `NULL` auto-detects
#'   (`"area"` if present and not all-`NA`, else `"intensity"`).
#' @param log_transform Log2-transform (with a +1 pseudo-count) before
#'   fitting -- MS intensity/area is heavy right-skewed, so this is `TRUE`
#'   by default.
#' @param scale Unit-variance-scale each metabolite before fitting
#'   (`stats::prcomp(scale. = )`) -- `TRUE` by default, the usual choice
#'   when variables (metabolites here) differ in scale.
#' @param min_samples Minimum study samples (after excluding calibration/
#'   QC/blanks) required to run PCA at all. Below this, an empty result
#'   with a `note` explaining why is returned instead of erroring.
#' @return A list: `scores` (`sample`, `PC1`, `PC2`, ..., plus `group`/
#'   `timepoint` if `sample_meta` supplied them), `loadings` (`met_id`,
#'   `met_name`, `PC1`, `PC2`, ...), `var_explained` (named numeric
#'   vector, proportion of variance per PC), `dropped_zero_variance`
#'   (metabolite IDs excluded for having no variance across the kept
#'   samples), and `note` (empty string on success, else why the result
#'   is empty).
#' @seealso [plot_pca_scores()], [run_hclust()] for the same matrix's
#'   dendrogram.
#' @export
run_pca <- function(batch_matches, sample_meta = NULL, met_ids = NULL, signal_col = NULL,
                     log_transform = TRUE, scale = TRUE, min_samples = 3) {
  b <- .build_multivariate_matrix(batch_matches, sample_meta, met_ids, signal_col,
                                   log_transform, min_samples)
  empty <- list(scores = data.frame(), loadings = data.frame(), var_explained = numeric(0),
                dropped_zero_variance = b$dropped_zero_variance, note = b$note)
  if (nzchar(b$note)) return(empty)

  pr <- stats::prcomp(b$mat, center = TRUE, scale. = scale)
  var_explained <- (pr$sdev^2) / sum(pr$sdev^2)
  names(var_explained) <- colnames(pr$x)

  scores <- data.frame(sample = rownames(b$mat), pr$x, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(sample_meta)) {
    ctx_cols <- intersect(c("sample", "group", "timepoint"), names(sample_meta))
    if (length(ctx_cols) > 1) {
      scores <- merge(scores, unique(sample_meta[, ctx_cols, drop = FALSE]), by = "sample", all.x = TRUE)
    }
  }
  loadings <- data.frame(met_id = b$met_info$met_id, met_name = b$met_info$met_name,
                          pr$rotation, stringsAsFactors = FALSE, check.names = FALSE)

  list(scores = scores, loadings = loadings, var_explained = var_explained,
       dropped_zero_variance = b$dropped_zero_variance, note = "")
}

## ---- PCA plot ------------------------------------------------------------------
#' Plot PCA sample scores (PC1 vs PC2)
#'
#' @param pca_result A [run_pca()] result.
#' @param color_by Column of `pca_result$scores` to color points by;
#'   `NULL` (the default) auto-picks `group` if present and non-blank,
#'   else `timepoint` if present, else leaves points uncolored.
#' @return A ggplot object.
#' @export
plot_pca_scores <- function(pca_result, color_by = NULL) {
  df <- pca_result$scores
  if (is.null(df) || nrow(df) == 0 || !all(c("PC1", "PC2") %in% names(df))) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = if (nzchar(pca_result$note %||% "")) pca_result$note else "No PCA result to plot"))
  }
  if (is.null(color_by)) {
    color_by <- if ("group" %in% names(df) && any(nzchar(df$group))) "group"
                else if ("timepoint" %in% names(df) && any(nzchar(as.character(df$timepoint)))) "timepoint"
                else NA_character_
  }
  has_color <- !is.na(color_by) && color_by %in% names(df)
  df$.color <- if (has_color) as.character(df[[color_by]]) else "all samples"
  pc1_pct <- round(100 * unname(pca_result$var_explained["PC1"]), 1)
  pc2_pct <- round(100 * unname(pca_result$var_explained["PC2"]), 1)
  n_col <- length(unique(df$.color))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$.color, label = .data$sample)) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::geom_text(vjust = -0.9, size = 3, show.legend = FALSE, check_overlap = TRUE) +
    ggplot2::scale_color_manual(values = grDevices::colorRampPalette(
      c("#0279EE", "#FF9400", "#75A025", "#FD9BED", "#E9ED4C"))(n_col),
      name = if (has_color) color_by else NULL, guide = if (has_color) "legend" else "none") +
    ggplot2::labs(x = paste0("PC1 (", pc1_pct, "%)"), y = paste0("PC2 (", pc2_pct, "%)"),
                  title = "PCA -- sample scores") +
    ggplot2::theme_minimal(base_size = 11, base_family = "Liberation Sans")
}

## ---- Hierarchical clustering ---------------------------------------------------
#' Hierarchical clustering of samples by metabolite profile
#'
#' Clusters samples on the same (optionally log-transformed, unit-variance
#' scaled) samples x metabolites matrix [run_pca()] uses -- distance
#' between samples' metabolite profiles, then agglomerative clustering.
#' Calibration standards/QC/blanks are excluded the same way (see
#' [run_pca()]/[degradation_summary()]).
#'
#' @inheritParams run_pca
#' @param dist_method Distance metric, passed to [stats::dist()].
#' @param clust_method Agglomeration method, passed to [stats::hclust()];
#'   `"ward.D2"` (the default) minimizes within-cluster variance, a
#'   common default for this kind of profiling data.
#' @param k Number of clusters to cut the tree into (`stats::cutree()`).
#'   `NULL` (the default) picks `min(4, floor(n_samples / 2))`, clamped to
#'   at least 2.
#' @return A list: `hclust` (the fitted `hclust` object, or `NULL` if
#'   there wasn't enough data), `clusters` (`sample`, `cluster`, plus
#'   `group`/`timepoint` if `sample_meta` supplied them -- for comparing
#'   the cut against the known study design), `k` (clusters actually
#'   used), `dropped_zero_variance`, and `note` (empty on success).
#' @seealso [plot_dendrogram()], [run_pca()] for the same matrix's PCA.
#' @export
run_hclust <- function(batch_matches, sample_meta = NULL, met_ids = NULL, signal_col = NULL,
                        log_transform = TRUE, scale = TRUE, dist_method = "euclidean",
                        clust_method = "ward.D2", k = NULL, min_samples = 3) {
  b <- .build_multivariate_matrix(batch_matches, sample_meta, met_ids, signal_col,
                                   log_transform, min_samples)
  empty <- list(hclust = NULL, clusters = data.frame(), k = 0L,
                dropped_zero_variance = b$dropped_zero_variance, note = b$note)
  if (nzchar(b$note)) return(empty)

  mat <- if (scale) scale(b$mat, center = TRUE, scale = TRUE) else b$mat
  d <- stats::dist(mat, method = dist_method)
  hc <- stats::hclust(d, method = clust_method)
  k <- k %||% max(2, min(4, floor(nrow(mat) / 2)))
  k <- min(k, nrow(mat) - 1)  # cutree() errors if k >= n_samples

  cl <- stats::cutree(hc, k = k)
  clusters <- data.frame(sample = names(cl), cluster = as.integer(cl), stringsAsFactors = FALSE)
  if (!is.null(sample_meta)) {
    ctx_cols <- intersect(c("sample", "group", "timepoint"), names(sample_meta))
    if (length(ctx_cols) > 1) {
      clusters <- merge(clusters, unique(sample_meta[, ctx_cols, drop = FALSE]), by = "sample", all.x = TRUE)
    }
  }
  clusters <- clusters[order(clusters$cluster, clusters$sample), ]

  list(hclust = hc, clusters = clusters, k = k,
       dropped_zero_variance = b$dropped_zero_variance, note = "")
}

## ---- Dendrogram plot ------------------------------------------------------------
#' Plot a dendrogram from a run_hclust() result
#'
#' Base-graphics dendrogram (not ggplot2 -- `hclust` objects have their
#' own `plot` method, and re-deriving that in ggplot2 buys nothing here),
#' with the `k`-cluster cut outlined if requested.
#'
#' @param hclust_result A [run_hclust()] result (or a raw `hclust`
#'   object, for direct use outside the Shiny app).
#' @param k Number of clusters to outline with [stats::rect.hclust()];
#'   `NULL` skips outlining. Defaults to the `k` a [run_hclust()] result
#'   carries, if given one of those.
#' @return `NULL`, invisibly -- called for its plotting side effect.
#' @export
plot_dendrogram <- function(hclust_result, k = NULL) {
  hc <- if (inherits(hclust_result, "hclust")) hclust_result else hclust_result$hclust
  if (is.null(k) && is.list(hclust_result) && !inherits(hclust_result, "hclust")) k <- hclust_result$k
  if (is.null(hc)) {
    graphics::plot.new()
    graphics::title(main = if (is.list(hclust_result) && nzchar(hclust_result$note %||% "")) {
      hclust_result$note
    } else "No clustering result to plot")
    return(invisible(NULL))
  }
  graphics::plot(hc, main = "Hierarchical clustering (samples)", xlab = "", sub = "", hang = -1)
  if (!is.null(k) && k >= 2) stats::rect.hclust(hc, k = k, border = "#0279EE")
  invisible(NULL)
}
