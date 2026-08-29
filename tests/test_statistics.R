# test_statistics.R -- validate the statistical comparison suite
# (R/statistics.R) against hand-built synthetic abundance data with known
# injected signal, the same known-answer style as tests/test_fragments.R
# and tests/test_mass_isotope.R.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
source(file.path(.pkg_root, "R", "statistics.R"))

cat("==== Statistics suite validation ====\n\n")

## ---- build_abundance_matrix() / abundance_long() ----------------------------
cat("--- build_abundance_matrix() / abundance_long() ---\n")
set.seed(7)
met_ids <- paste0("M", 1:6)
met_names <- paste("metabolite", 1:6)
samples_2g <- c("ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3")
group_2g <- c("control", "control", "control", "treated", "treated", "treated")

batch_matches <- do.call(rbind, lapply(seq_along(met_ids), function(i) {
  do.call(rbind, lapply(seq_along(samples_2g), function(j) {
    base <- 1e5
    # M1 has a strong, deliberate fold change; the rest are flat/noise
    fc <- if (met_ids[i] == "M1" && group_2g[j] == "treated") 4 else 1
    data.frame(met_id = met_ids[i], met_name = met_names[i], kind = "truncation",
               sample = samples_2g[j],
               intensity = base * fc * exp(rnorm(1, 0, 0.08)),
               stringsAsFactors = FALSE)
  }))
}))

abund_mat <- build_abundance_matrix(batch_matches)
cat("Abundance matrix:", nrow(abund_mat), "metabolites x", length(samples_2g), "samples\n")
stopifnot(nrow(abund_mat) == 6)
stopifnot(all(samples_2g %in% names(abund_mat)))

sample_meta_2g <- data.frame(sample = samples_2g, group = group_2g, stringsAsFactors = FALSE)
long_2g <- abundance_long(abund_mat, sample_meta_2g)
cat("Long format rows:", nrow(long_2g), "(expected", 6 * 6, ")\n")
stopifnot(nrow(long_2g) == 36)

## ---- compare_two_groups(): recovers the injected fold change ---------------
cat("\n--- compare_two_groups() ---\n")
two_grp <- compare_two_groups(long_2g, "control", "treated")
print(two_grp[, c("met_id", "mean_a", "mean_b", "log2fc", "p_value", "p_adj")])
m1 <- two_grp[two_grp$met_id == "M1", ]
stopifnot(m1$log2fc > 1.5)
stopifnot(m1$p_adj < 0.05)
cat("M1 (injected 4x fold change) recovered with log2fc =", round(m1$log2fc, 2),
    " p_adj =", signif(m1$p_adj, 3), ": PASS\n")
stopifnot(all(diff(sort(two_grp$p_value)) >= -1e-12))  # p.adjust(BH) preserves rank order
p_order <- order(two_grp$p_value)
stopifnot(all(diff(two_grp$p_adj[p_order]) >= -1e-9))
cat("BH-adjusted p-values are monotonic with raw p-values: PASS\n")

## ---- compare_multi_groups(): 3-group design ---------------------------------
cat("\n--- compare_multi_groups() ---\n")
samples_3g <- c("g1_1", "g1_2", "g2_1", "g2_2", "g3_1", "g3_2")
group_3g <- c("g1", "g1", "g2", "g2", "g3", "g3")
batch_matches_3g <- do.call(rbind, lapply(seq_along(met_ids), function(i) {
  do.call(rbind, lapply(seq_along(samples_3g), function(j) {
    fc <- if (met_ids[i] == "M2") switch(group_3g[j], g1 = 1, g2 = 2, g3 = 4) else 1
    data.frame(met_id = met_ids[i], met_name = met_names[i], kind = "truncation",
               sample = samples_3g[j], intensity = 1e5 * fc * exp(rnorm(1, 0, 0.08)),
               stringsAsFactors = FALSE)
  }))
}))
abund_mat_3g <- build_abundance_matrix(batch_matches_3g)
sample_meta_3g <- data.frame(sample = samples_3g, group = group_3g, stringsAsFactors = FALSE)
long_3g <- abundance_long(abund_mat_3g, sample_meta_3g)

multi <- compare_multi_groups(long_3g, groups = c("g1", "g2", "g3"))
cat("Omnibus rows:", nrow(multi$omnibus), " posthoc rows:", nrow(multi$posthoc), "\n")
stopifnot(nrow(multi$omnibus) == 6)
n_contrasts_per_met <- table(multi$posthoc$met_id)
cat("Posthoc contrasts per metabolite:", paste(unique(n_contrasts_per_met), collapse = ","),
    "(expected choose(3,2) = 3)\n")
stopifnot(all(n_contrasts_per_met == 3))
m2_omni <- multi$omnibus[multi$omnibus$met_id == "M2", ]
stopifnot(m2_omni$p_adj < 0.05)
cat("M2 (injected group trend) significant by ANOVA: PASS\n")

## ---- compare_time_series(): recovers an injected linear slope ---------------
cat("\n--- compare_time_series() ---\n")
samples_ts <- paste0("t", rep(1:4, each = 2), "_", rep(1:2, 4))
timepoint_ts <- rep(c(0, 1, 2, 3), each = 2)
batch_matches_ts <- do.call(rbind, lapply(seq_along(met_ids), function(i) {
  do.call(rbind, lapply(seq_along(samples_ts), function(j) {
    trend <- if (met_ids[i] == "M3") 1 + 0.6 * timepoint_ts[j] else 1
    data.frame(met_id = met_ids[i], met_name = met_names[i], kind = "truncation",
               sample = samples_ts[j], intensity = 1e5 * trend * exp(rnorm(1, 0, 0.05)),
               stringsAsFactors = FALSE)
  }))
}))
abund_mat_ts <- build_abundance_matrix(batch_matches_ts)
sample_meta_ts <- data.frame(sample = samples_ts, timepoint = timepoint_ts, stringsAsFactors = FALSE)
long_ts <- abundance_long(abund_mat_ts, sample_meta_ts)

trend_result <- compare_time_series(long_ts, time_var = "timepoint")
print(trend_result[, c("met_id", "slope", "p_value", "r_squared")])
m3_trend <- trend_result[trend_result$met_id == "M3", ]
stopifnot(m3_trend$slope > 0)
stopifnot(m3_trend$p_value < 0.05)
cat("M3 (injected positive linear trend) recovered with slope =", round(m3_trend$slope, 0),
    " p =", signif(m3_trend$p_value, 3), ": PASS\n")

## ---- Plot functions return ggplot objects -----------------------------------
cat("\n--- Plot functions ---\n")
p1 <- plot_volcano(two_grp)
p2 <- plot_group_boxplot(long_3g, "M2")
long_ts$group <- NULL  # plot_trend expects a `timepoint` column, already present
p3 <- plot_trend(long_ts, "M3")
stopifnot(inherits(p1, "gg"))
stopifnot(inherits(p2, "gg"))
stopifnot(inherits(p3, "gg"))
cat("plot_volcano/plot_group_boxplot/plot_trend all return ggplot objects: PASS\n")

cat("\n==== All statistics tests passed ====\n")
