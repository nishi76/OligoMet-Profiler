# test_multivariate.R -- validate R/multivariate.R: run_pca()/run_hclust()
# recover a known group structure from synthetic data with an injected
# separation, same known-answer style as tests/test_statistics.R.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("chemistry_dict.R", "degradation.R", "statistics.R", "multivariate.R")) {
  source(file.path(.pkg_root, "R", .f))
}

cat("==== multivariate.R validation ====\n\n")

set.seed(42)
met_ids <- paste0("M", 1:8)
met_names <- paste("metabolite", 1:8)
samples <- c(paste0("ctrl_", 1:4), paste0("treat_", 1:4))
group <- rep(c("control", "treated"), each = 4)

## Inject a real group separation on M1-M3 (treated >> control), the rest
## flat noise -- PCA/clustering should pick this up on PC1 and split the
## samples into their true groups.
batch_matches <- do.call(rbind, lapply(seq_along(met_ids), function(i) {
  do.call(rbind, lapply(seq_along(samples), function(j) {
    base <- 1e5
    fc <- if (met_ids[i] %in% c("M1", "M2", "M3") && group[j] == "treated") 8 else 1
    data.frame(met_id = met_ids[i], met_name = met_names[i], kind = "truncation",
               sample = samples[j],
               intensity = base * fc * exp(rnorm(1, 0, 0.05)),
               stringsAsFactors = FALSE)
  }))
}))

sample_meta <- data.frame(sample = samples, sample_type = "unknown",
                           group = group, timepoint = "", stringsAsFactors = FALSE)

## ---- run_pca(): backward-compat / no sample_meta ---------------------------
cat("--- run_pca() with no sample_meta ---\n")
pca0 <- run_pca(batch_matches)
stopifnot(nzchar(pca0$note) == FALSE)
stopifnot(nrow(pca0$scores) == 8)
stopifnot(all(c("PC1", "PC2") %in% names(pca0$scores)))
cat("PASS: runs with no sample_meta, every sample kept\n")

## ---- run_pca(): recovers the injected group separation on PC1 --------------
cat("\n--- run_pca() recovers group separation ---\n")
pca <- run_pca(batch_matches, sample_meta = sample_meta)
stopifnot(nzchar(pca$note) == FALSE)
stopifnot(nrow(pca$scores) == 8)
stopifnot("group" %in% names(pca$scores))
pc1_ctrl <- pca$scores$PC1[pca$scores$group == "control"]
pc1_treat <- pca$scores$PC1[pca$scores$group == "treated"]
cat("PC1 control:", round(pc1_ctrl, 2), "\n")
cat("PC1 treated:", round(pc1_treat, 2), "\n")
stopifnot(abs(mean(pc1_treat) - mean(pc1_ctrl)) > 2 * max(sd(pc1_ctrl), sd(pc1_treat)))
stopifnot(sum(pca$var_explained) <= 1 + 1e-9)
stopifnot(pca$var_explained[["PC1"]] > 0.3)  # the injected effect should dominate
cat("PASS: PC1 cleanly separates control vs treated, PC1 captures the dominant variance\n")

## ---- loadings: M1/M2/M3 dominate PC1 ----------------------------------------
cat("\n--- loadings correctly identify the injected metabolites ---\n")
load1 <- pca$loadings[order(-abs(pca$loadings$PC1)), ]
top3 <- load1$met_id[1:3]
stopifnot(setequal(top3, c("M1", "M2", "M3")))
cat("Top-3 |PC1 loading| metabolites:", paste(top3, collapse = ", "), "(expected M1/M2/M3): PASS\n")

## ---- run_pca(): excludes calibration standards/QC --------------------------
cat("\n--- run_pca() excludes standards/QC ---\n")
std_matches <- rbind(batch_matches,
  data.frame(met_id = "M1", met_name = "metabolite 1", kind = "truncation",
             sample = "std_1", intensity = 999999, stringsAsFactors = FALSE))
std_meta <- rbind(sample_meta,
  data.frame(sample = "std_1", sample_type = "standard", group = "", timepoint = "", stringsAsFactors = FALSE))
pca_std <- run_pca(std_matches, sample_meta = std_meta)
stopifnot(!"std_1" %in% pca_std$scores$sample)
stopifnot(nrow(pca_std$scores) == 8)
cat("PASS: standard-type sample excluded from PCA scores\n")

## ---- run_pca(): not enough samples -> empty result with a note -------------
cat("\n--- run_pca() with too few samples ---\n")
tiny <- batch_matches[batch_matches$sample %in% c("ctrl_1", "ctrl_2"), ]
pca_tiny <- run_pca(tiny, min_samples = 3)
stopifnot(nrow(pca_tiny$scores) == 0)
stopifnot(nzchar(pca_tiny$note))
cat("note:", pca_tiny$note, "\n")
cat("PASS: too few samples -> empty result with an explanatory note, not an error\n")

## ---- run_hclust(): recovers the same 2-cluster structure --------------------
cat("\n--- run_hclust() recovers control/treated as 2 clusters ---\n")
hc_res <- run_hclust(batch_matches, sample_meta = sample_meta, k = 2)
stopifnot(nzchar(hc_res$note) == FALSE)
stopifnot(inherits(hc_res$hclust, "hclust"))
stopifnot(nrow(hc_res$clusters) == 8)
ctrl_clusters <- unique(hc_res$clusters$cluster[hc_res$clusters$group == "control"])
treat_clusters <- unique(hc_res$clusters$cluster[hc_res$clusters$group == "treated"])
stopifnot(length(ctrl_clusters) == 1)
stopifnot(length(treat_clusters) == 1)
stopifnot(ctrl_clusters != treat_clusters)
cat("control samples all in cluster", ctrl_clusters, ", treated all in cluster", treat_clusters, ": PASS\n")

## ---- run_hclust(): default k, cutree() never errors for small n ------------
cat("\n--- run_hclust() default k on a small sample set ---\n")
hc_small <- run_hclust(batch_matches[batch_matches$sample %in% samples[1:3], ], min_samples = 3)
stopifnot(nzchar(hc_small$note) == FALSE)
stopifnot(hc_small$k < 3)
cat("k =", hc_small$k, "for n=3 samples (clamped below n): PASS\n")

## ---- plot_pca_scores() / plot_dendrogram() return without error ------------
cat("\n--- plotting functions ---\n")
p <- plot_pca_scores(pca)
stopifnot(inherits(p, "ggplot"))
p_empty <- plot_pca_scores(list(scores = data.frame(), note = "no data"))
stopifnot(inherits(p_empty, "ggplot"))
tmp_png <- tempfile(fileext = ".png")
grDevices::png(tmp_png)
plot_dendrogram(hc_res)
plot_dendrogram(list(hclust = NULL, note = "no data"))
grDevices::dev.off()
stopifnot(file.exists(tmp_png))
cat("PASS: PCA scores plot (populated + empty) and dendrogram (populated + empty) all render without error\n")

cat("\n==== All multivariate tests passed ====\n")
