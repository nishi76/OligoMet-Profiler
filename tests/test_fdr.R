# test_fdr.R -- validate R/fdr.R: .make_decoy_spec() and estimate_fdr()

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R", "ms_matching.R", "fdr.R")) {
  source(file.path(.pkg_root, "R", .f))
}

cat("==== fdr.R validation ====\n\n")

spec <- parse_input(INOTERSEN_TRIPLET)
mets <- generate_metabolites(spec, opts = list(
  oligo_name = "inotersen", max_3p = 5, max_5p = 5, endo = FALSE))
cat("Generated", length(mets), "metabolites\n\n")

## ---- .make_decoy_spec(): decoy-generation sanity -------------------------
cat("--- .make_decoy_spec() ---\n")
set.seed(7)
decoy_spec <- .make_decoy_spec(spec, seed = 7)
stopifnot(identical(decoy_spec$sugars, spec$sugars))
stopifnot(identical(decoy_spec$linkages, spec$linkages))
stopifnot(identical(decoy_spec$conj5, spec$conj5))
stopifnot(identical(decoy_spec$conj3, spec$conj3))
stopifnot(identical(decoy_spec$n, spec$n))
stopifnot(setequal(decoy_spec$bases, spec$bases))
# A random permutation of a low-diversity 4-symbol alphabet has a small but
# real chance of reproducing the original order -- for an 18-20mer this
# probability is negligible, but retry once rather than risk a flaky test.
if (identical(decoy_spec$bases, spec$bases)) {
  decoy_spec <- .make_decoy_spec(spec, seed = 8)
}
stopifnot(!identical(decoy_spec$bases, spec$bases))
cat("Decoy preserves sugars/linkages/conjugates/n, permutes bases: PASS\n")

## ---- estimate_fdr(): smoke test -------------------------------------------
cat("\n--- estimate_fdr() smoke test ---\n")
set.seed(123)
ms1_features <- list()
for (met in mets[1:8]) {
  info <- metabolite_mass_info(met)
  for (z in c(5, 6, 7, 8)) {
    mz <- (info$mono_mass - z * .PROTON) / z
    ppm_err <- runif(1, -3, 3)
    obs_mz <- mz * (1 + ppm_err / 1e6)
    ms1_features[[length(ms1_features) + 1]] <- data.frame(
      mz = obs_mz, rt = runif(1, 5, 25),
      max_intensity = runif(1, 1e4, 1e6),
      n_scans = sample(3:10, 1),
      stringsAsFactors = FALSE)
  }
}
for (i in 1:20) {
  ms1_features[[length(ms1_features) + 1]] <- data.frame(
    mz = runif(1, 400, 2000), rt = runif(1, 1, 30),
    max_intensity = runif(1, 100, 5000),
    n_scans = 1, stringsAsFactors = FALSE)
}
ms1_features <- do.call(rbind, ms1_features)
ms1_features <- ms1_features[order(ms1_features$mz), ]

opts <- list(oligo_name = "inotersen", max_3p = 5, max_5p = 5, endo = FALSE)
result <- estimate_fdr(spec, mets, ms1_features, ms2_data = NULL,
                        n_decoys = 2, opts = opts,
                        ppm_tol = 10, z_range = 5:8, adducts = c("H"),
                        max_oxid = 0, seed = 42)
cat("n_target_matches:", result$n_target_matches, "\n")
cat("n_decoy_matches (mean of", result$n_decoys, "decoys):", result$n_decoy_matches, "\n")
cat("fdr:", result$fdr, "\n")
stopifnot(is.list(result))
stopifnot(!is.null(result$fdr))
stopifnot(is.na(result$fdr) || result$fdr >= 0)
cat("estimate_fdr() runs without error, returns a sane fdr: PASS\n")

## ---- Bug-detector: decoys must not match at least as well as the target --
cat("\n--- Decoys must not out-match the real spec ---\n")
stopifnot(result$n_target_matches > 0)  # the synthetic features above are seeded at real hits
stopifnot(result$n_decoy_matches <= result$n_target_matches)
cat("n_decoy_matches (", result$n_decoy_matches, ") <= n_target_matches (",
    result$n_target_matches, "): PASS\n")

cat("\n==== All fdr tests passed ====\n")
