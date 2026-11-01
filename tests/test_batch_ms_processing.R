# test_batch_ms_processing.R -- validate the batch/multi-sample MS glue
# (R/batch_ms_processing.R) against hand-built synthetic multi-sample data.
# Needs no Python: features/ms2 tables are built directly in R to mirror
# the shape read_batch_features()/read_batch_ms2() would produce from the
# Python pipeline's TSV output.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
source(file.path(.pkg_root, "R", "about.R"))
source(file.path(.pkg_root, "R", "chemistry_dict.R"))
source(file.path(.pkg_root, "R", "oligo_io.R"))
source(file.path(.pkg_root, "R", "metabolites.R"))
source(file.path(.pkg_root, "R", "mass_isotope.R"))
source(file.path(.pkg_root, "R", "fragments.R"))
source(file.path(.pkg_root, "R", "ms_matching.R"))
source(file.path(.pkg_root, "R", "batch_ms_processing.R"))
source(file.path(.pkg_root, "R", "statistics.R"))
source(file.path(.pkg_root, "R", "degradation.R"))

cat("==== Batch MS processing validation ====\n\n")

cat("Python interpreter found:", !is.na(find_python()), "(informational only -- ",
    "the tests below exercise only the pure-R glue, not the Python pipeline)\n\n")

spec <- parse_input(INOTERSEN_TRIPLET)
mets <- generate_metabolites(spec, opts = list(
  oligo_name = "inotersen", max_3p = 5, max_5p = 5, endo = FALSE))
cat("Generated", length(mets), "metabolites\n\n")

## ---- Synthetic multi-sample feature table (mirrors combined_features.tsv) --
cat("--- Building synthetic multi-sample feature table ---\n")
set.seed(42)
samples <- c("ctrl_1", "treat_1")
rows <- list()
matched_feature_ids <- character(0)
unmatched_feature_ids <- character(0)

for (s in samples) {
  for (met in mets[1:6]) {
    info <- metabolite_mass_info(met)
    z <- sample(c(5, 6, 7, 8), 1)
    mz <- (info$mono_mass - z * .PROTON) / z
    ppm_err <- runif(1, -3, 3)
    fid <- paste0(s, "_F", length(rows) + 1)
    rows[[length(rows) + 1]] <- data.frame(
      sample = s, source_file = paste0(s, ".mzML"), feature_id = fid,
      mz = mz * (1 + ppm_err / 1e6), rt = runif(1, 5, 25),
      max_intensity = runif(1, 1e4, 1e6) * (if (s == "treat_1") 2.5 else 1),
      n_scans = sample(3:10, 1), charge = z,
      neutral_mass = info$mono_mass, n_charge_states = sample(2:4, 1),
      mass_cv_ppm = runif(1, 0, 5), rt_start = 0, rt_end = 0,
      area = runif(1, 1e3, 1e5), stringsAsFactors = FALSE)
    matched_feature_ids <- c(matched_feature_ids, fid)
  }
  # deliberately unmatchable noise features (far from any theoretical m/z)
  for (i in 1:5) {
    fid <- paste0(s, "_Fnoise", i)
    rows[[length(rows) + 1]] <- data.frame(
      sample = s, source_file = paste0(s, ".mzML"), feature_id = fid,
      mz = runif(1, 3500, 4000), rt = runif(1, 1, 30),
      max_intensity = runif(1, 100, 5000), n_scans = 1, charge = 3,
      neutral_mass = runif(1, 8000, 9000), n_charge_states = 1,
      mass_cv_ppm = 0, rt_start = 0, rt_end = 0, area = 0,
      stringsAsFactors = FALSE)
    unmatched_feature_ids <- c(unmatched_feature_ids, fid)
  }
}
features <- do.call(rbind, rows)
cat("Synthetic features:", nrow(features), "across", length(samples), "samples\n\n")

## ---- match_ms1_batch() -------------------------------------------------------
cat("--- match_ms1_batch() ---\n")
ms1_matches <- match_ms1_batch(mets, features, ppm_tol = 10, z_range = 3:12,
                                adducts = c("H"), max_oxid = 0, n_iso = 0)
cat("Matches:", nrow(ms1_matches), " samples represented:",
    paste(sort(unique(ms1_matches$sample)), collapse = ","), "\n")
stopifnot(nrow(ms1_matches) > 0)
stopifnot(all(samples %in% ms1_matches$sample))
stopifnot("sample" %in% names(ms1_matches))

## ---- area threading (M3): match_ms1_batch() carries peak area through ----
cat("\n--- area column threading ---\n")
stopifnot("area" %in% names(ms1_matches))
stopifnot(all(!is.na(ms1_matches$area)))
cat("area present on every matched row: PASS\n")

area_matrix <- build_abundance_matrix_area(ms1_matches)
cat("build_abundance_matrix_area():", nrow(area_matrix), "metabolites x",
    length(samples), "samples\n")
stopifnot(nrow(area_matrix) > 0)
stopifnot(all(samples %in% names(area_matrix)))
cat("build_abundance_matrix_area() produces a matrix: PASS\n")

# Features with no area at all (e.g. single-file R-native path) must not
# crash the matrix builder -- it should return an empty data.frame so
# degradation_summary() can fall back to intensity.
features_no_area <- features
features_no_area$area <- NULL
ms1_matches_no_area <- match_ms1_batch(mets, features_no_area, ppm_tol = 10,
                                        z_range = 3:12, adducts = c("H"),
                                        max_oxid = 0, n_iso = 0)
stopifnot("area" %in% names(ms1_matches_no_area))
stopifnot(all(is.na(ms1_matches_no_area$area)))
stopifnot(nrow(build_abundance_matrix_area(ms1_matches_no_area)) == 0)
cat("area-absent features: area column is always NA, matrix is empty: PASS\n")

## ---- unmatched_features_batch(): retained unidentified peaks ---------------
cat("\n--- unmatched_features_batch() ---\n")
unmatched <- unmatched_features_batch(features, ms1_matches)
cat("Unmatched features:", nrow(unmatched), " (expected ", length(unmatched_feature_ids), ")\n", sep = "")
stopifnot(all(unmatched_feature_ids %in% unmatched$feature_id))
stopifnot(length(intersect(matched_feature_ids, unmatched$feature_id)) == 0)
cat("All", length(unmatched_feature_ids), "deliberately-unmatchable features retained,",
    "none of the", length(matched_feature_ids), "matched ones leaked in: PASS\n")

## ---- confirm_ms2_batch(): MS2 confirmation of matched hits ------------------
cat("\n--- confirm_ms2_batch() ---\n")
parent <- mets[[1]]
parent_info <- metabolite_mass_info(parent)
frags <- generate_fragments(parent, z_range = 1:2)
ifrags <- generate_internal_fragments(parent, z_range = 1:2)
all_frags <- c(frags, ifrags)
frag_peaks <- do.call(rbind, lapply(all_frags[1:30], function(f) {
  data.frame(mz = (f$mono_mass - 1 * .PROTON) / 1, intensity = runif(1, 100, 10000),
             stringsAsFactors = FALSE)
}))
frag_peaks <- rbind(frag_peaks,
  data.frame(mz = c(94.9452, 192.9746), intensity = c(8000, 5000)))

parent_z <- ms1_matches$z[ms1_matches$met_id == parent$id & ms1_matches$sample == "ctrl_1"][1]
parent_mz <- ms1_matches$theo_mz[ms1_matches$met_id == parent$id & ms1_matches$sample == "ctrl_1"][1]
# MS2 data supplied for ctrl_1 only -- treat_1 should simply be skipped, not error
ms2_by_sample <- cbind(
  data.frame(sample = "ctrl_1", ms2_scan_id = "scan=ms2_1", rt = 15.0,
             precursor_mz = parent_mz, precursor_z = parent_z, stringsAsFactors = FALSE),
  frag_peaks)

ms2_conf <- confirm_ms2_batch(mets, ms1_matches, ms2_by_sample, frag_tol_ppm = 10, frag_z_range = 1:2)
cat("MS2 confirmations:", nrow(ms2_conf), "\n")
if (nrow(ms2_conf) > 0) {
  print(ms2_conf[, c("sample", "met_name", "confirmation_score", "coverage", "confident")])
}
stopifnot(nrow(ms2_conf) > 0)
stopifnot(all(ms2_conf$sample == "ctrl_1"))  # treat_1 had no MS2 data supplied
cat("MS2 confirmation only for samples with MS2 data: PASS\n")

## ---- Full batch annotation pipeline -----------------------------------------
cat("\n--- annotate_metabolites_batch() ---\n")
batch_results <- annotate_metabolites_batch(mets, features, ms2_by_sample,
                                             ppm_tol = 10, z_range = 3:12,
                                             adducts = c("H"), max_oxid = 0, n_iso = 0)
cat("ms1_matches:", nrow(batch_results$ms1_matches),
    " unmatched:", nrow(batch_results$unmatched),
    " ms2_confirmations:", nrow(batch_results$ms2_confirmations), "\n")
stopifnot(nrow(batch_results$ms1_matches) == nrow(ms1_matches))
stopifnot(nrow(batch_results$unmatched) == nrow(unmatched))

## ---- degradation summary is computed as part of the full pipeline (M4) ---
cat("\n--- annotate_metabolites_batch()$degradation ---\n")
stopifnot(!is.null(batch_results$degradation))
stopifnot(is.list(batch_results$degradation))
stopifnot(all(c("per_sample", "composition", "top_degradants", "signal_used") %in%
              names(batch_results$degradation)))
cat("degradation element present with expected shape: PASS\n")
# compute_degradation = FALSE must skip it (additive, opt-out-able)
batch_results_no_deg <- annotate_metabolites_batch(mets, features, ms2_by_sample,
                                                    ppm_tol = 10, z_range = 3:12,
                                                    adducts = c("H"), max_oxid = 0, n_iso = 0,
                                                    compute_degradation = FALSE)
stopifnot(is.null(batch_results_no_deg$degradation))
cat("compute_degradation = FALSE skips it: PASS\n")

## ---- read_batch_ms2() round-trip (semicolon-list parsing) -------------------
cat("\n--- read_batch_ms2() TSV round-trip ---\n")
tmp_ms2 <- tempfile(fileext = ".tsv")
raw_ms2 <- data.frame(
  sample = "ctrl_1", source_file = "ctrl_1.mzML", ms2_scan_id = "scan=1",
  precursor_mz = 1000.5, precursor_z = 7, rt = 12.3,
  mz_list = "100.1;200.2;300.3", intensity_list = "500;600;700",
  stringsAsFactors = FALSE)
utils::write.table(raw_ms2, tmp_ms2, sep = "\t", row.names = FALSE, quote = FALSE)
parsed_ms2 <- read_batch_ms2(tmp_ms2)
cat("Parsed rows:", nrow(parsed_ms2), "(expected 3)\n")
stopifnot(nrow(parsed_ms2) == 3)
stopifnot(all.equal(parsed_ms2$mz, c(100.1, 200.2, 300.3)))
stopifnot(all.equal(parsed_ms2$intensity, c(500, 600, 700)))
unlink(tmp_ms2)

## ---- Regression: numeric-looking sample names must stay character ----------
# Shiny renames uploads to numeric temp names (e.g. "0.mzML", "1.mzML"), so
# Python's sample column can end up literally "0"/"1". read.delim() auto-
# converts columns that look numeric -- if read_batch_features()/
# read_batch_ms2() didn't force `sample` back to character, a later
# name_map[feats$sample] lookup (see inst/app/app.R) would index the named
# lookup vector POSITIONALLY instead of by name, silently mislabeling every
# sample (caught via manual Shiny testing, not by the tests above, since
# they use non-numeric sample names like "ctrl_1").
cat("\n--- numeric-looking sample names stay character ---\n")
tmp_numeric <- tempfile(fileext = ".tsv")
raw_numeric <- rbind(
  data.frame(sample = "0", source_file = "a", feature_id = "0_F1", mz = 500, rt = 1,
             max_intensity = 1e5, n_scans = 5, charge = 6, neutral_mass = 3000,
             n_charge_states = 2, mass_cv_ppm = 1, rt_start = 0, rt_end = 0, area = 0),
  data.frame(sample = "1", source_file = "b", feature_id = "1_F1", mz = 500, rt = 1,
             max_intensity = 2e5, n_scans = 5, charge = 6, neutral_mass = 3000,
             n_charge_states = 2, mass_cv_ppm = 1, rt_start = 0, rt_end = 0, area = 0))
utils::write.table(raw_numeric, tmp_numeric, sep = "\t", row.names = FALSE, quote = FALSE)
parsed_numeric <- read_batch_features(tmp_numeric)
stopifnot(is.character(parsed_numeric$sample))
name_map <- stats::setNames(c("ctrl_1", "treat_1"), c("0", "1"))
remapped <- unname(name_map[parsed_numeric$sample])
cat("remapped samples:", paste(remapped, collapse = ", "), "\n")
stopifnot(identical(remapped, c("ctrl_1", "treat_1")))
unlink(tmp_numeric)

cat("\n==== All batch_ms_processing tests passed ====\n")
