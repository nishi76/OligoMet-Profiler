# test_ms_matching.R -- validate MS import + matching module
# Uses synthetic MS1 features and MS2 peaks derived from the inotersen library

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
source(file.path(.pkg_root, "R", "spectra_io.R"))

cat("==== MS matching validation ====\n\n")

# Check msconvert availability
msc <- find_msconvert()
cat("msconvert available:", !is.na(msc), "\n")

# Parse the inotersen reference case and generate metabolite library
spec <- parse_input(INOTERSEN_TRIPLET)
mets <- generate_metabolites(spec, opts = list(
  oligo_name = "inotersen", max_3p = 5, max_5p = 5, endo = FALSE))
cat("Generated", length(mets), "metabolites\n\n")

# ---- Create synthetic MS1 features from theoretical m/z values ----
cat("--- Creating synthetic MS1 features ---\n")
set.seed(123)
ms1_features <- list()
for (met in mets[1:8]) {  # parent + some truncations
  info <- metabolite_mass_info(met)
  for (z in c(5, 6, 7, 8)) {
    mz <- (info$mono_mass - z * .PROTON) / z
    # Add small random mass error (1-5 ppm)
    ppm_err <- runif(1, -3, 3)
    obs_mz <- mz * (1 + ppm_err / 1e6)
    ms1_features[[length(ms1_features) + 1]] <- data.frame(
      mz = obs_mz, rt = runif(1, 5, 25),
      max_intensity = runif(1, 1e4, 1e6),
      n_scans = sample(3:10, 1),
      stringsAsFactors = FALSE)
  }
}
# Add some noise features
for (i in 1:20) {
  ms1_features[[length(ms1_features) + 1]] <- data.frame(
    mz = runif(1, 400, 2000), rt = runif(1, 1, 30),
    max_intensity = runif(1, 100, 5000),
    n_scans = 1, stringsAsFactors = FALSE)
}
ms1_features <- do.call(rbind, ms1_features)
ms1_features <- ms1_features[order(ms1_features$mz), ]
cat("Synthetic MS1 features:", nrow(ms1_features), "\n\n")

# ---- Test MS1 matching ----
cat("--- MS1 targeted matching ---\n")
ms1_matches <- match_ms1(mets, ms1_features, ppm_tol = 10,
                          z_range = 3:12, adducts = c("H"),
                          max_oxid = 0, h_offset = 0, n_iso = 0)
cat("MS1 matches:", nrow(ms1_matches), "\n")
if (nrow(ms1_matches) > 0) {
  cat("Matched metabolites:", length(unique(ms1_matches$met_id)), "\n")
  cat("Charge states found:", paste(sort(unique(ms1_matches$z)), collapse=","), "\n")
  cat("Median ppm error:", round(median(ms1_matches$ppm_error), 3), "\n")
  cat("\nFirst 10 matches:\n")
  print(head(ms1_matches[, c("met_name", "z", "theo_mz", "obs_mz", "ppm_error")], 10))
}

# ---- Test envelope consistency ----
cat("\n--- Charge envelope consistency ---\n")
env <- envelope_consistency(ms1_matches)
cat("Envelope groups:", nrow(env), "\n")
if (nrow(env) > 0) {
  print(head(env, 5))
}

# ---- Test MS2 matching with synthetic spectra ----
cat("\n--- MS2 fragment matching ---\n")
# Create synthetic MS2 spectrum for the parent metabolite
parent <- mets[[1]]
frags <- generate_fragments(parent, z_range = 1:2)
ifrags <- generate_internal_fragments(parent, z_range = 1:2)
all_frags <- c(frags, ifrags)

# Build synthetic MS2 peaks from fragments (z=1)
ms2_peaks <- do.call(rbind, lapply(all_frags[1:30], function(f) {
  data.frame(mz = (f$mono_mass - 1 * .PROTON) / 1,
             intensity = runif(1, 100, 10000),
             stringsAsFactors = FALSE)
}))
# Add PS diagnostic ions
ms2_peaks <- rbind(ms2_peaks,
  data.frame(mz = 94.9452, intensity = 8000),
  data.frame(mz = 192.9746, intensity = 5000))
# Add noise
ms2_peaks <- rbind(ms2_peaks,
  data.frame(mz = runif(15, 100, 3000), intensity = runif(15, 50, 500)))
ms2_peaks <- ms2_peaks[order(ms2_peaks$mz), ]

# Match fragments
matched_frags <- match_fragments(all_frags, ms2_peaks, tol_ppm = 10, z_range = 1:2)
cat("Matched fragments:", nrow(matched_frags), "\n")
if (nrow(matched_frags) > 0) {
  cat("Ion types:", paste(unique(matched_frags$ion_type), collapse=", "), "\n")
  cat("Covered sites:", paste(sort(unique(matched_frags$cleavage_site)), collapse=","), "\n")
}
stopifnot(length(unique(matched_frags$obs_mz)) == nrow(matched_frags))
stopifnot(identical(matched_frags,
                     match_fragments(all_frags, ms2_peaks, tol_ppm = 10, z_range = 1:2,
                                      adducts = c("H"))))
cat("No multi-assigned peaks; adducts=c('H') default is byte-identical: PASS\n")

# Check diagnostics
diags <- check_ps_diagnostic(ms2_peaks)
cat("PS diagnostic ions:", nrow(diags), "\n")

# Confirmation score
score <- confirmation_score(matched_frags, parent$n, diags)
cat("\nConfirmation score:", score$total_score, "/ 100\n")
cat("Coverage:", sprintf("%.1f%%", score$coverage * 100), "\n")
cat("Confident:", score$confident, "\n")

# ---- Test full annotation pipeline ----
cat("\n--- Full annotation pipeline ---\n")
# Create MS2 data in the format expected by find_ms2_spectra
parent_mz <- (metabolite_mass_info(parent)$mono_mass - 6 * .PROTON) / 6
ms2_data <- cbind(data.frame(rt = 15.0, precursor_mz = parent_mz,
                              precursor_z = 6L), ms2_peaks)

results <- annotate_metabolites(mets, ms1_features, ms2_data,
                                  ppm_tol = 10, z_range = 3:12,
                                  adducts = c("H"), max_oxid = 0,
                                  frag_tol_ppm = 10, frag_z_range = 1:2,
                                  n_iso = 0, use_envipat = FALSE)
cat("Summary rows:", nrow(results$summary), "\n")
if (nrow(results$summary) > 0) {
  cat("\nSummary:\n")
  print(results$summary[, c("met_name", "n_ms1_matches", "best_ppm",
                             "n_charge_states", "has_ms2", "ms2_score",
                             "ms2_coverage", "confident")])
}

# ---- Test PRM inclusion list ----
cat("\n--- PRM inclusion list ---\n")
prm <- prm_inclusion_list(mets, z_range = 3:12, max_oxid = 2)
cat("PRM entries:", nrow(prm), "\n")
cat("Unique metabolites:", length(unique(prm$met_id)), "\n")
cat("Oxidation levels:", paste(sort(unique(prm$k_oxid)), collapse=","), "\n")
cat("m/z range:", sprintf("%.2f - %.2f\n", min(prm$precursor_mz), max(prm$precursor_mz)))

# ---- Test envelope_consistency() honors a non-zero h_offset ----
# Previously hardcoded to 0 internally with no way to override it; with
# the legacy workbook h_offset (3.0046), every reported mean_mass was off
# by exactly that constant (uniform across charge states -- it doesn't
# scale with z, since h_offset is subtracted once in the back-calculation
# formula, not divided by z -- so the *consistency/cv score* itself is
# actually unaffected, but the reported neutral mass is still wrong).
cat("\n--- envelope_consistency() h_offset ---\n")
true_M <- 5000; h_offset_legacy <- 3.0046
zs_ec <- c(5, 6, 7, 8)
obs_mz_ec <- (true_M + h_offset_legacy - zs_ec * .PROTON) / zs_ec
matches_ec <- data.frame(met_id = "M1", k_oxid = 0, adduct = "H", z = zs_ec,
                          obs_mz = obs_mz_ec, stringsAsFactors = FALSE)
env_correct <- envelope_consistency(matches_ec, h_offset = h_offset_legacy)
env_wrong <- envelope_consistency(matches_ec, h_offset = 0)
env_default <- envelope_consistency(matches_ec)
cat("true mass:", true_M, " | h_offset=3.0046:", env_correct$mean_mass,
    " | h_offset=0:", env_wrong$mean_mass, " | default (no arg):", env_default$mean_mass, "\n")
stopifnot(abs(env_correct$mean_mass - true_M) < 1e-6)
stopifnot(abs(env_wrong$mean_mass - (true_M + h_offset_legacy)) < 1e-6)
stopifnot(abs(env_default$mean_mass - env_wrong$mean_mass) < 1e-9)  # default stays h_offset=0, unchanged for standard mode
cat("envelope_consistency(h_offset=) parameter verified: PASS\n")

# ---- Test extract_ms1_features() clusters by RT, not m/z alone ----
# Previously grouped purely by m/z proximity across the WHOLE file,
# merging peaks that share an m/z but come from unrelated elution events
# anywhere in the run into one feature with a meaningless averaged RT.
cat("\n--- extract_ms1_features() RT tolerance ---\n")
same_mz <- 1000.0000
peaks_rt <- data.frame(
  mz = c(same_mz, same_mz, same_mz, same_mz + 0.0002, same_mz + 0.0002, same_mz + 0.0002),
  rt = c(10.00, 10.05, 10.10,    # one real elution event, tightly spaced
         2.00, 40.00, 40.03),    # rt=2.00 is an unrelated, far-away event;
                                 # rt=40.00/40.03 is a second real event
  intensity = c(5e4, 6e4, 5.5e4, 3e4, 4e4, 4.2e4),
  stringsAsFactors = FALSE)
feat_rt <- extract_ms1_features(peaks_rt, ppm = 10, min_intensity = 100, rt_tol = 0.15)
cat("features:", nrow(feat_rt), "\n")
print(feat_rt[, c("mz", "rt", "max_intensity", "n_scans")])
stopifnot(nrow(feat_rt) == 3)  # 2 events at same_mz's cluster + 1 lone point at same_mz+0.0002... see below
# the tight rt=10.00/10.05/10.10 triplet must merge into ONE feature
merged <- feat_rt[abs(feat_rt$rt - 10.05) < 0.1, ]
stopifnot(nrow(merged) == 1 && merged$n_scans == 3)
# rt=2.00 (unrelated, far from everything) must NOT be merged into any other feature
lone <- feat_rt[abs(feat_rt$rt - 2.00) < 0.01, ]
stopifnot(nrow(lone) == 1 && lone$n_scans == 1)
# rt=40.00/40.03 (its own close pair, same m/z bin) must merge into ONE feature, separate from the rt=10 cluster
pair <- feat_rt[abs(feat_rt$rt - 40.015) < 0.1, ]
stopifnot(nrow(pair) == 1 && pair$n_scans == 2)
cat("extract_ms1_features() RT clustering verified (no cross-event merging): PASS\n")

## ---- read_ms_file(): fallback path parity with parse_mzml() -------------
# Exercised unconditionally (no Spectra/mzR needed) -- prefer_spectra=FALSE
# forces read_ms_file() down the same parse_mzml() code path directly, so
# this just confirms the dispatcher wraps it faithfully (same shape, same
# content) rather than silently transforming anything.
cat("\n--- read_ms_file() fallback path ---\n")
fixture <- file.path(.pkg_root, "inst", "extdata", "batch_example", "ctrl_1.mzML")
ref <- parse_mzml(fixture)
fb  <- read_ms_file(fixture, prefer_spectra = FALSE)
stopifnot(identical(ref$ms1, fb$ms1))
stopifnot(identical(ref$ms2, fb$ms2))
stopifnot(identical(ref$info$n_ms1, fb$info$n_ms1))
stopifnot(identical(ref$info$n_ms2, fb$info$n_ms2))
cat("read_ms_file(prefer_spectra=FALSE) matches parse_mzml() exactly: PASS\n")

## ---- read_ms_file(): Spectra/mzR path, only if installed -----------------
# Same fixture, real Spectra::Spectra()/MsBackendMzR() read this time --
# skipped (not failed) when Spectra/mzR aren't available, so CI without
# Bioconductor still exercises the fallback above and stays green.
if (.have_spectra()) {
  cat("\n--- read_ms_file() Spectra/mzR path ---\n")
  sp <- read_ms_file(fixture, prefer_spectra = TRUE)
  stopifnot(nrow(sp$ms1) == nrow(ref$ms1))
  stopifnot(nrow(sp$ms2) == nrow(ref$ms2))
  stopifnot(max(abs(sort(sp$ms1$mz) - sort(ref$ms1$mz))) < 1e-6)
  stopifnot(max(abs(sort(sp$ms1$intensity) - sort(ref$ms1$intensity))) < 1e-6)
  stopifnot(all(sort(round(sp$ms2$mz, 4)) == sort(round(ref$ms2$mz, 4))))
  stopifnot(unique(sp$ms2$precursor_z) == unique(ref$ms2$precursor_z))
  cat("Spectra/mzR peaks agree with parse_mzml() within float tolerance: PASS\n")
} else {
  cat("\n--- read_ms_file() Spectra/mzR path: SKIPPED (Spectra/mzR not installed) ---\n")
}

## ---- .isotope_fit() penalizes partial isotope-pattern matches ------------
cat("\n--- .isotope_fit() partial-match penalty ---\n")
theo_iso_test <- data.frame(mz = c(1000.00, 1000.50, 1001.00, 1001.50, 1002.00),
                             abundance = c(1.0, 0.6, 0.3, 0.1, 0.05))
# Full match: every theoretical peak has a corresponding observed feature.
obs_full <- data.frame(mz = theo_iso_test$mz,
                        max_intensity = c(1.0, 0.6, 0.3, 0.1, 0.05) * 1e5)
fit_full <- .isotope_fit(theo_iso_test, obs_full, ppm_tol = 10)
# Partial match: only the first 2 of 5 theoretical peaks are observed.
obs_partial <- obs_full[1:2, ]
fit_partial <- .isotope_fit(theo_iso_test, obs_partial, ppm_tol = 10)
cat(sprintf("  full match fit: %.4f   2/5 partial match fit: %.4f\n", fit_full, fit_partial))
stopifnot(fit_full > 0.99)                 # full match still ~1
stopifnot(fit_partial < 0.5)               # 2/5 matched must score well below a naive 0.857
stopifnot(fit_partial < fit_full)
cat(".isotope_fit() partial-match penalty verified: PASS\n")

## ---- MS1-gating invariant: MS2 lookup only for MS1-confirmed hits --------
# annotate_metabolites() must never call find_ms2_spectra() for a
# metabolite that had no MS1 match -- locks in the "match MS1 first, only
# build MS2 for detected precursors" invariant the mirror-plot workflow
# relies on. mets[[9]] has no MS1 feature at all (ms1_features above only
# covers mets[1:8]), so even with a plausible-looking MS2 spectrum sitting
# right at its z=6 precursor mz, it must never produce an ms2_results entry.
cat("\n--- MS1-gating: MS2 lookup skipped without an MS1 match ---\n")
ungated_met <- mets[[9]]
ungated_mz <- (metabolite_mass_info(ungated_met)$mono_mass - 6 * .PROTON) / 6
ungated_ms2 <- data.frame(rt = 15.0, precursor_mz = ungated_mz, precursor_z = 6L,
                           mz = c(200.1, 300.2), intensity = c(500, 400),
                           stringsAsFactors = FALSE)
ms2_data_gating <- rbind(ms2_data, ungated_ms2)
results_gating <- annotate_metabolites(mets, ms1_features, ms2_data_gating,
                                        ppm_tol = 10, z_range = 3:12,
                                        adducts = c("H"), max_oxid = 0,
                                        frag_tol_ppm = 10, frag_z_range = 1:2,
                                        n_iso = 0, use_envipat = FALSE)
stopifnot(!(ungated_met$id %in% names(results_gating$ms2_results)))
stopifnot(!(ungated_met$id %in% results_gating$ms1_matches$met_id))
cat("mets[[9]] has no MS1 match and correctly gets no MS2 lookup: PASS\n")

cat("\n==== All MS matching tests passed ====\n")
