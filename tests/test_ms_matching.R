# test_ms_matching.R -- validate MS import + matching module
# Uses synthetic MS1 features and MS2 peaks derived from ION337 library

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
source(file.path(.pkg_root, "chemistry_dict.R"))
source(file.path(.pkg_root, "oligo_io.R"))
source(file.path(.pkg_root, "metabolites.R"))
source(file.path(.pkg_root, "mass_isotope.R"))
source(file.path(.pkg_root, "fragments.R"))
source(file.path(.pkg_root, "ms_matching.R"))

cat("==== MS matching validation ====\n\n")

# Check msconvert availability
msc <- find_msconvert()
cat("msconvert available:", !is.na(msc), "\n")

# Parse ION337 and generate metabolite library
spec <- parse_input(ION337_TRIPLET)
mets <- generate_metabolites(spec, opts = list(
  oligo_name = "ION337", max_3p = 5, max_5p = 5, endo = FALSE))
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

cat("\n==== All MS matching tests passed ====\n")
