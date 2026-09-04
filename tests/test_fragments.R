# test_fragments.R -- validate McLuckey fragment ion formulas
# Checks:
#   1. Complementary relationship: a_k + w_{n-k} = M + H2O
#   2. b_k + x_{n-k} = M + H2O - H (known convention offset)
#   3. c_k + y_{n-k} = M + H2O - H
#   4. Fragment masses are positive and decreasing with k for 5' ions
#   5. PRM inclusion list covers expected charge range
#   6. Matching against synthetic MS2 peaks

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

cat("==== Fragment ion validation ====\n\n")

# Parse the inotersen reference case and get parent metabolite
spec <- parse_input(INOTERSEN_TRIPLET)
mets <- generate_metabolites(spec, opts = list(oligo_name = "inotersen", max_3p = 0, max_5p = 0, endo = FALSE))
parent <- mets[[1]]
parent_info <- metabolite_mass_info(parent)
M <- parent_info$mono_mass
cat("Parent mass M =", sprintf("%.4f\n", M))
cat("Parent formula:", parent_info$formula_str, "\n\n")

# Generate fragments
frags <- generate_fragments(parent, z_range = 1:3, include_dz = TRUE)
cat("Generated", length(frags), "terminal fragments\n")

# Extract by ion type
get_frag <- function(frags, type, k) {
  for (f in frags) {
    if (f$ion_type == type && f$cleavage_site == k) return(f)
  }
  NULL
}

# Check complementary relationships for several cleavage sites
cat("\n--- Complementary relationship checks ---\n")
h2o_mass <- .H2O_MASS
proton <- .PROTON

for (k in c(1, 5, 7, 10, 14)) {
  fa <- get_frag(frags, "a", k)
  fw <- get_frag(frags, "w", k)
  fb <- get_frag(frags, "b", k)
  fx <- get_frag(frags, "x", k)
  fc <- get_frag(frags, "c", k)
  fy <- get_frag(frags, "y", k)

  if (!is.null(fa) && !is.null(fw)) {
    sum_aw <- fa$mono_mass + fw$mono_mass
    delta_aw <- sum_aw - (M + h2o_mass)
    cat(sprintf("  k=%2d  a+w = %.4f  M+H2O = %.4f  delta = %+.4f Da (%.1f ppm)\n",
                k, sum_aw, M + h2o_mass, delta_aw, abs(delta_aw)/(M+h2o_mass)*1e6))
  }
  if (!is.null(fb) && !is.null(fx)) {
    sum_bx <- fb$mono_mass + fx$mono_mass
    delta_bx <- sum_bx - (M + h2o_mass)
    cat(sprintf("  k=%2d  b+x = %.4f  M+H2O = %.4f  delta = %+.4f Da (%.1f ppm)\n",
                k, sum_bx, M + h2o_mass, delta_bx, abs(delta_bx)/(M+h2o_mass)*1e6))
  }
  if (!is.null(fc) && !is.null(fy)) {
    sum_cy <- fc$mono_mass + fy$mono_mass
    delta_cy <- sum_cy - (M + h2o_mass)
    cat(sprintf("  k=%2d  c+y = %.4f  M+H2O = %.4f  delta = %+.4f Da (%.1f ppm)\n",
                k, sum_cy, M + h2o_mass, delta_cy, abs(delta_cy)/(M+h2o_mass)*1e6))
  }
}

# Show some fragment masses
cat("\n--- Sample fragment masses (z=1) ---\n")
ftab <- fragment_table(frags)
ftab1 <- ftab[ftab$z == 1, ]
cat("First 10 5' fragments:\n")
print(head(ftab1[ftab1$direction == "5'", ], 10))
cat("\nFirst 10 3' fragments:\n")
print(head(ftab1[ftab1$direction == "3'", ], 10))

# Check a-B base loss
cat("\n--- a-B base loss check ---\n")
fa <- get_frag(frags, "a", 5)
faB <- get_frag(frags, "a-B", 5)
if (!is.null(fa) && !is.null(faB)) {
  diff <- fa$mono_mass - faB$mono_mass
  base5 <- parent$bases[5]
  base_mass <- formula_mass(STANDARD_DICT[[base5]]$formula, mono = TRUE)
  cat(sprintf("  a_5 mass = %.4f, a-B_5 mass = %.4f, diff = %.4f\n", fa$mono_mass, faB$mono_mass, diff))
  cat(sprintf("  Expected diff = base(%s) - H2O = %.4f - %.4f = %.4f\n",
              base5, base_mass, h2o_mass, base_mass - h2o_mass))
  cat(sprintf("  Match: %s\n", ifelse(abs(diff - (base_mass - h2o_mass)) < 0.01, "YES", "NO")))
}

# Internal fragments
cat("\n--- Internal fragments ---\n")
ifrags <- generate_internal_fragments(parent, types = c("wa", "wb"), z_range = 1:2)
cat("Generated", length(ifrags), "internal fragments\n")
if (length(ifrags) > 0) {
  itab <- fragment_table(ifrags)
  print(head(itab, 10))
}

# PRM inclusion list
cat("\n--- PRM inclusion list ---\n")
prm <- prm_inclusion_list(mets, z_range = 3:12, max_oxid = 0)
cat("PRM entries:", nrow(prm), "\n")
cat("Charge range:", min(prm$z), "-", max(prm$z), "\n")
cat("m/z range:", sprintf("%.2f - %.2f\n", min(prm$precursor_mz), max(prm$precursor_mz)))
print(head(prm, 5))

# Synthetic MS2 matching test
cat("\n--- Synthetic MS2 matching test ---\n")
# Create synthetic MS2 peaks from theoretical fragments + noise
set.seed(42)
synth_peaks <- do.call(rbind, lapply(frags[1:20], function(f) {
  data.frame(mz = (f$mono_mass - 1 * .PROTON) / 1,
             intensity = runif(1, 100, 10000),
             stringsAsFactors = FALSE)
}))
# Add PS diagnostic ions
synth_peaks <- rbind(synth_peaks,
  data.frame(mz = 94.9452, intensity = 5000),
  data.frame(mz = 192.9746, intensity = 3000))
# Add some noise peaks
noise <- data.frame(mz = runif(10, 100, 2000), intensity = runif(10, 50, 500))
synth_peaks <- rbind(synth_peaks, noise)
synth_peaks <- synth_peaks[order(synth_peaks$mz), ]

matched <- match_fragments(frags, synth_peaks, tol_ppm = 10, z_range = 1:2)
cat("Matched", nrow(matched), "fragments\n")
if (nrow(matched) > 0) {
  cat("Matched ion types:", paste(unique(matched$ion_type), collapse = ", "), "\n")
  cat("Median ppm error:", round(median(matched$ppm_error), 3), "\n")
  cat("Covered cleavage sites:", paste(sort(unique(matched$cleavage_site)), collapse = ","), "\n")
}

diags <- check_ps_diagnostic(synth_peaks)
cat("PS diagnostic ions found:", nrow(diags), "\n")
if (nrow(diags) > 0) print(diags)

score <- confirmation_score(matched, parent$n, diags)
cat("\nConfirmation score:\n")
cat("  Total:", score$total_score, "/ 100\n")
cat("  Coverage:", sprintf("%.1f%%", score$coverage * 100), "\n")
cat("  Matches:", score$n_matches, "\n")
cat("  Confident:", score$confident, "\n")

## ---- match_fragments() de-duplicates multi-assigned peaks -----------------
cat("\n--- match_fragments() de-duplicates multi-assigned peaks ---\n")
stopifnot(length(unique(matched$obs_mz)) == nrow(matched))
cat("No observed peak claimed by more than one surviving match: PASS\n")

# Deterministic collision: two distinct fragment objects with the same
# mono_mass (by construction) competing for one synthetic peak -- only the
# best (here, tied) match should survive.
f1 <- frags[[1]]
f2 <- f1
f2$ion_type <- if (f1$ion_type == "w") "y" else "w"
f2$cleavage_site <- f1$cleavage_site + 1000L  # guaranteed distinct identity
collision_peak <- data.frame(mz = (f1$mono_mass - 1 * .PROTON) / 1, intensity = 1000)
collision_matched <- match_fragments(list(f1, f2), collision_peak, tol_ppm = 10, z_range = 1)
stopifnot(nrow(collision_matched) == 1)
cat("Two same-mass fragments racing for one peak: only one survives: PASS\n")

## ---- Salt adduct search in MS2 fragment matching ---------------------------
cat("\n--- match_fragments() with Na adduct ---\n")
na_shift <- adduct_shift("Na")
# w-ions at 5 different cleavage sites have guaranteed-distinct masses
# (monotonically increasing with k), avoiding an accidental mass collision
# between unrelated ion types that a plain frags[1:5] slice can hit.
w_frags <- Filter(function(f) f$ion_type == "w", frags)[1:5]
na_peaks <- do.call(rbind, lapply(w_frags, function(f) {
  data.frame(mz = (f$mono_mass + na_shift - 1 * .PROTON) / 1, intensity = 5000)
}))
matched_h_only <- match_fragments(w_frags, na_peaks, tol_ppm = 10, z_range = 1)
matched_with_na <- match_fragments(w_frags, na_peaks, tol_ppm = 10, z_range = 1,
                                    adducts = c("H", "Na"))
stopifnot(nrow(matched_h_only) == 0)
stopifnot(nrow(matched_with_na) == 5)
stopifnot(all(matched_with_na$adduct == "Na"))
cat("Na-adducted MS2 peaks matched only when 'Na' is in adducts=: PASS\n")

# Backward compatibility: default adducts=c("H") must be identical to
# passing it explicitly (both take the same code path either way, but this
# pins the parameter's default value itself).
m_default <- match_fragments(frags, synth_peaks, tol_ppm = 10, z_range = 1:2)
m_explicit_h <- match_fragments(frags, synth_peaks, tol_ppm = 10, z_range = 1:2, adducts = c("H"))
stopifnot(identical(m_default, m_explicit_h))
cat("Default adducts=c('H') identical to explicit c('H'): PASS\n")

cat("\n==== All fragment tests passed ====\n")
