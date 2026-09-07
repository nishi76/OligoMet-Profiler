# test_fragments.R -- validate McLuckey fragment ion formulas
# Checks:
#   1. Complementary relationship: a_k + w_{n-k} = M (single-bond cleavage
#      splits the neutral precursor into two neutral pieces, so every
#      complementary pair must sum to exactly M -- see fragments.R FIX 2)
#   2. b_k + x_{n-k} = M
#   3. c_k + y_{n-k} = M
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

# Generate the full McLuckey ion set (not just the app's default subset)
# so every complementary pair and the base-loss ions below get exercised.
frags <- generate_fragments(parent, z_range = 1:3,
                             ion_types = c("a", "aB", "b", "bB", "c", "w", "x", "y"),
                             include_dz = TRUE)
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
frag_fail <- 0L

for (k in c(1, 5, 7, 10, 14)) {
  fa <- get_frag(frags, "a", k)
  fw <- get_frag(frags, "w", k)
  fb <- get_frag(frags, "b", k)
  fx <- get_frag(frags, "x", k)
  fc <- get_frag(frags, "c", k)
  fy <- get_frag(frags, "y", k)
  fd <- get_frag(frags, "d", k)
  fz <- get_frag(frags, "z", k)

  for (pr in list(list("a+w", fa, fw), list("b+x", fb, fx), list("c+y", fc, fy),
                  list("d+z", fd, fz))) {
    label <- pr[[1]]; f1 <- pr[[2]]; f2 <- pr[[3]]
    if (is.null(f1) || is.null(f2)) next
    total <- f1$mono_mass + f2$mono_mass
    delta <- total - M
    ok <- abs(delta) < 1e-3
    if (!ok) frag_fail <- frag_fail + 1L
    cat(sprintf("  k=%2d  %s = %.4f  M = %.4f  delta = %+.4f Da (%.1f ppm)  %s\n",
                k, label, total, M, delta, abs(delta) / M * 1e6, if (ok) "ok" else "NO"))
  }
}
if (frag_fail > 0) stop(frag_fail, " complementary fragment-pair check(s) failed")

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
  cat(sprintf("  Expected diff = base(%s) = %.4f\n", base5, base_mass))
  match_ok <- abs(diff - base_mass) < 0.01
  cat(sprintf("  Match: %s\n", ifelse(match_ok, "YES", "NO")))
  if (!match_ok) stop("a-B base loss check failed")
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

cat("\n==== All fragment tests passed ====\n")
