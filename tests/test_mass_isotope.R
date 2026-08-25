# test_mass_isotope.R -- validate mass_isotope.R against published formulas
#
# Anchors: nusinersen and inotersen, whose free-acid molecular formulas are
# on their approved product labelling. Prints computed values for reading;
# the formula comparison is the one hard pass/fail here.
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

cat("=== 1. Formula engine vs published molecular formulas ===\n")
ref <- validate_reference(verbose = TRUE)
if (!isTRUE(ref$ok)) stop("Formula engine no longer reproduces the published formulas")

sp <- parse_input(INOTERSEN_TRIPLET)
mets <- generate_metabolites(sp, list(oligo_name = "inotersen", endo = FALSE))
par <- mets[[1]]
info <- metabolite_mass_info(par)

cat("=== 2. Parent mass info (inotersen) ===\n")
cat("formula:", info$formula_str,
    "  mono:", sprintf("%.6f", info$mono_mass),
    "  avg:", sprintf("%.4f", info$avg_mass), "\n")

cat("\n=== 3. Standard envelope [M-zH]^z- (h_offset=0), z=5..12 ===\n")
print(charge_envelope(info$mono_mass, 5:12, h_offset = 0), digits = 6)

cat("\n=== 4. Isotope pattern (parent, top 8) ===\n")
print(isotope_pattern(info$formula_str, n_top = 8), digits = 6)

# Two invariants that hold for any correct isotope engine:
#   (a) no isotopologue is lighter than the monoisotopic mass, and
#   (b) the monoisotopic peak is in the pattern when enough peaks are asked
#       for. isotope_pattern() returns the top-N *by abundance*, and for a
#       7 kDa oligo the monoisotopic peak is not the tallest -- at n_top=8 it
#       legitimately falls outside the set, so this asks for a generous 60.
# Together these catch an engine that raises the per-atom pattern to the
# wrong convolution power: the built-in fallback once returned ~15571 Da for
# this 7178 Da molecule, and every m/z derived from it (envelope sheet,
# isotope plot, MS1 isotope fit, spectral libraries) was wrong on any
# machine without enviPat installed.
.iso_check <- function(formula_vec, mono, label) {
  for (engine in c(TRUE, FALSE)) {
    p <- isotope_pattern(formula_vec, n_top = 60, use_envipat = engine)
    if (is.null(p) || nrow(p) == 0) next
    closest <- min(abs(p$mass - mono))
    cat(sprintf("  %-12s use_envipat=%-5s lightest %.4f, mono %.4f, closest peak %+.4f Da\n",
                label, engine, min(p$mass), mono, closest))
    if (min(p$mass) < mono - 0.05)
      stop("Isotope pattern contains a peak lighter than the monoisotopic ",
           "mass (", label, ", use_envipat=", engine, ")")
    if (closest > 0.05)
      stop("Isotope pattern does not contain the monoisotopic peak (",
           label, ", use_envipat=", engine, ")")
  }
}
cat("\n=== 4b. Isotope pattern anchors on the monoisotopic mass ===\n")
.iso_check(info$formula_vec, info$mono_mass, "inotersen")
for (fs in c("C8H10N4O2", "C6H12O6", "S8")) {
  .iso_check(parse_formula(fs), formula_mass(parse_formula(fs)), fs)
}

cat("\n=== 5. Isotope m/z cluster at z=8 ===\n")
cl <- isotope_mz_cluster(info$formula_vec, z = 8, n_top = 5)
print(cl[, c("iso", "mz", "abundance")], digits = 6)

cat("\n=== 6. Charge-state consistency check ===\n")
# Every charge state must back-calculate to the same neutral mass.
zs <- 4:12
neutral <- vapply(zs, function(z) {
  mz <- charge_envelope(info$mono_mass, z)$mz[1]
  mz * z + z * .PROTON
}, numeric(1))
cat("max deviation across z=4..12:",
    sprintf("%.3e Da\n", max(abs(neutral - info$mono_mass))))
if (max(abs(neutral - info$mono_mass)) > 1e-6) {
  stop("Charge envelope is not self-consistent across charge states")
}

cat("\n=== 7. PS->PO oxidation series (k=0..6) ===\n")
ox <- ps_oxidation_series(par, max_oxid = 6)
print(ox, digits = 8)
# Each desulfurization must shift the mass by exactly one S->O swap.
deltas <- diff(ox$mono_mass)
cat("observed per-event shift:", sprintf("%.6f Da", unique(round(deltas, 6))),
    " expected:", sprintf("%.6f Da\n", -.PS_TO_PO_SHIFT))

cat("\n=== 8. Adducts ===\n")
for (a in c("Na", "K", "NH4")) cat(sprintf("  +%s: %+.4f Da\n", a, adduct_shift(a)))

cat("\n=== 9. Depurination variants ===\n")
dep <- depurination_variants(par)
if (!is.null(dep)) print(head(dep, 6), digits = 6, row.names = FALSE)

cat("\n==== mass_isotope tests complete ====\n")
