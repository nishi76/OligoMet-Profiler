# test_mass_isotope.R -- validate mass_isotope.R + reverse-engineer workbook offset
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

sp <- parse_input(ION337_TRIPLET)
mets <- generate_metabolites(sp, list(oligo_name = "ION337", endo = FALSE))
par <- mets[[1]]
info <- metabolite_mass_info(par)

cat("=== 1. Parent mass info ===\n")
cat("formula:", info$formula_str, "  mono:", sprintf("%.6f", info$mono_mass), "\n")

# 2. Standard envelope
cat("\n=== 2. Standard envelope [M-zH]^z- (h_offset=0), z=5..12 ===\n")
env <- charge_envelope(info$mono_mass, 5:12, h_offset = 0)
print(env, digits = 4)

# 3. Reverse-engineer workbook offset
cat("\n=== 3. Reverse-engineer ION337 workbook offset ===\n")
library(openxlsx)
REFERENCE_WORKBOOK <- Sys.getenv("ION337_WORKBOOK_PATH", file.path(.pkg_root, "tests", "ION337_Workbook.xlsx"))
if (!file.exists(REFERENCE_WORKBOOK)) {
  stop("Reference workbook not found at ", REFERENCE_WORKBOOK, ". ",
       "Place ION337_Workbook.xlsx in tests/, or set the ",
       "ION337_WORKBOOK_PATH environment variable to its location. ",
       "This file is not distributed with the package.")
}
wb <- readWorkbook(REFERENCE_WORKBOOK, sheet = 1)
# Parent block: rows 3,6,9,12,15,18,21,24,27,30 have charge + m/z
# charge in Sequence col, m/z in Sequence.Triplet col
wb_z <- as.numeric(wb$Sequence[c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30)])
wb_mz <- as.numeric(wb$Sequence.Triplet[c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30)])
wb_z <- abs(wb_z)  # charge states stored as negative
cat(sprintf("%-4s %-14s %-14s %-10s\n", "z", "wb_mz", "offset", "our_mz(0)"))
offsets <- c()
for (i in seq_along(wb_z)) {
  z <- wb_z[i]; mz_wb <- wb_mz[i]
  # offset = mz*z + z*mp - M
  off <- mz_wb * z + z * .PROTON - info$mono_mass
  our0 <- charge_envelope(info$mono_mass, z, h_offset = 0)$mz[1]
  cat(sprintf("%-4d %-14.6f %-14.6f %-10.6f\n", z, mz_wb, off, our0))
  offsets <- c(offsets, off)
}
cat(sprintf("\nMean offset: %.6f  SD: %.6f\n", mean(offsets), sd(offsets)))
cat(sprintf("3*proton = %.6f   3*H_atom = %.6f\n", 3 * .PROTON, 3 * .atomic_mass_mono[["H"]]))

# 4. Verify with empirical offset
emp_off <- mean(offsets)
cat(sprintf("\n=== 4. Envelope with empirical offset %.4f ===\n", emp_off))
env2 <- charge_envelope(info$mono_mass, wb_z, h_offset = emp_off)
cmp <- data.frame(z = wb_z, wb_mz = wb_mz, our_mz = env2$mz,
                  ppm = abs(env2$mz - wb_mz) / wb_mz * 1e6)
print(cmp, digits = 5)
cat("Max ppm:", sprintf("%.4f\n", max(cmp$ppm)))

# 5. Isotope pattern (enviPat)
cat("\n=== 5. Isotope pattern (parent, enviPat, top 8) ===\n")
pat <- isotope_pattern(info$formula_str, n_top = 8)
print(pat, digits = 4)

# 6. Isotope m/z at z=12 vs workbook
cat("\n=== 6. Isotope m/z cluster at z=12 vs workbook ===\n")
cl <- isotope_mz_cluster(info$formula_vec, z = 12, n_top = 5, h_offset = emp_off)
cat("Workbook isotopes (z=12): 515.33, 515.413, 515.497, 515.247, 515.58\n")
cat("Our top-5 by abundance:\n")
print(cl[, c("iso", "mz", "abundance")], digits = 4)

# 7. PS oxidation series
cat("\n=== 7. PS->PO oxidation series (k=0..6) ===\n")
ox <- ps_oxidation_series(par, max_oxid = 6)
print(ox, digits = 4)

# 8. Adducts + depurination
cat("\n=== 8. Adducts ===\n")
for (a in c("Na", "K", "NH4")) cat(sprintf("  +%s: %+.4f Da\n", a, adduct_shift(a)))
cat("\n=== 9. Depurination variants ===\n")
dep <- depurination_variants(par)
if (!is.null(dep)) print(dep, digits = 4, row.names = FALSE)
