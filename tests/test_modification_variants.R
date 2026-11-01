# test_modification_variants.R -- in-place chemical modification variants
# (R/mass_isotope.R): deamination, 2',3'-cyclic phosphate, the consolidated
# metabolite_mod_variants() enumerator, and its wiring into match_ms1().

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R", "ms_matching.R")) {
  source(file.path(.pkg_root, "R", .f))
}
check <- function(label, ok) {
  cat(sprintf("%-82s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
approx <- function(a, b, tol = 1e-3) abs(a - b) < tol
cat("==== Modification variants suite ====\n\n")

## ---- deamination ----------------------------------------------------------
spec <- list(bases = c("A", "C", "G", "T", "U", "S"),
             sugars = rep("e", 6), linkages = c(rep("s", 5), NA),
             conj5 = "none", conj3 = "none", n = 6)
par <- generate_metabolites(spec, list(oligo_name = "x", endo = FALSE,
                                        max_3p = 0, max_5p = 0))[[1]]
dv <- deamination_variants(par)
check("deamination hits only aminated bases (A,C,G,S) not T/U",
      setequal(dv$base, c("A", "C", "G", "S")))
check("every deamination shift is +0.98402 Da",
      all(approx(dv$mass_shift, 0.9840156)))
check("cytosine deaminates to U, adenine to hypoxanthine",
      dv$product[dv$base == "C"] == "U" && dv$product[dv$base == "A"] == "hypoxanthine")

## ---- cyclic phosphate -----------------------------------------------------
# A 3'-phosphate fragment ending on ribose can cyclize; DNA/2'-mod cannot.
rna <- list(bases = c("A", "C", "G", "U"), sugars = rep("r", 4),
            linkages = c("s", "s", "s", NA), conj5 = "none", conj3 = "none", n = 4)
rmets <- generate_metabolites(rna, list(oligo_name = "r", endo = TRUE,
                                        endo_sites = 2, max_3p = 1, max_5p = 1,
                                        min_frag_len = 1))
cyc <- Filter(can_form_cyclic_phosphate, rmets)
check("ribose 3'-phosphate fragment can cyclize", length(cyc) >= 1)
check("cyclic phosphate shift is -H2O (-18.0106)",
      approx(cyclic_phosphate_variant(cyc[[1]])$mass_shift, -18.010565))
# the all-MOE parent (3'-OH, no ribose) cannot cyclize
check("2'-MOE parent (3'-OH) cannot cyclize", !can_form_cyclic_phosphate(par))

## ---- consolidated enumerator ----------------------------------------------
v_none <- metabolite_mod_variants(par, mods = character(0))
check("no mods -> only the unmodified species", length(v_none) == 1 && v_none[[1]]$mod == "none")

v_all <- metabolite_mod_variants(par, mods = "all")
mods_seen <- vapply(v_all, function(x) x$mod, "")
check("mods='all' includes depurination, oxidation, deamination",
      all(c("depurination", "oxidation", "deamination") %in% mods_seen))
# same-mass sites collapse to mod_site = NA (MS1 can't localize)
deam_rows <- Filter(function(x) x$mod == "deamination", v_all)
check("deamination collapses to one unresolved row (mod_site = NA)",
      length(deam_rows) == 1 && is.na(deam_rows[[1]]$mod_site))
# depurination of A vs G differ by one O, so two distinct masses survive
dep_rows <- Filter(function(x) x$mod == "depurination", v_all)
check("depurination keeps distinct A-loss and G-loss masses",
      length(dep_rows) == 2)

## ---- wiring into match_ms1 ------------------------------------------------
info <- metabolite_mass_info(par)
prot <- 1.007276466
z <- 3
mk <- function(m) (m - z * prot) / z
feats <- data.frame(
  mz = c(mk(info$mono_mass), mk(info$mono_mass + 0.9840156)),
  rt = c(5, 5), max_intensity = c(1e6, 3e5))
base <- match_ms1(list(par), feats, ppm_tol = 5, z_range = 3:4, adducts = "H", n_iso = 0)
check("unmodified search ignores the deaminated feature",
      all(base$mod == "none") && nrow(base) == 1)
withm <- match_ms1(list(par), feats, ppm_tol = 5, z_range = 3:4, adducts = "H",
                   n_iso = 0, mods = "deamination")
check("mods='deamination' finds the +0.984 feature",
      "deamination" %in% withm$mod)
check("match table carries mod and mod_site columns",
      all(c("mod", "mod_site") %in% names(withm)))

cat("\nAll modification-variant checks passed.\n")
