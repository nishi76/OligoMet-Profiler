# test_pmo_metabolism.R -- PMO metabolism mode: a morpholino oligo is not a
# nuclease substrate, so its metabolite set must drop the exo/endonuclease
# semantics and PS->PO/cyclic-phosphate routes, and add the phosphorodiamidate
# -> phosphoramidate charge-conversion mass variant.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (m in c("about", "progress_utils", "chemistry_dict", "oligo_io",
            "metabolites", "mass_isotope", "fragments")) {
  source(file.path(.pkg_root, "R", paste0(m, ".R")))
}
check <- function(label, ok) {
  cat(sprintf("%-74s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== PMO metabolism suite ====\n\n")

# A 6-mer PMO (uniform morpholino sugar, phosphorodiamidate linkage).
pmo <- parse_input("AMo-pdCMo-pdGMo-pdTMo-pdAMo-pdCMo")
dna <- parse_input("Ad-sCd-sGd-sTd-sAd-sCd")  # a PS DNA oligo, for contrast

check("is_pmo() detects a uniform morpholino backbone", is_pmo(pmo))
check("is_pmo() is FALSE for a DNA/PS oligo", !is_pmo(dna))

mp <- generate_metabolites(pmo, list(oligo_name = "PMO", max_3p = 3, max_5p = 3, endo = FALSE))
md <- generate_metabolites(dna, list(oligo_name = "DNA", max_3p = 3, max_5p = 3, endo = FALSE))

mech <- vapply(mp, function(m) m$mechanism %||% "", "")
kinds <- vapply(mp, function(m) m$kind, "")
check("PMO parent carries mechanism = 'parent'", mech[kinds == "parent"] == "parent")
check("PMO chain fragments carry mechanism = 'acid_hydrolysis'",
      all(mech[kinds %in% c("exo_3p", "exo_5p")] == "acid_hydrolysis"))
check("PMO fragment descriptions are reworded off 'exonuclease'",
      !any(grepl("exonuclease", vapply(mp, function(m) m$modification %||% "", ""))))
check("a DNA oligo keeps mechanism = 'nuclease'",
      all(vapply(md, function(m) m$mechanism %||% "", "")[
        vapply(md, function(m) m$kind, "") %in% c("exo_3p","exo_5p")] == "nuclease"))

# Weights: a PMO has NO 3'>5' nuclease preference -- same-depth 3' and 5'
# fragments get equal weight; a DNA oligo keeps 3' > 5'.
w_of <- function(mets, kind, k) {
  for (m in mets) if (identical(m$kind, kind) && !is.na(m$site) && m$site == k) return(m$rel_weight)
  NA_real_
}
check("PMO: 3' and 5' fragments of the same depth get EQUAL weight (no nuclease bias)",
      isTRUE(all.equal(w_of(mp, "exo_3p", 2), w_of(mp, "exo_5p", 2))))
check("PMO: parent weight is 1 and dominates the fragments",
      w_of_parent <- {
        pw <- vapply(mp, function(m) if (identical(m$kind,"parent")) m$rel_weight else NA_real_, 0)
        fw <- vapply(mp, function(m) if (m$kind %in% c("exo_3p","exo_5p")) m$rel_weight else NA_real_, 0)
        max(pw, na.rm = TRUE) == 1 && max(fw, na.rm = TRUE) < 0.2
      })
check("DNA: 3' exo is weighted above 5' exo of the same depth (nuclease prior intact)",
      w_of(md, "exo_3p", 2) > w_of(md, "exo_5p", 2))

# No sulfur in a PMO -> no PS->PO oxidation series possible.
check("PMO metabolites carry zero phosphorothioate bonds (PS->PO auto-off)",
      all(vapply(mp, function(m) m$n_ps %||% 0, 0) == 0))

# Phosphorodiamidate -> phosphoramidate charge-conversion mass variant.
parent <- mp[[which(kinds == "parent")]]
info <- metabolite_mass_info(parent, STANDARD_DICT)
v_pd <- metabolite_mod_variants(parent, info, STANDARD_DICT, mods = "pd_hydrolysis")
pdv <- Filter(function(v) v$mod == "pd_hydrolysis", v_pd)
check("pd_hydrolysis variant is produced for a PMO", length(pdv) == 1)
check("pd_hydrolysis mass shift is -27.0473 Da (NMe2 -> OH)",
      abs((pdv[[1]]$mono_mass - info$mono_mass) - (-27.0473)) < 1e-3)
check("pd_hydrolysis is MS1-unlocalized (mod_site NA)", is.na(pdv[[1]]$mod_site))

# A DNA oligo has no 'pd' linkage, so pd_hydrolysis must NOT fire for it.
pm <- md[[which(vapply(md, function(m) m$kind, "") == "parent")]]
info_d <- metabolite_mass_info(pm, STANDARD_DICT)
v_d <- metabolite_mod_variants(pm, info_d, STANDARD_DICT, mods = c("pd_hydrolysis", "depurination"))
check("pd_hydrolysis does NOT fire on a non-PMO oligo",
      !any(vapply(v_d, function(v) v$mod, "") == "pd_hydrolysis"))
check("depurination still fires on the DNA oligo (purine present)",
      any(vapply(v_d, function(v) v$mod, "") == "depurination"))

cat("\nAll PMO-metabolism checks passed.\n")
