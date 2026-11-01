# test_inclusion_prioritization.R -- gapmer-aware endo default, relevance
# tiers, mass-scaled charge window, and tier/cap filtering of the inclusion
# list (R/metabolites.R, R/mass_isotope.R, R/fragments.R).

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R")) {
  source(file.path(.pkg_root, "R", .f))
}
check <- function(label, ok) {
  cat(sprintf("%-78s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
cat("==== Inclusion-list prioritisation suite ====\n\n")

gap <- list(bases = strsplit("TCGTCATGACTTGCAGAG", "")[[1]],
            sugars = c(rep("e", 5), rep("d", 8), rep("e", 5)),
            linkages = c(rep("s", 17), NA), conj5 = "none", conj3 = "none", n = 18)
sso <- list(bases = gap$bases, sugars = rep("e", 18),
            linkages = c(rep("s", 17), NA), conj5 = "none", conj3 = "none", n = 18)
dna <- list(bases = gap$bases, sugars = rep("d", 18),
            linkages = c(rep("s", 17), NA), conj5 = "none", conj3 = "none", n = 18)

## ---- gapmer detection --------------------------------------------------------
check("5-10-5 MOE/DNA is a gapmer", is_gapmer(gap))
check("uniform 2'-MOE SSO is not a gapmer", !is_gapmer(sso))
check("uniform DNA (first-gen PS) is not a gapmer", !is_gapmer(dna))

## ---- auto endo default trims the gapmer library ------------------------------
n_auto <- length(generate_metabolites(gap, list(oligo_name = "g")))             # auto -> gap, mf 5
n_all  <- length(generate_metabolites(gap, list(oligo_name = "g",
                                               endo_sites = "all", min_frag_len = 3)))
check("auto (gap) + min_frag 5 yields fewer species than all + min_frag 3",
      n_auto < n_all)
# a non-gapmer under auto still searches all internal bonds
check("auto on an SSO searches all positions (no DNA gap to restrict to)",
      length(generate_metabolites(sso, list(oligo_name = "s"))) ==
        length(generate_metabolites(sso, list(oligo_name = "s", endo_sites = "all"))))

## ---- the opts$endo partial-match guard ---------------------------------------
# $ partial-matches endo -> endo_sites; [[ (used internally) must not.
check("passing endo_sites without endo does not crash (exact-match guard)",
      is.list(generate_metabolites(gap, list(oligo_name = "g", endo_sites = "all"))))

## ---- relevance tiers ---------------------------------------------------------
tb <- metabolite_table(generate_metabolites(gap, list(oligo_name = "g")))
check("metabolite table carries a tier column", "tier" %in% names(tb))
check("parent and 3'-exo are Tier 1", all(tb$tier[tb$kind %in% c("parent", "exo_3p")] == 1))
check("endonuclease fragments are Tier 2", all(tb$tier[grepl("^endo", tb$kind)] == 2))

## ---- mass-scaled charge window -----------------------------------------------
z7 <- suggest_z_range(7122)
check("a ~7 kDa oligo suggests a higher charge floor than 3", z7[1] >= 4)
check("suggested window widens with mass", suggest_z_range(14000)[2] > z7[2])
check("tiny mass falls back sanely", all(suggest_z_range(0) == c(3L, 12L)))

## ---- tier / cap filtering of the inclusion list ------------------------------
mets <- generate_metabolites(gap, list(oligo_name = "g"))
all_prm <- prm_inclusion_list(mets, z_range = 4:11, max_oxid = 4)
t1_prm  <- prm_inclusion_list(mets, z_range = 4:11, max_oxid = 4, scope = "tier1")
cap_prm <- prm_inclusion_list(mets, z_range = 4:11, max_oxid = 4, max_targets = 50)
check("tier1 scope keeps only Tier-1 rows and is smaller", all(t1_prm$tier == 1) && nrow(t1_prm) < nrow(all_prm))
check("max_targets caps the row count", nrow(cap_prm) == 50)
check("the cap keeps Tier-1 rows first", all(cap_prm$tier[order(cap_prm$tier)][1] == 1))

## ---- nuclease-preference weighting (item 3) ----------------------------------
wtb <- metabolite_table(generate_metabolites(gap, list(oligo_name = "g")))
check("metabolite table carries a rel_weight column", "rel_weight" %in% names(wtb))
check("parent rel_weight is 1", wtb$rel_weight[wtb$kind == "parent"] == 1)
check("3'-exo outranks 5'-exo (dominant plasma pathway)",
      max(wtb$rel_weight[wtb$kind == "exo_3p"]) > max(wtb$rel_weight[wtb$kind == "exo_5p"]))
check("3'-exo weight decays with truncation depth (N-1 > N-2)",
      wtb$rel_weight[wtb$kind == "exo_3p" & wtb$site == "1"] >
        wtb$rel_weight[wtb$kind == "exo_3p" & wtb$site == "2"])
check("endonuclease weight does not exceed the dominant 3'-exo route",
      max(wtb$rel_weight[grepl("^endo", wtb$kind)]) <= max(wtb$rel_weight[wtb$kind == "exo_3p"]))
# pyrimidine terminus is removed faster than a purine terminus at equal depth
seqAC <- list(bases = c("A","A","A","A","C"), sugars = rep("d", 5),
              linkages = c(rep("o", 4), NA), conj5 = "none", conj3 = "none", n = 5)
wAC <- assign_nuclease_weights(list(
  .make_met(NA, "pyr", "exo_3p", "", 1L, truncate_3p(seqAC, 1)),   # removed C (pyrimidine)
  .make_met(NA, "pur", "exo_5p", "", 1L, truncate_5p(seqAC, 1))),  # removed A (purine)
  seqAC)
check("removing a pyrimidine 3' terminus scores above removing a purine 5' terminus",
      wAC[[1]]$rel_weight > wAC[[2]]$rel_weight)

## ---- GalNAc conjugate catabolism (item 4) ------------------------------------
sirna <- list(bases = strsplit("UCGUCAUGACUUGCAGAG", "")[[1]], sugars = rep("f", 18),
              linkages = c(rep("o", 17), NA), conj5 = "none", conj3 = "GalNAc3", n = 18)
check("GalNAc3-conjugated strand is detected", is_galnac_conjugated(sirna))
check("a bare oligo is not GalNAc-conjugated", !is_galnac_conjugated(sso))
catab <- generate_conjugate_catabolites(sirna, "si")
check("GalNAc3 yields a 4-step ladder (-1/-2/-3 sugar + de-conjugation)", length(catab) == 4)
cmass <- vapply(catab, function(x) metabolite_mass_info(x)$mono_mass, numeric(1))
intact_mass <- metabolite_mass_info(.make_met(NA, "p", "parent", "", NA, sirna))$mono_mass
steps <- diff(c(intact_mass, cmass))
check("each sugar-loss step drops one GalNAc residue (-203.0794 Da)",
      all(abs(steps[1:3] - (-203.0794)) < 1e-3))
gtab <- metabolite_table(generate_metabolites(sirna, list(oligo_name = "si")))
check("catabolites enter the library and are Tier 1 (early-dominant)",
      sum(gtab$kind == "conj_catab") == 4 && all(gtab$tier[gtab$kind == "conj_catab"] == 1))
check("conj_catab off suppresses the ladder",
      sum(metabolite_table(generate_metabolites(sirna,
        list(oligo_name = "si", conj_catab = FALSE)))$kind == "conj_catab") == 0)

cat("\nAll inclusion-prioritisation checks passed.\n")
