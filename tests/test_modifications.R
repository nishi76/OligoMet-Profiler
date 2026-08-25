# test_modifications.R -- anchor every sugar and linkage residue
#
# The residue formulas in chemistry_dict.R are only meaningful if they
# reassemble into the right molecule. Two checks:
#
#   1. Sugars. nucleoside = base + sugar - H2O. Each sugar below is
#      checked against a nucleoside whose formula is known independently
#      (deoxyadenosine, adenosine, 2',3'-seco-uridine, and so on), so a
#      typo in a residue formula fails here rather than silently shifting
#      every mass in a run.
#   2. Linkages. Every backbone chemistry is the phosphodiester bridge
#      with its non-bridging -OH swapped for something else, so each one
#      is checked as an explicit mass difference from PO.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (m in c("progress_utils", "chemistry_dict", "oligo_io", "metabolites",
            "mass_isotope")) {
  source(file.path(.pkg_root, "R", paste0(m, ".R")))
}

fail <- 0L

cat("=== 1. Sugar residues rebuild a known nucleoside ===\n")
cat(sprintf("  %-7s %-5s %-22s %-22s %s\n",
            "sugar", "base", "computed", "expected", "ok"))

# base + sugar - H2O, compared against the published nucleoside formula.
nucleoside_anchors <- list(
  list(sugar = "d",     base = "A", expect = "C10H13N5O3",
       what = "2'-deoxyadenosine"),
  list(sugar = "r",     base = "A", expect = "C10H13N5O4",
       what = "adenosine"),
  list(sugar = "m",     base = "A", expect = "C11H15N5O4",
       what = "2'-O-methyladenosine"),
  list(sugar = "f",     base = "A", expect = "C10H12N5O3F",
       what = "2'-fluoro-2'-deoxyadenosine"),
  list(sugar = "MOE",   base = "A", expect = "C13H19N5O5",
       what = "2'-O-MOE-adenosine"),
  list(sugar = "UNA",   base = "U", expect = "C9H14N2O6",
       what = "2',3'-seco-uridine (UNA-U)"),
  list(sugar = "GNA",   base = "A", expect = "C8H11N5O2",
       what = "(S)-9-(2,3-dihydroxypropyl)adenine (GNA-A)"),
  list(sugar = "allyl", base = "A", expect = "C13H17N5O4",
       what = "2'-O-allyladenosine"),
  list(sugar = "NP",    base = "T", expect = "C10H15N3O4",
       what = "3'-amino-3'-deoxythymidine"),
  list(sugar = "NMA",   base = "A", expect = "C13H18N6O5",
       what = "2'-O-NMA-adenosine"),
  list(sugar = "ENA",   base = "T", expect = "C12H16N2O6",
       what = "ENA-thymidine")
)

for (a in nucleoside_anchors) {
  f <- add_formulas(STANDARD_DICT[[a$base]]$formula,
                    STANDARD_DICT[[a$sugar]]$formula)
  f <- f - .as_formula(.H2O)
  got <- format_formula(f)
  ok <- identical(got, a$expect)
  if (!ok) fail <- fail + 1L
  cat(sprintf("  %-7s %-5s %-22s %-22s %s   %s\n",
              a$sugar, a$base, got, a$expect, if (ok) "yes" else "NO",
              a$what))
}

cat("\n=== 2. Sugar mass differences from their parent ===\n")
.sugar_mass <- function(code) formula_mass(STANDARD_DICT[[code]]$formula)
diffs <- list(
  list("NMA",   "MOE", 12.9953, "amide N-H in place of MOE's ether CH2"),
  list("UNA",   "r",    2.0157, "ribose + H2 (C2'-C3' bond absent)"),
  list("GNA",   "r",  -58.0055, "glycerol in place of ribose"),
  list("ENA",   "LNA", 14.0157, "one more CH2 in the bridge"),
  list("LNA",   "r",   12.0000, "ribose + CH2 - H2 = ribose + C"),
  list("allyl", "r",   40.0313, "2'-O-allyl"),
  list("AP",    "r",   57.0578, "2'-O-aminopropyl"),
  list("NP",    "d",   -0.9840, "3'-OH -> 3'-NH2"),
  list("FANA",  "f",    0.0000, "arabino vs ribo -- isobaric")
)
for (d in diffs) {
  got <- .sugar_mass(d[[1]]) - .sugar_mass(d[[2]])
  ok <- abs(got - d[[3]]) < 5e-4
  if (!ok) fail <- fail + 1L
  cat(sprintf("  %-6s vs %-4s %+10.4f Da (expected %+9.4f) %s   %s\n",
              d[[1]], d[[2]], got, d[[3]], if (ok) "ok " else "NO ", d[[4]]))
}

cat("\n=== 3. Linkage residues are PO with the -OH swapped ===\n")
po <- formula_mass(STANDARD_DICT[["o"]]$formula)
oh <- formula_mass(parse_formula("OH"))
# code, substituent replacing the -OH, description
swaps <- list(
  list("s",     "SH",        "phosphorothioate"),
  list("mp",    "CH3",       "methylphosphonate"),
  list("pace",  "C2H3O2",    "phosphonoacetate (CH2COOH)"),
  list("msp",   "CH4NO2S",   "mesyl phosphoramidate (NHSO2CH3)"),
  list("pgo",   "C5H10N3",   "phosphoryl guanidine (imidazolidin-2-imine)"),
  list("tmg",   "C5H12N3",   "phosphoryl guanidine (tetramethylguanidine)"),
  list("prp",   "C3H7",      "propyl phosphonate"),
  list("ibu",   "C4H9",      "isobutyl phosphonate"),
  list("chx",   "C6H11",     "cyclohexyl phosphonate"),
  list("mop",   "C4H9O",     "methoxypropyl phosphonate")
)
for (s in swaps) {
  got <- formula_mass(STANDARD_DICT[[s[[1]]]]$formula)
  expect <- po - oh + formula_mass(parse_formula(s[[2]]))
  ok <- abs(got - expect) < 5e-4
  if (!ok) fail <- fail + 1L
  cat(sprintf("  %-6s %10.4f Da  = PO - OH + %-9s (%10.4f) %s   %s\n",
              s[[1]], got, s[[2]], expect, if (ok) "ok " else "NO ", s[[3]]))
}
# thioPACE is PACE with the P=O replaced by P=S.
tp <- formula_mass(STANDARD_DICT[["tpace"]]$formula) -
      formula_mass(STANDARD_DICT[["pace"]]$formula)
ok <- abs(tp - 15.9772) < 5e-4
if (!ok) fail <- fail + 1L
cat(sprintf("  %-6s %+10.4f Da vs pace (expected +15.9772) %s   O -> S\n",
            "tpace", tp, if (ok) "ok " else "NO "))

cat("\n=== 4. Per-linkage mass differences from PO ===\n")
for (code in c("s", "mp", "prp", "ibu", "chx", "mop", "pace", "tpace",
               "msp", "pgo", "tmg")) {
  cat(sprintf("  %-6s %+9.4f Da vs PO   %s\n", code,
              formula_mass(STANDARD_DICT[[code]]$formula) - po,
              STANDARD_DICT[[code]]$name))
}

cat("\n=== 5. New codes are usable in a real sequence ===\n")
# A gapmer whose wings are NMA and whose backbone is mesyl phosphoramidate
# in the wings and PS in the gap -- exercises parsing, validation and the
# formula engine together.
sp <- parse_three_line(
  bases    = "TSASTTTSATAATGSTGG",
  sugars   = c(rep("NMA", 5), rep("d", 8), rep("NMA", 5)),
  linkages = c(rep("msp", 5), rep("s", 7), rep("msp", 5)))
info <- metabolite_mass_info(sp)
cat("  mixed NMA/DNA + msPA/PS 18-mer\n")
cat("   ", info$formula_str, sprintf(" mono %.4f Da\n", info$mono_mass))
mets <- generate_metabolites(sp, list(oligo_name = "test", max_3p = 2,
                                      max_5p = 2, endo = FALSE))
cat("  metabolites generated:", length(mets), "\n")
if (length(mets) != 5) { fail <- fail + 1L; cat("  NO: expected 5\n") }

# Every new code must round-trip through triplet notation.
cat("\n=== 6. New codes round-trip through triplet notation ===\n")
sugar_codes <- c("NMA", "UNA", "GNA", "ENA", "allyl", "AP", "FANA", "NP")
link_codes  <- c("msp", "pace", "tpace", "pgo", "tmg", "mp", "prp", "ibu",
                 "chx", "mop")
for (sg in sugar_codes) {
  for (lk in link_codes) {
    s1 <- parse_three_line("AGT", sg, lk)
    s2 <- parse_input(format_triplet(s1))
    if (!identical(s1$sugars, s2$sugars) || !identical(s1$linkages, s2$linkages)) {
      fail <- fail + 1L
      cat("  NO: round trip failed for", sg, "/", lk, "\n")
    }
  }
}
cat("  all", length(sugar_codes) * length(link_codes),
    "sugar/linkage combinations round-tripped\n")

cat("\n=== 6b. Fatty acid conjugates add the acyl group ===\n")
base_spec <- parse_three_line("AGST", "d", "s")
plain <- metabolite_mass_info(base_spec)$mono_mass
acyl_expect <- list(myristoyl = 210.1984, palmitoyl = 238.2297,
                    stearoyl = 266.2610, docosanoyl = 322.3236)
for (nm in names(acyl_expect)) {
  sp_c <- parse_three_line("AGST", "d", "s", conj5 = nm)
  got <- metabolite_mass_info(sp_c)$mono_mass - plain
  ok <- abs(got - acyl_expect[[nm]]) < 5e-4
  if (!ok) fail <- fail + 1L
  cat(sprintf("  %-11s %+10.4f Da (expected %+9.4f) %s\n",
              nm, got, acyl_expect[[nm]], if (ok) "ok" else "NO"))
}

cat("\n=== 7. The published reference formulas still reproduce ===\n")
ref <- validate_reference(verbose = FALSE)
cat("  validate_reference():", if (isTRUE(ref$ok)) "ok" else "FAILED", "\n")
if (!isTRUE(ref$ok)) fail <- fail + 1L

if (fail > 0) stop(fail, " modification check(s) failed")
cat("\n==== All modification tests passed ====\n")
