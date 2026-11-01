# test_structure_annotation.R -- sequence-level cleavage-map model and SVG
# (R/structure_annotation.R): fragment span, termini, modification badges,
# sugar coloring, and well-formed SVG output.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "structure_annotation.R")) {
  source(file.path(.pkg_root, "R", .f))
}
check <- function(label, ok) {
  cat(sprintf("%-82s %s\n", label, if (isTRUE(ok)) "PASS" else "FAIL"))
  if (!isTRUE(ok)) stop("FAILED: ", label)
}
wellformed <- function(svg) {
  if (!requireNamespace("xml2", quietly = TRUE)) return(TRUE)  # skip if absent
  tryCatch({ xml2::read_xml(svg); TRUE }, error = function(e) FALSE)
}
cat("==== Structure annotation suite ====\n\n")

# 5-10-5 MOE/DNA/MOE gapmer, fully phosphorothioate
N <- 16
spec <- list(bases = strsplit("TCGTCATGACTTGCAG", "")[[1]],
             sugars = c(rep("e", 5), rep("d", 6), rep("e", 5)),
             linkages = c(rep("s", 15), NA), conj5 = "none", conj3 = "none", n = N)
mets <- generate_metabolites(spec, list(oligo_name = "demo", endo = TRUE,
                                        endo_sites = "gap", max_3p = 3, max_5p = 2))
pick <- function(kind) Filter(function(m) m$kind == kind, mets)[[1]]

## ---- fragment span --------------------------------------------------------
m_par <- metabolite_annotation_model(mets[[1]], parent = spec)
check("parent spans the whole sequence", identical(m_par$span, c(1L, N)))

e3 <- pick("exo_3p")                          # removes k from the 3' end
m_e3 <- metabolite_annotation_model(e3, parent = spec)
check("3' exonuclease span starts at 1", m_e3$span[1] == 1L && m_e3$span[2] == N - e3$site)
check("3' exonuclease fades the lost 3' residues",
      all(!m_e3$residues$in_frag[(m_e3$span[2] + 1):N]))

e5 <- pick("exo_5p")                          # removes k from the 5' end
m_e5 <- metabolite_annotation_model(e5, parent = spec)
check("5' exonuclease span ends at N", m_e5$span[2] == N && m_e5$span[1] == e5$site + 1L)

## ---- termini --------------------------------------------------------------
m_frag_p <- metabolite_annotation_model(pick("endo_5frag_p"), parent = spec)
check("endonuclease 3'-phosphate fragment labels a 3'-p terminus",
      m_frag_p$term3 == "3'-p" || m_frag_p$term3 == "3'-ps")
check("cyclic phosphate relabels the 3' terminus when applicable",
      {
        # force a ribose 3'-phosphate fragment
        rna <- list(bases = c("A","C","G","U"), sugars = rep("r",4),
                    linkages = c("s","s","s",NA), conj5="none", conj3="none", n=4)
        rm <- Filter(can_form_cyclic_phosphate,
                     generate_metabolites(rna, list(oligo_name="r", endo=TRUE,
                       endo_sites=2, max_3p=1, max_5p=1, min_frag_len=1)))[[1]]
        mm <- metabolite_annotation_model(rm, mod = "cyclic_phosphate")
        grepl("cyclic", mm$term3)
      })

## ---- modification badges --------------------------------------------------
m_deam <- metabolite_annotation_model(mets[[1]], parent = spec,
                                      mod = "deamination", mod_site = 2)
check("localized deamination badges exactly one position (not a candidate)",
      nrow(m_deam$mods) == 1 && !m_deam$mods$candidate[1] && m_deam$mods$position == 2)

m_dep <- metabolite_annotation_model(mets[[1]], parent = spec,
                                     mod = "depurination", mod_site = NA)
check("unlocalized depurination badges all purines as candidates",
      nrow(m_dep$mods) == sum(mets[[1]]$bases %in% c("A", "G")) &&
        all(m_dep$mods$candidate))

## ---- sugar coloring -------------------------------------------------------
check("DNA gap and MOE wings get different colors",
      length(unique(m_par$residues$color[m_par$residues$sugar %in% c("d", "e")])) == 2)

## ---- SVG output -----------------------------------------------------------
for (m in list(mets[[1]], e3, pick("endo_5frag_p"))) {
  svg <- annotate_metabolite_svg(m, parent = spec)
  check(sprintf("SVG renders well-formed for %s", m$kind),
        grepl("^<svg", svg) && wellformed(svg))
}
svg_mod <- annotate_metabolite_svg(mets[[1]], parent = spec,
                                   mod = "deamination", mod_site = 2, k_oxid = 3)
check("SVG with modification + desulfurization note is well-formed",
      wellformed(svg_mod))

cat("\nAll structure-annotation checks passed.\n")
