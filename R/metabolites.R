# =============================================================================
# metabolites.R
# Metabolite library generation for oligonucleotide metabolite identification.
#
# Generates the theoretical metabolite series from a canonical oligo_spec:
#   - Parent (full-length)
#   - 3' exonuclease truncations  (n-1 .. n-k, removing from 3' end)
#   - 5' exonuclease truncations  (n-1 .. n-k, removing from 5' end)
#   - Endonuclease internal fragments (5' and 3' fragments at each cleavage)
#
# Biotransformations that alter mass but not sequence (PS->PO oxidation series,
# +Na/+K/+NH4 adducts, depurination, +O) are recorded as per-metabolite metadata
# here and computed as mass-shift variants in mass_isotope.R.
#
# Grounded in:
#   - Eluforsen (Kim et al.): 3'/5' shortmer series n-1..n-k; PS->PO oxidation
#     columns 0-6; pyrimidines cleaved faster than purines.
#   - FMVS (Ye et al.): endonuclease (S1/RNase) internal cleavage then exo trim.
#   - OligoDistiller: gapmers cleaved in DNA gap first, then exonuclease.
# =============================================================================

## ---- Truncation helpers ----------------------------------------------------
#' Truncate an oligo_spec from the 3' end (exonuclease truncation)
#'
#' @param spec An `oligo_spec` list (see [parse_input()]).
#' @param k Number of nucleotides to remove, `1` to `n-1`.
#' @return A truncated spec-like list (`bases`, `sugars`, `linkages`,
#'   `conj5`, `conj3 = "none"` (lost with the 3' end), `n`).
#' @seealso [truncate_5p()], [generate_metabolites()]
#' @export
truncate_3p <- function(spec, k) {
  if (k < 1 || k >= spec$n) stop("3' truncation k must be in 1..(n-1)")
  n2 <- spec$n - k
  list(
    bases    = spec$bases[1:n2],
    sugars   = spec$sugars[1:n2],
    linkages = {lk <- spec$linkages[1:n2]; lk[n2] <- NA_character_; lk},
    conj5    = spec$conj5,
    conj3    = "none",                 # 3' conjugate lost with the 3' end
    n        = n2
  )
}

#' Truncate an oligo_spec from the 5' end (exonuclease truncation)
#'
#' @param spec An `oligo_spec` list (see [parse_input()]).
#' @param k Number of nucleotides to remove, `1` to `n-1`.
#' @return A truncated spec-like list (`bases`, `sugars`, `linkages`,
#'   `conj5 = "none"` (lost with the 5' end), `conj3`, `n`).
#' @seealso [truncate_3p()], [generate_metabolites()]
#' @export
truncate_5p <- function(spec, k) {
  if (k < 1 || k >= spec$n) stop("5' truncation k must be in 1..(n-1)")
  n2 <- spec$n - k
  idx <- (k + 1):spec$n
  list(
    bases    = spec$bases[idx],
    sugars   = spec$sugars[idx],
    linkages = {lk <- spec$linkages[idx]; lk[n2] <- NA_character_; lk},
    conj5    = "none",                 # 5' conjugate lost with the 5' end
    conj3    = spec$conj3,
    n        = n2
  )
}

# Endonuclease cleavage between position i and i+1 (1-indexed, i in 1..n-1).
# A single phosphodiester bond hydrolysis leaves the bridging phosphate on
# exactly one product, never neither -- giving both products a free
# hydroxyl (as this function did previously) drops the entire bridging
# phosphate/phosphorothioate group from the pair's combined mass. Which
# side keeps the phosphate depends on the nuclease, so both mass-balanced
# variants are returned: the free-3'-OH/free-5'-OH pair (as before, still
# the majority mechanism for most DNases) *and* the complementary pair
# where the phosphate stays on the other product (RNase A/RNase H-like).
# The cap itself must match the cleaved bond's own chemistry -- a PS
# linkage leaves a terminal THIOphosphate, not a plain phosphate, or the
# pair is still short by exactly one S->O swap (15.9772 Da); sulfur
# content (not the literal "s"/"u" code) decides which cap to use, same
# rationale as count_linkages() above.
#' Split an oligo_spec at an endonuclease cleavage site
#'
#' Cleaves a single phosphodiester (or phosphorothioate) bond between
#' position `i` and `i+1`. Hydrolysis leaves the bridging phosphate on
#' exactly one product, never neither and never both -- which side keeps
#' it depends on the nuclease, so both mass-balanced outcomes are
#' returned: the free-3'-OH/free-5'-OH pair (the majority mechanism for
#' most DNases) and the complementary pair where the phosphate stays on
#' the other product (RNase A/RNase H-like). The terminal cap matches the
#' cleaved bond's own chemistry -- a phosphorothioate linkage leaves a
#' terminal thiophosphate, decided by the linkage's sulfur content, not
#' by its literal code.
#'
#' @param spec An `oligo_spec` list (see [parse_input()]).
#' @param i Cleavage site: the bond between position `i` and `i+1`,
#'   `1` to `n-1`.
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return A list of four spec-like fragments: `frag5`/`frag3` (free
#'   3'-OH / free 5'-OH -- the bridging phosphate went to the OTHER
#'   fragment), `frag5_p`/`frag3_p` (carrying the terminal (thio)phosphate
#'   instead). Use `frag5` with `frag3_p`, or `frag5_p` with `frag3`, as
#'   a mass-balanced pair -- never `frag5` with `frag3` together (that
#'   drops the bridging phosphate entirely).
#' @seealso [generate_metabolites()]
#' @export
endo_cleave <- function(spec, i, dict = STANDARD_DICT) {
  if (i < 1 || i >= spec$n) stop("endo cleavage site i must be in 1..(n-1)")
  lk_formula <- dict[[spec$linkages[i]]]$formula
  is_ps <- !is.null(lk_formula) && isTRUE(unname(lk_formula[["S"]]) > 0)
  cap5 <- if (is_ps) "3'-thiophosphate" else "3'-phosphate"
  cap3 <- if (is_ps) "5'-thiophosphate" else "5'-phosphate"
  # 5' fragment: positions 1..i, new 3' end at position i
  mk5 <- function(conj3) list(
    bases    = spec$bases[1:i],
    sugars   = spec$sugars[1:i],
    linkages = {lk <- spec$linkages[1:i]; lk[i] <- NA_character_; lk},
    conj5    = spec$conj5,
    conj3    = conj3,
    n        = i)
  # 3' fragment: positions (i+1)..n, new 5' end at position i+1
  idx <- (i + 1):spec$n
  mk3 <- function(conj5) list(
    bases    = spec$bases[idx],
    sugars   = spec$sugars[idx],
    linkages = {lk <- spec$linkages[idx]; lk[length(idx)] <- NA_character_; lk},
    conj5    = conj5,
    conj3    = spec$conj3,
    n        = length(idx))
  list(
    frag5   = mk5("none"),   # free 3'-OH (phosphate stayed on frag3_p)
    frag3   = mk3("none"),   # free 5'-OH (phosphate stayed on frag5_p)
    frag5_p = mk5(cap5),     # 3'-(thio)phosphate (phosphate stayed on frag5)
    frag3_p = mk3(cap3))     # 5'-(thio)phosphate (phosphate stayed on frag3)
}

## ---- Linkage accounting ----------------------------------------------------
# Count PS vs PO bonds in a spec (for the oxidation-series dimension).
# PS is identified by sulfur content in the linkage's own formula, not by
# matching the literal codes "s"/"u" -- fragments.R's .is_ps_linkage()
# already does this correctly; matching codes alone silently gave n_ps = 0
# (no desulfurization series at all) for any other sulfur-bearing backbone
# in the dictionary (mesyl-phosphoramidate, thio-PACE, a user-defined PS
# analogue added via the Custom Chemistry table).
#' Count phosphorothioate vs. phosphodiester linkages in an oligo_spec
#'
#' A linkage counts as phosphorothioate by sulfur content in its own
#' dictionary formula, not by matching the literal codes `"s"`/`"u"` --
#' so a custom sulfur-bearing backbone (mesyl phosphoramidate, thioPACE,
#' a user-defined PS analogue) is counted correctly too.
#'
#' @param spec An `oligo_spec` list (see [parse_input()]).
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return A list: `n_ps`, `n_po`, `n_bonds` (total internal linkages).
#' @seealso [generate_metabolites()] (drives the PS -> PO oxidation series)
#' @export
count_linkages <- function(spec, dict = STANDARD_DICT) {
  lk <- spec$linkages[!is.na(spec$linkages)]
  n_ps <- sum(vapply(lk, function(code) {
    f <- dict[[code]]$formula
    !is.null(f) && isTRUE(unname(f[["S"]]) > 0)
  }, logical(1)))
  n_po <- length(lk) - n_ps
  list(n_ps = n_ps, n_po = n_po, n_bonds = length(lk))
}

## ---- Gap detection (for gapmer endonuclease cleavage) ----------------------
#' Locate a gapmer's DNA gap
#'
#' @param spec An `oligo_spec` list (see [parse_input()]).
#' @return Integer positions where the sugar is deoxyribose (`"d"`) --
#'   the DNA gap. `integer(0)` for a non-gapmer (no deoxyribose
#'   positions).
#' @seealso [generate_metabolites()] (`endo_sites = "gap"`)
#' @export
find_gap <- function(spec) {
  which(spec$sugars == "d")
}

#' Does this oligo look like a gapmer?
#'
#' A gapmer carries an internal 2'-deoxyribose (DNA) gap flanked by
#' 2'-modified wings. This is `TRUE` when the spec has at least one DNA
#' sugar and its 5'- and 3'-terminal sugars are both non-DNA, i.e. the
#' DNA sits in the middle, not at the ends (a uniformly-DNA oligo is a
#' first-generation PS DNA, not a gapmer, and returns `FALSE`).
#'
#' @param spec An `oligo_spec` list (see [parse_input()]).
#' @return A single logical.
#' @seealso [find_gap()], [generate_metabolites()]
#' @export
is_gapmer <- function(spec) {
  s <- spec$sugars
  n <- length(s)
  if (is.null(n) || n < 3) return(FALSE)
  any(s == "d") && !identical(s[1], "d") && !identical(s[n], "d")
}

# Pyrimidines (including the common methylated/modified forms). Pyrimidine
# termini are cleaved faster than purines by the exo/endonucleases.
.PYRIMIDINES <- c("C", "T", "U", "S", "m5C", "hm5C", "m5U", "psi")
.is_pyrimidine <- function(base) !is.null(base) && base %in% .PYRIMIDINES

# Heuristic relative-abundance prior (a prior, NOT validated kinetics):
# 3'-exonuclease is the dominant plasma pathway, 5'-exonuclease is minor,
# endonuclease is minor except within a gapmer's DNA gap. Each successive
# exonuclease removal is less abundant (geometric decay), and a pyrimidine
# terminus is removed faster than a purine.
.NUCLEASE_PRIOR <- list(exo3_base = 1.0, exo5_base = 0.35, endo_base = 0.12,
                        endo_gap_base = 0.5, decay = 0.8, pyr_mult = 1.3,
                        conj_base = 0.8)

#' Attach a heuristic nuclease-preference weight to each metabolite
#'
#' Scores every metabolite with a relative-abundance prior (parent = 1):
#' 3'-exonuclease > 5'-exonuclease > endonuclease, with geometric decay by
#' truncation depth and a boost for pyrimidine termini and for cleavage
#' inside a gapmer's DNA gap. This is an ordering prior for prioritising the
#' inclusion list, not a kinetic prediction, and carries no claim of
#' quantitative accuracy.
#'
#' @param mets A list of metabolite objects (see [generate_metabolites()]).
#' @param spec The parent `oligo_spec` (for sequence/gap lookups).
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return `mets` with a numeric `rel_weight` (0-1) added to each entry.
#' @seealso [generate_metabolites()], [prm_inclusion_list()]
#' @export
assign_nuclease_weights <- function(mets, spec, dict = STANDARD_DICT) {
  n <- spec$n
  gap <- if (is_gapmer(spec)) find_gap(spec) else integer(0)
  P <- .NUCLEASE_PRIOR
  pyr <- function(b) if (.is_pyrimidine(b)) P$pyr_mult else 1
  for (i in seq_along(mets)) {
    m <- mets[[i]]
    kinds <- m$kind_all %||% m$kind
    k <- suppressWarnings(as.integer(m$site))
    w <- if ("parent" %in% kinds) {
      1.0
    } else if ("conj_catab" %in% kinds) {
      P$conj_base * P$decay^(if (is.na(k)) 1 else k)
    } else if ("exo_3p" %in% kinds && !is.na(k)) {
      P$exo3_base * P$decay^(k - 1) * pyr(spec$bases[n - k + 1])
    } else if ("exo_5p" %in% kinds && !is.na(k)) {
      P$exo5_base * P$decay^(k - 1) * pyr(spec$bases[k])
    } else if (any(grepl("^endo", kinds)) && !is.na(k)) {
      base <- if (length(gap) >= 2 && k >= min(gap) && k < max(gap)) P$endo_gap_base else P$endo_base
      base * pyr(spec$bases[k])
    } else 0.1
    mets[[i]]$rel_weight <- round(min(1, w), 4)
  }
  mets
}

# Relevance tier for the inclusion list. Tier 1 is the set that dominates
# real plasma/tissue/urine profiles (parent, the full 3'-exonuclease ladder,
# and short 5'-exonuclease products); Tier 2 is the wider net (deeper 5'-exo
# and endonuclease internal fragments), usually minor. Computed from the
# collapsed kinds so a species reachable by a Tier-1 route is Tier 1.
.met_tier <- function(m, max_exo5_tier1 = 3L) {
  kinds <- m$kind_all %||% m$kind
  if ("parent" %in% kinds) return(1L)
  if ("conj_catab" %in% kinds) return(1L)
  if ("exo_3p" %in% kinds) return(1L)
  if ("exo_5p" %in% kinds && !is.null(m$site) && !is.na(m$site) &&
      as.integer(m$site) <= max_exo5_tier1) return(1L)
  2L
}

## ---- Metabolite object builder ---------------------------------------------
.make_met <- function(id, name, kind, modification, site, sp, parent_id = NA,
                       dict = STANDARD_DICT) {
  lk <- count_linkages(sp, dict)
  list(
    id = id, name = name, kind = kind, modification = modification,
    site = site, n = sp$n,
    bases = sp$bases, sugars = sp$sugars, linkages = sp$linkages,
    conj5 = sp$conj5 %||% "none", conj3 = sp$conj3 %||% "none",
    n_ps = lk$n_ps, n_po = lk$n_po, n_bonds = lk$n_bonds,
    parent_id = parent_id
  )
}

## ---- Conjugate (GalNAc) catabolism -----------------------------------------
# GalNAc catabolic ladders: ordered conjugate codes from the intact cluster
# down to the fully de-conjugated terminus ("none"). Each step past the
# intact cluster is one more hydrolysed GalNAc sugar (-C8H13NO5,
# -203.0794 Da); the scaffold is all three sugars gone with the linker
# still attached; "none" is the naked oligo released by linker cleavage.
.GALNAC_LADDER <- list(
  GalNAc               = c("none"),
  GalNAc3              = c("GalNAc3_m1", "GalNAc3_m2", "GalNAc3_m3", "none"),
  GalNAc3_triantennary = c("GalNAc3t_m1", "GalNAc3t_m2", "GalNAc3t_m3", "none")
)
.galnac_step_label <- function(step, n_steps) {
  if (step == n_steps) return("linker cleavage (de-conjugated oligo)")
  if (step == n_steps - 1L) return(sprintf("-%d GalNAc (scaffold)", step))
  sprintf("-%d GalNAc", step)
}

#' Is this oligo carrying a (recognised) GalNAc conjugate?
#'
#' `TRUE` when either terminus carries a GalNAc-type conjugate for which a
#' catabolic ladder is defined (`GalNAc`, `GalNAc3`, or
#' `GalNAc3_triantennary`).
#'
#' @param spec An `oligo_spec` list (see [parse_input()]).
#' @return A single logical.
#' @seealso [generate_conjugate_catabolites()]
#' @export
is_galnac_conjugated <- function(spec) {
  conj <- c(spec$conj5 %||% "none", spec$conj3 %||% "none")
  any(conj %in% names(.GALNAC_LADDER))
}

#' Generate the conjugate-catabolism metabolite series
#'
#' For a GalNAc-conjugated strand, builds the characteristic catabolite
#' ladder that dominates early DMPK profiles: stepwise loss of each GalNAc
#' sugar (-203.0794 Da per residue), the de-glycosylated linker scaffold,
#' and the fully de-conjugated (naked) oligo released by linker cleavage.
#' Each catabolite is a full-length species with a mass-correct partial
#' conjugate (see the `GalNAc3_m*`/`GalNAc3t_m*` codes in
#' `CONJUGATE_FORMULAS`), so it flows through mass, structure-annotation,
#' and fragment tooling unchanged.
#'
#' @param spec An `oligo_spec` list (see [parse_input()]).
#' @param oligo_name Prefix for metabolite names.
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return A list of metabolite objects (`kind = "conj_catab"`), one per
#'   ladder step on each conjugated terminus; empty if no GalNAc conjugate
#'   is present.
#' @seealso [is_galnac_conjugated()], [generate_metabolites()]
#' @export
generate_conjugate_catabolites <- function(spec, oligo_name = "OLIGO",
                                           dict = STANDARD_DICT) {
  out <- list()
  for (end in c("conj5", "conj3")) {
    code <- spec[[end]] %||% "none"
    ladder <- .GALNAC_LADDER[[code]]
    if (is.null(ladder)) next
    end_lab <- if (end == "conj5") "5'" else "3'"
    n_steps <- length(ladder)
    for (step in seq_len(n_steps)) {
      sp <- spec
      sp[[end]] <- ladder[step]
      lab <- .galnac_step_label(step, n_steps)
      out[[length(out) + 1L]] <- .make_met(
        NA, sprintf("%s %s %s", oligo_name, end_lab, lab),
        "conj_catab", sprintf("%s GalNAc catabolism: %s", end_lab, lab),
        step, sp, dict = dict)
    }
  }
  out
}

## ---- Main generator --------------------------------------------------------
#' Generate the theoretical metabolite library for an oligonucleotide
#'
#' The main library-generation entry point: builds the parent, the 3'
#' and 5' exonuclease truncation series, and (optionally) endonuclease
#' internal fragments (both mass-balanced termini per cleavage site --
#' see [endo_cleave()]), then collapses structurally identical species
#' produced by more than one route (see [dedupe_metabolites()]).
#'
#' @param spec An `oligo_spec` list (see [parse_input()]).
#' @param opts A list of options:
#'   \describe{
#'     \item{oligo_name}{Prefix for metabolite names, e.g. `"inotersen"`.}
#'     \item{max_3p, max_5p}{Max exonuclease truncations from each end
#'       (default 10 each; capped at `spec$n - 1`).}
#'     \item{endo}{Include endonuclease internal fragments? (default `TRUE`)}
#'     \item{endo_sites}{`"all"` (every internal bond), `"gap"` (only
#'       within a gapmer's DNA gap -- see [find_gap()]), or an integer
#'       vector of specific cleavage positions.}
#'     \item{min_frag_len}{Minimum fragment length to keep (default 3).}
#'     \item{dedupe}{Collapse structurally identical species across
#'       routes? (default `TRUE`; see [dedupe_metabolites()])}
#'     \item{conj_catab}{Include the GalNAc conjugate-catabolism ladder
#'       (stepwise sugar loss + linker cleavage) when a GalNAc conjugate
#'       is present? (default `TRUE`; see
#'       [generate_conjugate_catabolites()])}
#'   }
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return A list of metabolite objects, each with `id` (`"M01"`,
#'   `"M02"`, ...), `name`, `kind` (`"parent"`/`"conj_catab"`/`"exo_3p"`/
#'   `"exo_5p"`/`"endo_5frag"`/`"endo_3frag"`/`"endo_5frag_p"`/
#'   `"endo_3frag_p"`), `modification`, `site`, `n`,
#'   `bases`/`sugars`/`linkages`, `conj5`/`conj3`, `n_ps`/`n_po`/`n_bonds`,
#'   `tier`, `rel_weight` (nuclease-preference prior), `parent_id`, and
#'   (after deduplication) `kind_all`/`n_routes`.
#' @seealso [metabolite_mass_info()], [metabolite_table()],
#'   [match_ms1()]
#' @export
generate_metabolites <- function(spec, opts = list(), dict = STANDARD_DICT) {
  # Use [[ for exact matching: $ does partial matching, so opts$endo would
  # silently resolve to opts$endo_sites when only the latter is supplied.
  oligo_name <- opts[["oligo_name"]] %||% "OLIGO"
  max_3p     <- opts[["max_3p"]] %||% 10
  max_5p     <- opts[["max_5p"]] %||% 10
  endo       <- opts[["endo"]] %||% TRUE
  endo_sites <- opts[["endo_sites"]] %||% "auto"
  min_frag   <- opts[["min_frag_len"]] %||% 5
  dedupe     <- opts[["dedupe"]] %||% TRUE
  conj_catab <- opts[["conj_catab"]] %||% TRUE

  # "auto": restrict endonuclease cleavage to the DNA gap for gapmers (the
  # dominant endo route for that architecture), and search all internal bonds
  # otherwise. Keeps the default inclusion list from exploding with short,
  # biologically negligible internal fragments on a gapmer.
  if (identical(endo_sites, "auto")) endo_sites <- if (is_gapmer(spec)) "gap" else "all"

  max_3p <- min(max_3p, spec$n - 1)
  max_5p <- min(max_5p, spec$n - 1)
  mets <- list()
  ctr <- 0
  add <- function(m) { ctr <<- ctr + 1; m$id <- sprintf("M%02d", ctr); mets[[ctr]] <<- m }

  # Parent
  add(.make_met(NA, oligo_name, "parent", "parent", NA, spec, dict = dict))

  # Conjugate catabolism (stepwise GalNAc loss + linker cleavage). For a
  # GalNAc-conjugated strand these dominate the early profile, so they are
  # on by default whenever a GalNAc conjugate is present.
  if (conj_catab && is_galnac_conjugated(spec)) {
    for (m in generate_conjugate_catabolites(spec, oligo_name, dict)) add(m)
  }

  # 3' exonuclease series
  for (k in seq_len(max_3p)) {
    sp <- truncate_3p(spec, k)
    add(.make_met(NA, sprintf("%s 3' N-%d", oligo_name, k),
                  "exo_3p", sprintf("3' exonuclease (-%d nt)", k), k, sp, dict = dict))
  }

  # 5' exonuclease series
  for (k in seq_len(max_5p)) {
    sp <- truncate_5p(spec, k)
    add(.make_met(NA, sprintf("%s 5' N-%d", oligo_name, k),
                  "exo_5p", sprintf("5' exonuclease (-%d nt)", k), k, sp, dict = dict))
  }

  # Endonuclease internal fragments
  if (endo) {
    sites <- if (is.character(endo_sites)) {
      if (endo_sites == "all") seq_len(spec$n - 1)
      else if (endo_sites == "gap") {
        g <- find_gap(spec)
        if (length(g) < 2) integer(0) else seq(min(g), max(g) - 1)
      } else integer(0)
    } else as.integer(endo_sites)
    sites <- sites[sites >= 1 & sites < spec$n]
    for (i in sites) {
      cl <- endo_cleave(spec, i, dict)
      if (cl$frag5$n >= min_frag) {
        add(.make_met(NA, sprintf("%s Endo 5'frag @%d", oligo_name, i),
                      "endo_5frag", sprintf("endonuclease cleavage after pos %d (5' fragment, free 3'-OH)", i),
                      i, cl$frag5, dict = dict))
        add(.make_met(NA, sprintf("%s Endo 5'frag @%d (3'-phosphate)", oligo_name, i),
                      "endo_5frag_p", sprintf("endonuclease cleavage after pos %d (5' fragment, 3'-phosphate)", i),
                      i, cl$frag5_p, dict = dict))
      }
      if (cl$frag3$n >= min_frag) {
        add(.make_met(NA, sprintf("%s Endo 3'frag @%d", oligo_name, i),
                      "endo_3frag", sprintf("endonuclease cleavage after pos %d (3' fragment, free 5'-OH)", i),
                      i, cl$frag3, dict = dict))
        add(.make_met(NA, sprintf("%s Endo 3'frag @%d (5'-phosphate)", oligo_name, i),
                      "endo_3frag_p", sprintf("endonuclease cleavage after pos %d (3' fragment, 5'-phosphate)", i),
                      i, cl$frag3_p, dict = dict))
      }
    }
  }
  if (dedupe) mets <- dedupe_metabolites(mets)
  # Relevance tier for inclusion-list prioritisation (1 = dominant species,
  # 2 = wider net). See .met_tier().
  for (i in seq_along(mets)) mets[[i]]$tier <- .met_tier(mets[[i]])
  # Heuristic nuclease-preference weight (ordering prior, not kinetics).
  mets <- assign_nuclease_weights(mets, spec, dict)
  mets
}

## ---- De-duplicate structurally identical species ---------------------------
#' Collapse structurally identical metabolites produced by multiple routes
#'
#' With `endo = TRUE, endo_sites = "all"`, a 3' exonuclease truncation
#' and an endonuclease 5' fragment (or a 5' truncation and an
#' endonuclease 3' fragment) can be the exact same molecule -- same
#' sequence, same termini, same mass -- differing only in the `kind`
#' label attached by whichever route produced it. Left uncollapsed,
#' [degradation_summary()] (which sums signal by `kind`) would count
#' that one species twice on the degradant side, inflating % degradation.
#' Collapses on a structural key (sequence + sugars + linkages + terminal
#' conjugates), keeping every route that produces the surviving species
#' in `kind_all`/`n_routes` so the reporting stays informative.
#'
#' @param mets A list of metabolite objects (see [generate_metabolites()]).
#' @return The deduplicated list -- one entry per distinct structure,
#'   each with `kind_all` (every `kind` that produced this structure) and
#'   `n_routes` added.
#' @seealso [generate_metabolites()]
#' @export
dedupe_metabolites <- function(mets) {
  key <- vapply(mets, function(m) paste(paste(m$bases, collapse = ""),
                                        paste(m$sugars, collapse = ""),
                                        paste(m$linkages, collapse = ""),
                                        m$conj5, m$conj3, sep = "|"), "")
  keep <- !duplicated(key)
  out <- mets[keep]
  for (i in seq_along(out)) {
    routes <- unique(vapply(mets[key == key[keep][i]], function(m) m$kind, ""))
    out[[i]]$kind_all <- routes
    out[[i]]$n_routes <- length(routes)
  }
  out
}

## ---- Flatten to a display table -------------------------------------------
#' Flatten a metabolite list to a display data.frame
#'
#' @param mets A list of metabolite objects (see [generate_metabolites()]).
#' @return A data.frame, one row per metabolite: `id`, `name`, `kind`,
#'   `n`, `n_ps`, `n_po`, `bases` (collapsed to a single string),
#'   `modification`, `site`.
#' @seealso [generate_metabolites()]
#' @export
metabolite_table <- function(mets) {
  do.call(rbind, lapply(mets, function(m) data.frame(
    id = m$id, name = m$name, kind = m$kind, n = m$n,
    n_ps = m$n_ps, n_po = m$n_po,
    bases = paste(m$bases, collapse = ""),
    modification = m$modification, site = ifelse(is.null(m$site) || is.na(m$site), "", m$site),
    tier = m$tier %||% NA_integer_,
    rel_weight = m$rel_weight %||% NA_real_,
    stringsAsFactors = FALSE
  )))
}
