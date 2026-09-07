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
# Remove k nucleotides from the 3' end.
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

# Remove k nucleotides from the 5' end.
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
# Returns list(frag5, frag3, frag5_p, frag3_p) -- frag5/frag3 keep the old
# free-terminus behaviour, frag5_p/frag3_p carry the terminal (thio)phosphate.
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
# Return positions where the sugar is deoxyribose (the DNA gap).
# For non-gapmers (all-modified) this returns integer(0).
find_gap <- function(spec) {
  which(spec$sugars == "d")
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

## ---- Main generator --------------------------------------------------------
# opts:
#   oligo_name   : prefix for metabolite names (e.g. "inotersen")
#   max_3p       : max 3' exonuclease truncations (default 10)
#   max_5p       : max 5' exonuclease truncations (default 10)
#   endo         : include endonuclease fragments? (default TRUE)
#   endo_sites   : "all" | "gap" | integer vector of cleavage positions
#   min_frag_len : minimum fragment length to keep (default 3)
#   dedupe       : collapse structurally identical species? (default TRUE;
#                  see dedupe_metabolites() below)
generate_metabolites <- function(spec, opts = list(), dict = STANDARD_DICT) {
  oligo_name <- opts$oligo_name %||% "OLIGO"
  max_3p     <- opts$max_3p %||% 10
  max_5p     <- opts$max_5p %||% 10
  endo       <- opts$endo %||% TRUE
  endo_sites <- opts$endo_sites %||% "all"
  min_frag   <- opts$min_frag_len %||% 3
  dedupe     <- opts$dedupe %||% TRUE

  max_3p <- min(max_3p, spec$n - 1)
  max_5p <- min(max_5p, spec$n - 1)
  mets <- list()
  ctr <- 0
  add <- function(m) { ctr <<- ctr + 1; m$id <- sprintf("M%02d", ctr); mets[[ctr]] <<- m }

  # Parent
  add(.make_met(NA, oligo_name, "parent", "parent", NA, spec, dict = dict))

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
  mets
}

## ---- De-duplicate structurally identical species ---------------------------
# With endo = TRUE, endo_sites = "all", a 3' exonuclease truncation and an
# endonuclease 5' fragment (or a 5' truncation and an endonuclease 3'
# fragment) can be the exact same molecule -- same sequence, same termini,
# same mass -- differing only in the `kind` label attached by whichever
# route produced it. Left uncollapsed, degradation_summary() (which sums
# signal by `kind`) counts that one species twice on the degradant side,
# inflating % degradation. Collapse on a structural key (sequence + sugars
# + linkages + terminal conjugates) and keep every route that produces the
# surviving species in `kind_all`/`n_routes`, so the reporting stays
# informative even after the count no longer double-counts mass.
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
# Returns a data.frame with one row per metabolite (spec vectors as strings).
metabolite_table <- function(mets) {
  do.call(rbind, lapply(mets, function(m) data.frame(
    id = m$id, name = m$name, kind = m$kind, n = m$n,
    n_ps = m$n_ps, n_po = m$n_po,
    bases = paste(m$bases, collapse = ""),
    modification = m$modification, site = ifelse(is.null(m$site) || is.na(m$site), "", m$site),
    stringsAsFactors = FALSE
  )))
}
