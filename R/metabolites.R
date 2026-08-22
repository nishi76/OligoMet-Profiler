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
# Returns list(frag5 = positions 1..i, frag3 = positions i+1..n).
endo_cleave <- function(spec, i) {
  if (i < 1 || i >= spec$n) stop("endo cleavage site i must be in 1..(n-1)")
  # 5' fragment: positions 1..i, new 3' end at position i
  f5 <- list(
    bases    = spec$bases[1:i],
    sugars   = spec$sugars[1:i],
    linkages = {lk <- spec$linkages[1:i]; lk[i] <- NA_character_; lk},
    conj5    = spec$conj5,
    conj3    = "none",                 # cleavage exposes a new 3'-OH
    n        = i)
  # 3' fragment: positions (i+1)..n, new 5' end at position i+1
  idx <- (i + 1):spec$n
  f3 <- list(
    bases    = spec$bases[idx],
    sugars   = spec$sugars[idx],
    linkages = {lk <- spec$linkages[idx]; lk[length(idx)] <- NA_character_; lk},
    conj5    = "none",                 # cleavage exposes a new 5'-OH
    conj3    = spec$conj3,
    n        = length(idx))
  list(frag5 = f5, frag3 = f3)
}

## ---- Linkage accounting ----------------------------------------------------
# Count PS vs PO bonds in a spec (for the oxidation-series dimension).
count_linkages <- function(spec) {
  lk <- spec$linkages[!is.na(spec$linkages)]
  n_ps <- sum(lk %in% c("s", "u"))     # 'u' is the PS stereochem variant
  n_po <- sum(lk == "o")
  list(n_ps = n_ps, n_po = n_po, n_bonds = length(lk))
}

## ---- Gap detection (for gapmer endonuclease cleavage) ----------------------
# Return positions where the sugar is deoxyribose (the DNA gap).
# For non-gapmers (all-modified) this returns integer(0).
find_gap <- function(spec) {
  which(spec$sugars == "d")
}

## ---- Metabolite object builder ---------------------------------------------
.make_met <- function(id, name, kind, modification, site, sp, parent_id = NA) {
  lk <- count_linkages(sp)
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
generate_metabolites <- function(spec, opts = list()) {
  oligo_name <- opts$oligo_name %||% "OLIGO"
  max_3p     <- opts$max_3p %||% 10
  max_5p     <- opts$max_5p %||% 10
  endo       <- opts$endo %||% TRUE
  endo_sites <- opts$endo_sites %||% "all"
  min_frag   <- opts$min_frag_len %||% 3

  max_3p <- min(max_3p, spec$n - 1)
  max_5p <- min(max_5p, spec$n - 1)
  mets <- list()
  ctr <- 0
  add <- function(m) { ctr <<- ctr + 1; m$id <- sprintf("M%02d", ctr); mets[[ctr]] <<- m }

  # Parent
  add(.make_met(NA, oligo_name, "parent", "parent", NA, spec))

  # 3' exonuclease series
  for (k in seq_len(max_3p)) {
    sp <- truncate_3p(spec, k)
    add(.make_met(NA, sprintf("%s 3' N-%d", oligo_name, k),
                  "exo_3p", sprintf("3' exonuclease (-%d nt)", k), k, sp))
  }

  # 5' exonuclease series
  for (k in seq_len(max_5p)) {
    sp <- truncate_5p(spec, k)
    add(.make_met(NA, sprintf("%s 5' N-%d", oligo_name, k),
                  "exo_5p", sprintf("5' exonuclease (-%d nt)", k), k, sp))
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
      cl <- endo_cleave(spec, i)
      if (cl$frag5$n >= min_frag)
        add(.make_met(NA, sprintf("%s Endo 5'frag @%d", oligo_name, i),
                      "endo_5frag", sprintf("endonuclease cleavage after pos %d (5' fragment)", i),
                      i, cl$frag5))
      if (cl$frag3$n >= min_frag)
        add(.make_met(NA, sprintf("%s Endo 3'frag @%d", oligo_name, i),
                      "endo_3frag", sprintf("endonuclease cleavage after pos %d (3' fragment)", i),
                      i, cl$frag3))
    }
  }
  mets
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
