# =============================================================================
# oligo_io.R
# Input parsing for oligonucleotide metabolite identification.
#
# Supports three notations, all reduced to a canonical `oligo_spec`:
#   bases[i]    5'->3' nucleobase code   (length n)
#   sugars[i]   5'->3' sugar code        (length n)
#   linkages[i] bond from pos i to i+1   (length n; NA at 3' end)
#   conj5/conj3 terminal conjugate codes
#
# Notations:
#   1. Triplet        "Ge-uAn-sGn-sSn-...-sTn"
#      [linkage][base][sugar] per token; 5' token has no linkage prefix.
#      Linkage prefix on token i = INCOMING bond (i-1 -> i).
#   2. OligoDistiller "OH-Am*-Af*-Cm-...-Cm-OH"
#      <base><sugar>[*] per token; * = PS outgoing bond, else PO.
#      Sugar suffix: m=2'OMe, f=2'F, d=deoxy, r=RNA, e=MOE (configurable).
#   3. Structured     list(bases=, sugars=, linkages=, conj5=, conj3=)
# =============================================================================

## ---- Sugar suffix map for OligoDistiller notation --------------------------
# Maps the single-letter suffix to the dictionary sugar code.
# Users can extend via parse_input(..., sugar_map = ...).
.default_sugar_map <- c(
  m = "m",   # 2'-O-methyl
  f = "f",   # 2'-fluoro
  d = "d",   # 2'-deoxy (DNA)
  r = "r",   # ribose (RNA)
  e = "MOE"  # 2'-O-methoxyethyl
)

## ---- Triplet parser (re-exported from chemistry_dict.R) --------------------
# parse_triplet() is defined in chemistry_dict.R and sourced before this file.

## ---- OligoDistiller parser -------------------------------------------------
parse_oligodistiller <- function(seq, dict = STANDARD_DICT,
                                 sugar_map = .default_sugar_map) {
  s <- trimws(seq)
  if (!grepl("^OH-", s)) stop("OligoDistiller sequence must start with 'OH-'")
  if (!grepl("-OH$", s)) stop("OligoDistiller sequence must end with '-OH'")
  s <- sub("^OH-", "", s)
  s <- sub("-OH$", "", s)
  toks <- strsplit(s, "-")[[1]]
  toks <- toks[nzchar(toks)]
  n <- length(toks)
  bases <- sugars <- character(n)
  linkages <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    tk <- toks[i]
    ps <- endsWith(tk, "*")
    if (ps) tk <- sub("\\*$", "", tk)
    if (nchar(tk) < 2)
      stop("Malformed OligoDistiller token '", toks[i], "' at position ", i)
    base <- substr(tk, 1, 1)
    suf  <- substr(tk, 2, nchar(tk))   # allow multi-char suffix if needed
    if (!(base %in% names(dict) && dict[[base]]$kind == "base"))
      stop("Unknown base '", base, "' at position ", i)
    if (!(suf %in% names(sugar_map)))
      stop("Unknown sugar suffix '", suf, "' at position ", i,
           " (extend via sugar_map=)")
    bases[i] <- base
    sugars[i] <- sugar_map[[suf]]
    # * = PS outgoing bond from this position; else PO
    if (i < n) linkages[i] <- if (ps) "s" else "o"
  }
  linkages[n] <- NA_character_
  list(bases = bases, sugars = sugars, linkages = linkages)
}

## ---- Structured input parser ----------------------------------------------
parse_structured <- function(spec, dict = STANDARD_DICT) {
  need <- c("bases", "sugars", "linkages")
  miss <- setdiff(need, names(spec))
  if (length(miss)) stop("Structured input missing fields: ",
                         paste(miss, collapse = ", "))
  bases <- as.character(spec$bases)
  sugars <- as.character(spec$sugars)
  linkages <- as.character(spec$linkages)
  n <- length(bases)
  stopifnot(length(sugars) == n, length(linkages) == n)
  linkages[is.na(linkages) | linkages == "NA"] <- NA_character_
  conj5 <- spec$conj5 %||% "none"
  conj3 <- spec$conj3 %||% "none"
  list(bases = bases, sugars = sugars, linkages = linkages,
       conj5 = conj5, conj3 = conj3)
}

## ---- Three-line (bases / sugars / linkages) parser -------------------------
# The layout a chemical analysis file or a BioPharma Finder sequence entry
# gives you directly:
#   bases     "TSASTTTSATAATGSTGG"   one code per position, 5'->3'
#   sugars    "eeeeeeeeeeeeeeeeee"   one code per position
#   linkages  "sssssssssssssssss"    one code per *bond*, so n-1 of them
#
# Tokenizing: a string containing a separator (space, comma, dash, slash or
# pipe) is split on it, which is how multi-character codes are written --
# "A,G,S" or "MOE-MOE-d". Otherwise the string is split into single
# characters, the compact form above.
#
# Recycling: a sugar or linkage field holding exactly one code is applied to
# every position ("e" = uniform MOE, "s" = fully phosphorothioate). This is
# the common case for antisense drugs and saves typing out 18 identical
# letters.
#
# Linkage length: n-1 is the correct count (bonds between positions). A
# vector of n is also accepted -- the trailing entry is the 3'-terminal
# "bond", which does not exist and is dropped -- because sequence tables
# often carry a placeholder there.
.tokenize_codes <- function(x, field) {
  if (length(x) > 1) return(as.character(x))
  s <- trimws(as.character(x))
  if (!nzchar(s)) stop("The ", field, " field is empty")
  if (grepl("[ ,/|-]", s)) {
    toks <- strsplit(s, "[ ,/|-]+")[[1]]
    toks <- toks[nzchar(toks)]
  } else {
    toks <- strsplit(s, "")[[1]]
  }
  toks
}

parse_three_line <- function(bases, sugars, linkages,
                             conj5 = "none", conj3 = "none",
                             dict = STANDARD_DICT) {
  b <- .tokenize_codes(bases, "bases")
  s <- .tokenize_codes(sugars, "sugars")
  l <- .tokenize_codes(linkages, "linkages")
  n <- length(b)
  if (n < 2) stop("Need at least 2 bases; got ", n)

  if (length(s) == 1) s <- rep(s, n)
  if (length(s) != n)
    stop("Got ", n, " bases but ", length(s), " sugars. There must be one ",
         "sugar per base (or a single sugar code to apply to all of them).")

  if (length(l) == 1) l <- rep(l, n - 1)
  if (length(l) == n) l <- l[-n]          # drop the non-existent 3' bond
  if (length(l) != n - 1)
    stop("Got ", n, " bases but ", length(l), " linkages. Linkages sit ",
         "*between* positions, so a sequence of ", n, " bases has ", n - 1,
         " of them (or give a single linkage code to apply to all of them).")

  parse_input(list(bases = b, sugars = s,
                   linkages = c(l, NA_character_),
                   conj5 = conj5, conj3 = conj3),
              dict = dict, notation = "structured")
}

## ---- Notation auto-detection + dispatcher ----------------------------------
# Returns a canonical oligo_spec list.
parse_input <- function(x, dict = STANDARD_DICT, sugar_map = .default_sugar_map,
                        notation = c("auto", "triplet", "oligodistiller",
                                     "structured")) {
  notation <- match.arg(notation)

  if (notation == "structured" || is.list(x)) {
    p <- parse_structured(x, dict)
    used <- "structured"
  } else if (notation == "oligodistiller" ||
             (notation == "auto" && grepl("^OH-", trimws(x)))) {
    p <- parse_oligodistiller(x, dict, sugar_map)
    used <- "oligodistiller"
  } else {
    # triplet (default for any other dash-separated string)
    p <- parse_triplet(x, dict)
    used <- "triplet"
  }

  n <- length(p$bases)
  spec <- list(
    bases    = p$bases,
    sugars   = p$sugars,
    linkages = p$linkages,
    conj5    = p$conj5 %||% "none",
    conj3    = p$conj3 %||% "none",
    n        = n,
    notation = used,
    raw      = if (is.character(x) && length(x) == 1) x else deparse(substitute(x))
  )
  validate_spec(spec, dict)
  spec
}

## ---- Spec validation -------------------------------------------------------
# Check every base/sugar/linkage/conjugate code exists in the dictionary.
validate_spec <- function(spec, dict = STANDARD_DICT) {
  chk <- function(codes, kind, ctx) {
    for (cd in unique(codes[!is.na(codes)])) {
      e <- dict[[cd]]
      if (is.null(e))
        stop("Unknown ", kind, " code '", cd, "' in ", ctx)
      if (!is.null(e$kind) && e$kind != kind)
        stop("Code '", cd, "' is a ", e$kind, " but used as ", kind,
             " in ", ctx)
    }
  }
  chk(spec$bases, "base", "bases")
  chk(spec$sugars, "sugar", "sugars")
  chk(spec$linkages, "linkage", "linkages")
  for (cc in c(spec$conj5, spec$conj3)) {
    if (!is.na(cc) && nzchar(cc) && cc != "none") {
      e <- dict[[cc]]
      if (is.null(e)) stop("Unknown conjugate code '", cc, "'")
      if (!is.null(e$kind) && e$kind != "conjugate")
        stop("Code '", cc, "' is a ", e$kind, " but used as conjugate")
    }
  }
  invisible(TRUE)
}

## ---- Renderers (round-trip back to notation) -------------------------------
# Render an oligo_spec to triplet notation.
format_triplet <- function(spec) {
  n <- spec$n
  toks <- character(n)
  for (i in seq_len(n)) {
    core <- paste0(spec$bases[i], spec$sugars[i])
    if (i == 1) {
      toks[i] <- core                       # 5' terminal: no linkage prefix
    } else {
      lk <- spec$linkages[i - 1]            # incoming bond = outgoing of i-1
      toks[i] <- if (is.na(lk)) core else paste0(lk, core)
    }
  }
  paste(toks, collapse = "-")
}

# Render an oligo_spec to OligoDistiller notation.
format_oligodistiller <- function(spec, sugar_map = .default_sugar_map) {
  n <- spec$n
  # invert sugar_map: dict code -> suffix (single-bracket lookup -> NA if absent)
  inv <- setNames(names(sugar_map), unname(sugar_map))
  toks <- character(n)
  for (i in seq_len(n)) {
    suf <- inv[spec$sugars[i]]
    if (is.na(suf)) suf <- spec$sugars[i]   # fall back to raw code
    star <- if (i < n && !is.na(spec$linkages[i]) && spec$linkages[i] == "s") "*" else ""
    toks[i] <- paste0(spec$bases[i], suf, star)
  }
  paste0("OH-", paste(toks, collapse = "-"), "-OH")
}

# Render an oligo_spec as the three plain lines a sequence table (or a
# BioPharma Finder sequence entry form) is filled in with. Returns a named
# character vector: bases (n codes), sugars (n codes), linkages (n-1 codes).
# Multi-character codes are dash-separated so the string round-trips through
# parse_three_line(); single-character codes are written compactly.
format_three_line <- function(spec) {
  pack <- function(v) {
    if (any(nchar(v) > 1)) paste(v, collapse = "-") else paste(v, collapse = "")
  }
  lk <- spec$linkages[seq_len(spec$n - 1)]
  c(bases = pack(spec$bases),
    sugars = pack(spec$sugars),
    linkages = pack(lk))
}

# Render an oligo_spec as a FASTA record for BioPharma Finder's Sequence
# Manager ("Import FASTA File"). The sequence line is triplet notation,
# which is the same 5'-1'-2' building-block notation BPF itself uses --
# [linkage][base][sugar] per residue, 'p' for phosphodiester and 's' for
# phosphorothioate (see the LINKAGE_FORMULAS note in chemistry_dict.R).
#
# Two things to check before importing, because they depend on your BPF
# build rather than on this package:
#   * Sugar codes. This writes the dictionary's own codes (d/r/m/f/e, or
#     MOE/cEt/LNA spelled out). BPF resolves sugar codes against its
#     Building Block editor, so a modification your site has defined under a
#     different code needs renaming here or adding there.
#   * Terminal conjugates. Triplet notation has nowhere to put them, so a
#     5'/3' conjugate on the spec is written into the header as a comment
#     and must be set on the BPF side as a terminus modification.
#
# po_code = "p" writes BPF's phosphodiester prefix; set it to "o" to keep
# this package's own notation instead.
format_biopharma_fasta <- function(spec, name = "oligo", po_code = "p") {
  s <- spec
  s$linkages <- ifelse(!is.na(s$linkages) & s$linkages == "o",
                       po_code, s$linkages)
  conj <- character(0)
  if (!is.null(s$conj5) && nzchar(s$conj5) && s$conj5 != "none")
    conj <- c(conj, paste0("5'-", s$conj5))
  if (!is.null(s$conj3) && nzchar(s$conj3) && s$conj3 != "none")
    conj <- c(conj, paste0("3'-", s$conj3))
  header <- paste0(">", name, " | ", s$n, "-mer",
                   if (length(conj))
                     paste0(" | set as terminus modifications in BioPharma ",
                            "Finder: ", paste(conj, collapse = ", ")) else "")
  paste0(header, "\n", format_triplet(s), "\n")
}

# Human-readable one-line summary.
format_spec <- function(spec) {
  paste0("[", spec$notation, "] n=", spec$n,
         "  bases=", paste(spec$bases, collapse = ""),
         "  sugars=", paste(spec$sugars, collapse = ""),
         "  links=", paste(ifelse(is.na(spec$linkages), ".", spec$linkages), collapse = ""),
         "  conj5=", spec$conj5, "  conj3=", spec$conj3)
}
