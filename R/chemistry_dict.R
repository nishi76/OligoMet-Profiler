# =============================================================================
# chemistry_dict.R
# Oligonucleotide metabolite identification -- chemistry dictionary & formula engine
#
# Grounded in:
#   - SynONIM (Lippens et al., JASMS 2024): synthetic oligo modifications/impurities
#   - Eluforsen metabolite profiling (Kim et al., Mol Ther Nucleic Acids 2019)
#   - OligoDistiller (Liu et al., Anal Chem 2025): McLuckey fragments, isotope fit
#   - FMVS automatic metabolite ID (Ye et al., J Chromatogr B 2025)
#
# Notation conventions supported (see oligo_io.R):
#   * Triplet:        [linkage][base][sugar] per token, dash-separated
#                     e.g. "Ge-uAn-sGn-sSn-...-sTn"  (5' token has no linkage prefix)
#                     Also accepts BioPharma Finder-style triplet notation,
#                     e.g. "Ad-pTd-pCd-pAd" ('p' = phosphodiester, 's' =
#                     phosphorothioate -- see LINKAGE_FORMULAS below).
#   * OligoDistiller: OH-<base><sugar>[*]-...-<base><sugar>-OH  (* = PS bond)
#                     e.g. "OH-Am*-Af*-Cm-...-Cm-OH"  (m=2'OMe, f=2'F, d=deoxy)
#   * Structured:     parallel vectors bases/sugars/linkages + conjugate spec
# =============================================================================

## ---- Atomic masses ----------------------------------------------------------
# Monoisotopic (IUPAC) and average atomic weights.
# Elements vector fixes the canonical order used by all formula vectors.
.ELEMENTS <- c("C","H","N","O","P","S","F","Na","K","I","Cl","Br")

.atomic_mass_mono <- c(
  C = 12.000000000, H = 1.007825032, N = 14.003074005, O = 15.994914619,
  P = 30.973761998, S = 31.972071174, F = 18.998403163, Na = 22.989769281,
  K = 38.963706687, I = 126.904473, Cl = 34.968852682, Br = 78.9183376)

.atomic_mass_avg <- c(
  C = 12.011, H = 1.008, N = 14.007, O = 15.999, P = 30.974, S = 32.06,
  F = 18.998, Na = 22.990, K = 39.098, I = 126.904, Cl = 35.45, Br = 79.904)

# Proton / electron masses for charge-state math
.PROTON_MONO <- 1.007276466
.ELECTRON_MONO <- 0.00054858

## ---- Shorthand operator ----------------------------------------------------
# Defined early so every downstream function can rely on it.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## ---- Formula helpers --------------------------------------------------------
# A formula is a named numeric vector over .ELEMENTS (missing elements = 0).
.empty_formula <- function() setNames(rep(0, length(.ELEMENTS)), .ELEMENTS)

# Coerce a compact named vector (e.g. c(C=5,H=5,N=5)) to full element order.
.as_formula <- function(x) {
  f <- .empty_formula()
  if (is.null(x) || length(x) == 0) return(f)
  x <- x[names(x) %in% .ELEMENTS]
  f[names(x)] <- as.numeric(x)
  f
}

# Parse a formula string like "C5H5N5O" into a full element vector.
parse_formula <- function(s) {
  f <- .empty_formula()
  if (is.null(s) || !nzchar(s)) return(f)
  m <- gregexpr("([A-Z][a-z]?)([0-9]*)", s)[[1]]
  if (m[1] < 0) return(f)
  for (i in seq_along(m)) {
    start <- m[i]; end <- start + attr(m, "match.length")[i] - 1
    tok <- substr(s, start, end)
    el <- sub("[0-9]+$", "", tok)
    cnt <- suppressWarnings(as.numeric(sub("^[A-Za-z]+", "", tok)))
    if (is.na(cnt)) cnt <- 1
    if (el %in% .ELEMENTS) f[el] <- f[el] + cnt
  }
  f
}

# Formula -> string (only non-zero elements, canonical order).
format_formula <- function(f) {
  f <- .as_formula(f)
  parts <- vapply(seq_along(f), function(i) {
    if (f[i] == 0) return("")
    nm <- names(f)[i]
    if (f[i] == 1) nm else paste0(nm, as.integer(f[i]))
  }, "")
  paste0(parts[parts != ""], collapse = "")
}

# Add two formulas (coerce to full order first).
add_formulas <- function(a, b) .as_formula(a) + .as_formula(b)

# Scale a formula by an integer (e.g. 14 * PS linkage).
scale_formula <- function(a, k) .as_formula(a) * k

# Formula -> monoisotopic or average mass.
formula_mass <- function(f, mono = TRUE) {
  f <- .as_formula(f)
  m <- if (mono) .atomic_mass_mono else .atomic_mass_avg
  sum(f * m[names(f)])
}

# Water formula (for condensation accounting).
.H2O <- c(H = 2, O = 1)

## ---- Chemistry dictionary --------------------------------------------------
# Each entry: list(formula = <compact>, name = <human readable>, verify = <flag>)
# verify = TRUE marks formulas that are best-estimate and should be confirmed
# against a known mass before relying on them.

# Nucleobases (free base). 'S' is the ION337 code for 5-methylcytosine (m5C),
# reverse-engineered from the ION337 formula and validated to 0 ppm.
BASE_FORMULAS <- list(
  A   = list(formula = c(C=5,H=5,N=5),            name = "adenine"),
  G   = list(formula = c(C=5,H=5,N=5,O=1),        name = "guanine"),
  C   = list(formula = c(C=4,H=5,N=3,O=1),        name = "cytosine"),
  T   = list(formula = c(C=5,H=6,N=2,O=2),        name = "thymine"),
  U   = list(formula = c(C=4,H=4,N=2,O=2),        name = "uracil"),
  S   = list(formula = c(C=5,H=7,N=3,O=1),        name = "5-methylcytosine (m5C)"),
  m5C = list(formula = c(C=5,H=7,N=3,O=1),        name = "5-methylcytosine"),
  hm5C= list(formula = c(C=5,H=7,N=3,O=2),        name = "5-hydroxymethylcytosine"),
  m5U = list(formula = c(C=5,H=6,N=2,O=2),        name = "5-methyluracil (thymine base)"),
  psi = list(formula = c(C=4,H=4,N=2,O=2),        name = "pseudouracil (pseudoU base)"),
  D   = list(formula = c(C=5,H=6,N=6),            name = "2,6-diaminopurine"),
  I   = list(formula = c(C=5,H=4,N=4,O=1),        name = "hypoxanthine (inosine base)")
)

# Sugars (free sugar moiety). Standard codes plus ION337 reverse-engineered codes.
# ION337 'n' (C8H15NO6) and 'e' (C9H19NO8) are N-containing 2'-modified riboses
# (likely 2'-O-aminoalkyl-PEG variants); 'a' (C7H12O7) is the 3'-terminal variant.
# These three were derived exactly from the ION337 formula/variant diff and
# reproduce C196H270N72O105S14P14 + 6193.0418 Da to 0 ppm.
SUGAR_FORMULAS <- list(
  d   = list(formula = c(C=5,H=10,O=4),           name = "2'-deoxyribose (DNA)"),
  r   = list(formula = c(C=5,H=10,O=5),           name = "ribose (RNA)"),
  m   = list(formula = c(C=6,H=12,O=5),           name = "2'-O-methylribose"),
  f   = list(formula = c(C=5,H=9,F=1,O=4),        name = "2'-fluororibose"),
  e   = list(formula = c(C=9,H=19,N=1,O=8),       name = "ION337 5'-terminal 2'-mod sugar (C9H19NO8)"),
  n   = list(formula = c(C=8,H=15,N=1,O=6),       name = "ION337 2'-mod sugar (C8H15NO6)"),
  a   = list(formula = c(C=7,H=12,O=7),           name = "ION337 3'-terminal sugar variant (C7H12O7)"),
  MOE = list(formula = c(C=8,H=16,O=6),           name = "2'-O-methoxyethylribose (standard MOE)"),
  cEt = list(formula = c(C=7,H=12,O=4),           name = "constrained ethyl (cEt)", verify = TRUE),
  LNA = list(formula = c(C=6,H=10,O=5),           name = "locked nucleic acid sugar", verify = TRUE),
  NH2 = list(formula = c(C=5,H=11,N=1,O=4),       name = "2'-amino-2'-deoxyribose")
)

# Internucleotide linkages (bridge residue added per internal bond).
# Oligo = sum(bases) + sum(sugars) + sum(linkages) - (2n-1)*H2O + conjugates.
# PO bridge = HPO3 ; PS bridge = HPSO2 (one non-bridging O -> S, +15.977 Da vs PO).
# 'p' added for BioPharma Finder triplet notation compatibility -- BPF's
# building-block editor uses "p" as the explicit phosphodiester linkage
# prefix (e.g. "Ad-pTd-pCd-pAd"), where this pipeline's own notation
# instead used "o". Same formula, both accepted going forward.
# Source: Thermo BioPharma Finder 5.2 Oligonucleotide Analysis User Guide,
# "Manually create a new oligonucleotide sequence" and "Modification
# notation" topics.
LINKAGE_FORMULAS <- list(
  o  = list(formula = c(H=1,P=1,O=3),             name = "phosphodiester (PO)"),
  p  = list(formula = c(H=1,P=1,O=3),             name = "phosphodiester (PO) -- BioPharma Finder notation"),
  s  = list(formula = c(H=1,P=1,S=1,O=2),         name = "phosphorothioate (PS) -- also BioPharma Finder notation"),
  u  = list(formula = c(H=1,P=1,S=1,O=2),         name = "ION337 linkage (PS, stereochem variant)"),
  mp = list(formula = c(C=1,H=3,P=1,O=2),         name = "methylphosphonate", verify = TRUE)
)

# Terminal conjugates / linkers / caps. These are added as residues; the
# attachment chemistry (replace 5'/3'-H vs -OH) is handled in assemble step via
# the conjugate's 'attach' field: "replace_H" (-H) or "replace_OH" (-H2O) or "add" (+0).
CONJUGATE_FORMULAS <- list(
  none          = list(formula = c(),                       name = "no conjugate",          attach = "add"),
  `5'-phosphate`= list(formula = c(H=1,P=1,O=3),            name = "5'-phosphate cap",      attach = "replace_H"),
  `3'-phosphate`= list(formula = c(H=1,P=1,O=3),            name = "3'-phosphate cap",      attach = "replace_H"),
  `3'-cyclophos`= list(formula = c(H=1,P=1,O=3),            name = "3'-cyclic phosphate",   attach = "replace_OH"),
  GalNAc        = list(formula = c(C=8,H=15,N=1,O=6),        name = "N-acetylgalactosamine (mono)", attach = "replace_H"),
  GalNAc3       = list(formula = c(C=31,H=51,N=3,O=23),      name = "trivalent GalNAc cluster",     attach = "replace_H", verify = TRUE),
  cholesterol   = list(formula = c(C=27,H=46,O=1),           name = "cholesterol",                  attach = "replace_H"),
  C6            = list(formula = c(C=6,H=12,N=1),            name = "aminohexyl (C6) linker",       attach = "replace_H"),
  TEG           = list(formula = c(C=8,H=17,O=4),            name = "tetraethylene glycol (TEG)",   attach = "replace_H"),
  C12           = list(formula = c(C=12,H=24,N=1),           name = "aminododecyl (C12) linker",    attach = "replace_H"),
  FAM           = list(formula = c(C=21,H=12,O=5),           name = "fluorescein (FAM)",            attach = "replace_H", verify = TRUE),
  Cy3           = list(formula = c(C=30,H=37,N=3,O=1),       name = "Cy3 dye",                      attach = "replace_H", verify = TRUE),

  # -- BioPharma Finder terminal modification codes -------------------------
  # From the "Modification notation" topic of the Thermo BioPharma Finder
  # 5.2 Oligonucleotide Analysis User Guide (3'-/5'-terminal modification
  # tables). BPF's own single-letter codes are p, s, g (3'), and p, s, b,
  # a, u, r, c (5'); those letters already mean other things in this
  # dictionary's flat namespace (sugars 'a'/'r', linkage 'u'), so these are
  # stored under distinct, descriptive keys instead -- select them from the
  # conj5/conj3 dropdown (app.R) or opts$conj5/conj3 (CLI) rather than
  # inline in the triplet string. 'p' and 's' terminal caps reuse the
  # existing 5'-/3'-phosphate entries above (same net addition once
  # condensation is accounted for); the rest are new.
  # Formulas below are the BPF-reported whole-group masses; treat as
  # best-estimate (verify = TRUE) until checked against a BPF-computed
  # mass for an actual sequence carrying each modification.
  `5'-thiophosphate` = list(formula = c(H=1,P=1,S=1,O=2),    name = "5'-phosphorothioate cap (BioPharma Finder 's')", attach = "replace_H"),
  `3'-thiophosphate` = list(formula = c(H=1,P=1,S=1,O=2),    name = "3'-phosphorothioate cap (BioPharma Finder 's')", attach = "replace_H"),
  biotin        = list(formula = c(C=10,H=16,N=2,O=3,S=1),   name = "5'-biotin (BioPharma Finder 'b')", attach = "replace_H"),
  cAG_cap       = list(formula = c(C=32,H=44,N=15,O=27,P=5), name = "5'-cAG cap analog (BioPharma Finder 'a')", attach = "replace_H", verify = TRUE),
  cAU_cap       = list(formula = c(C=31,H=43,N=12,O=28,P=5), name = "5'-cAU cap analog (BioPharma Finder 'u')", attach = "replace_H", verify = TRUE),
  ARCA_cap      = list(formula = c(C=22,H=32,N=10,O=21,P=4), name = "5'-ARCA cap analog (BioPharma Finder 'r')", attach = "replace_H", verify = TRUE),
  mCAP          = list(formula = c(C=21,H=30,N=10,O=21,P=4), name = "5'-mCAP cap analog (BioPharma Finder 'c')", attach = "replace_H", verify = TRUE),
  GalNAc3_triantennary = list(formula = c(C=78,H=140,N=11,O=34,P=1), name = "3'-triantennary GalNAc (BioPharma Finder 'g'; distinct from GalNAc3 above)", attach = "replace_H", verify = TRUE)
)

## ---- Dictionary assembly / override ----------------------------------------
# Build a single lookup list combining base/sugar/linkage/conjugate entries,
# each value = list(formula=<full vector>, name, verify, kind).
# User can pass overrides as a named list of compact formulas to replace/extend.
build_dictionary <- function(overrides = list()) {
  d <- list()
  add <- function(tbl, kind) for (nm in names(tbl)) {
    e <- tbl[[nm]]
    d[[nm]] <<- list(
      formula = .as_formula(e$formula),
      name = e$name, verify = isTRUE(e$verify), kind = kind,
      attach = e$attach)
  }
  add(BASE_FORMULAS, "base")
  add(SUGAR_FORMULAS, "sugar")
  add(LINKAGE_FORMULAS, "linkage")
  add(CONJUGATE_FORMULAS, "conjugate")
  # apply user overrides (compact formula or full vector)
  for (nm in names(overrides)) {
    ov <- overrides[[nm]]
    if (is.character(ov)) ov <- parse_formula(ov)
    if (is.list(ov)) {
      d[[nm]] <- list(formula = .as_formula(ov$formula), name = ov$name %||% nm,
                      verify = isTRUE(ov$verify), kind = ov$kind %||% d[[nm]]$kind,
                      attach = ov$attach %||% d[[nm]]$attach)
    } else {
      kind <- if (!is.null(d[[nm]]$kind)) d[[nm]]$kind else "custom"
      d[[nm]] <- list(formula = .as_formula(ov), name = nm, verify = FALSE,
                      kind = kind, attach = d[[nm]]$attach)
    }
  }
  d
}

# The default chemistry dictionary (standard bases/sugars/linkages/conjugates,
# no overrides). This is the dictionary every module in the pipeline falls
# back to when the caller doesn't supply one -- it covers DNA, RNA, 2'OMe,
# 2'F, MOE, cEt, LNA sugars, PO/PS linkages, and common conjugates, and is
# not specific to any one oligonucleotide. Build a project-specific
# dictionary instead with build_dictionary(overrides = list(...)) -- see
# run_custom_oligo.R.
STANDARD_DICT <- build_dictionary()

# Backward-compatible alias. The name predates this dictionary being used as
# the general-purpose default for every sequence, not just the ION337
# reference case -- kept so existing scripts/tests referencing ION337_DICT
# keep working unchanged.
ION337_DICT <- STANDARD_DICT

## ---- Oligo formula assembly -------------------------------------------------
# Given canonical per-position vectors (bases, sugars, linkages) and optional
# 5'/3' conjugate codes, assemble the neutral molecular formula.
#   bases, sugars: length n (1..n, 5'->3')
#   linkages: length n; linkages[i] = bond between pos i and i+1; NA at 3' end (pos n)
#   conj5, conj3: conjugate codes (default "none")
# Returns a full element vector.
assemble_oligo_formula <- function(bases, sugars, linkages, conj5 = "none",
                                    conj3 = "none", dict = STANDARD_DICT) {
  n <- length(bases)
  stopifnot(length(sugars) == n, length(linkages) == n)
  f <- .empty_formula()
  for (i in seq_len(n)) {
    f <- add_formulas(f, dict[[bases[i]]]$formula)
    f <- add_formulas(f, dict[[sugars[i]]]$formula)
  }
  # internal linkages (i = 1..n-1)
  for (i in seq_len(n - 1)) {
    lk <- linkages[i]
    if (!is.na(lk) && nzchar(lk)) f <- add_formulas(f, dict[[lk]]$formula)
  }
  # condensation: n glycosidic bonds + (n-1) phosphodiester bonds = (2n-1) H2O
  f <- add_formulas(f, scale_formula(.H2O, -(2 * n - 1)))
  # terminal conjugates
  f <- attach_conjugate(f, conj5, dict)
  f <- attach_conjugate(f, conj3, dict)
  f
}

attach_conjugate <- function(f, code, dict) {
  if (is.null(code) || is.na(code) || code == "none" || !nzchar(code)) return(f)
  e <- dict[[code]]
  if (is.null(e)) stop("Unknown conjugate code: ", code)
  f <- add_formulas(f, e$formula)
  att <- e$attach %||% "add"
  if (att == "replace_H")     f <- add_formulas(f, c(H = -1))
  else if (att == "replace_OH") f <- add_formulas(f, c(H = -2, O = -1))
  f
}

## ---- parse_triplet (general-purpose; works with any dictionary) -----------
# Parse a triplet-notation string into canonical vectors. Used by any
# sequence, not just the reference example below -- also called from
# oligo_io.R.
# Convention: linkages[i] = bond from position i to position i+1 (outgoing);
#             linkages[n] = NA (3' terminal has no outgoing bond).
# In the triplet string the linkage prefix on token i is the INCOMING bond
# (from i-1 to i), so it maps to linkages[i-1].
parse_triplet <- function(triplet, dict = STANDARD_DICT) {
  toks <- trimws(strsplit(triplet, "-")[[1]])
  toks <- toks[nzchar(toks)]
  n <- length(toks)
  bases <- sugars <- character(n)
  linkages <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    tk <- toks[i]
    ch1 <- substr(tk, 1, 1)
    if (nchar(tk) >= 3 && ch1 %in% names(dict) &&
        dict[[ch1]]$kind == "linkage") {
      # prefix = incoming bond (i-1 -> i) = outgoing bond of position i-1
      if (i >= 2) linkages[i - 1] <- ch1
      base <- substr(tk, 2, 2); sugar <- substr(tk, 3, 3)
    } else {
      base <- substr(tk, 1, 1); sugar <- substr(tk, 2, 2)
    }
    bases[i] <- base; sugars[i] <- sugar
  }
  linkages[n] <- NA_character_
  list(bases = bases, sugars = sugars, linkages = linkages)
}

## ---- Reference example / self-test -----------------------------------------
# ION337 is a published gapmer ASO used here purely as a known-mass reference
# to sanity-check the formula engine above -- it is one example sequence
# among any the pipeline can run, not a required input or a default target.
# See run_custom_oligo.R (repository root) to run the pipeline on your own
# sequence; validate_ION337() below re-derives this reference case whenever
# a regression check on the formula engine is needed.
ION337_TRIPLET <- "Ge-uAn-sGn-sSn-sAn-sAn-sGn-sAn-sTn-sTn-sAn-sTn-sSn-sSn-sTn"
ION337_TARGET_FORMULA <- "C196H270N72O105S14P14"
ION337_TARGET_MASS    <- 6193.0417731619

validate_ION337 <- function(dict = STANDARD_DICT, verbose = TRUE) {
  p <- parse_triplet(ION337_TRIPLET, dict)
  f <- assemble_oligo_formula(p$bases, p$sugars, p$linkages, dict = dict)
  got_formula <- format_formula(f)
  got_mass <- formula_mass(f, mono = TRUE)
  ppm <- abs(got_mass - ION337_TARGET_MASS) / ION337_TARGET_MASS * 1e6
  # Compare by atom counts (order-independent), not by string.
  want_vec <- parse_formula(ION337_TARGET_FORMULA)
  ok_f <- all(f == want_vec)
  ok_m <- ppm < 1
  if (verbose) {
    cat("ION337 reference self-test (sanity-checks the formula engine;\n",
        "does not validate whatever sequence you're actually running)\n",
        "  parsed bases : ", paste(p$bases, collapse = ""), "\n",
        "  parsed sugars: ", paste(p$sugars, collapse = ""), "\n",
        "  parsed links : ", paste(ifelse(is.na(p$linkages), ".", p$linkages), collapse = ""), "\n",
        "  formula got : ", got_formula, "\n",
        "  formula want: ", ION337_TARGET_FORMULA, "  ", ifelse(ok_f, "OK", "MISMATCH"), "\n",
        "  mono mass   : ", sprintf("%.6f", got_mass), "\n",
        "  want mass   : ", sprintf("%.6f", ION337_TARGET_MASS), "\n",
        "  delta       : ", sprintf("%.6f Da (%.4f ppm)", got_mass - ION337_TARGET_MASS, ppm),
        "  ", ifelse(ok_m, "OK (<1 ppm)", "FAIL"), "\n", sep = "")
  }
  invisible(list(formula = got_formula, formula_vec = f, mass = got_mass,
                 ppm = ppm, ok = ok_f && ok_m, parsed = p))
}
