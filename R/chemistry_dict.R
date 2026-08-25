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
#                     e.g. "Te-sSe-sAe-sSe-...-sGe"  (5' token has no linkage prefix)
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

# Nucleobases (free base). 'S' is the single-letter code for
# 5-methylcytosine (m5C), which appears at every cytosine position in most
# 2'-MOE antisense drugs; 'T' doubles as 5-methyluracil, since the two are
# the same nucleobase.
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

# Sugars (free sugar moiety). 'e' and 'MOE' are the same 2'-MOE sugar --
# 'e' is the single-letter form used in triplet notation, 'MOE' the
# spelled-out form; both resolve to C8H16O6.
# All sugar residues follow one rule: take the parent sugar as its free
# form (ribose C5H10O5, 2'-deoxyribose C5H10O4) and substitute. A 2'-O-R
# ether is ribose with the 2'-OH hydrogen replaced by R, i.e.
# C5H9O5 + R -- which is where MOE (R = C3H7O), NMA (R = C3H6NO),
# allyl (R = C3H5) and aminopropyl (R = C3H8N) all come from. A bridged
# sugar adds the bridge and loses 2 H for the two new bonds (LNA = ribose
# + CH2 - H2; ENA = ribose + C2H4 - H2).
#
# Each residue below is anchored on a nucleoside whose formula is known
# independently -- see test_modifications.R, which rebuilds the nucleoside
# as base + sugar - H2O and checks it.
SUGAR_FORMULAS <- list(
  d   = list(formula = c(C=5,H=10,O=4),           name = "2'-deoxyribose (DNA)"),
  r   = list(formula = c(C=5,H=10,O=5),           name = "ribose (RNA)"),
  m   = list(formula = c(C=6,H=12,O=5),           name = "2'-O-methylribose"),
  f   = list(formula = c(C=5,H=9,F=1,O=4),        name = "2'-fluororibose"),
  e   = list(formula = c(C=8,H=16,O=6),           name = "2'-O-methoxyethylribose (MOE)"),
  MOE = list(formula = c(C=8,H=16,O=6),           name = "2'-O-methoxyethylribose (MOE)"),
  cEt = list(formula = c(C=7,H=12,O=4),           name = "constrained ethyl (cEt)", verify = TRUE),
  LNA = list(formula = c(C=6,H=10,O=5),           name = "locked nucleic acid sugar", verify = TRUE),
  NH2 = list(formula = c(C=5,H=11,N=1,O=4),       name = "2'-amino-2'-deoxyribose"),

  # --- Second-generation and next-generation sugar modifications ---
  # 2'-O-NMA: 2'-O-[2-(methylamino)-2-oxoethyl], the amide-containing
  # MOE successor (Prakash et al., J Med Chem 2008; Ionis). Substituent
  # -CH2-C(=O)-NH-CH3 = C3H6NO, so it is exactly N-H heavier than MOE
  # (+12.9953 Da per residue). MOE-like affinity, better metabolic
  # stability.
  NMA = list(formula = c(C=8,H=15,N=1,O=6),       name = "2'-O-NMA (2'-O-[2-(methylamino)-2-oxoethyl])"),

  # UNA: the C2'-C3' bond of ribose is absent, giving an acyclic
  # (2',3'-seco) residue -- ribose + H2, so +2.0157 Da vs RNA. A duplex
  # destabilizer, used positionally to blunt siRNA off-target effects.
  UNA = list(formula = c(C=5,H=12,O=5),           name = "unlocked nucleic acid (UNA, 2',3'-seco-ribose)"),

  # GNA: the sugar is replaced outright by a three-carbon glycerol unit.
  # Much lighter than ribose (-46.0055 Da).
  GNA = list(formula = c(C=3,H=8,O=3),            name = "glycol nucleic acid (GNA, glycerol backbone)"),

  # ENA: 2'-O,4'-C-ethylene bridge -- LNA's methylene bridge with one
  # more CH2 (+14.0157 Da vs LNA).
  ENA = list(formula = c(C=7,H=12,O=5),           name = "2'-O,4'-C-ethylene bridged (ENA)", verify = TRUE),

  # 2'-O-allyl and 2'-O-aminopropyl: earlier 2'-O-alkyl chemistries,
  # still common in probes and in aptamers.
  allyl = list(formula = c(C=8,H=14,O=5),         name = "2'-O-allylribose"),
  AP    = list(formula = c(C=8,H=17,N=1,O=5),     name = "2'-O-aminopropylribose", verify = TRUE),

  # FANA: 2'-deoxy-2'-fluoro-ARABINOnucleic acid. Same atoms as 2'-F
  # ribo ('f'), opposite 2' stereochemistry -- so it is exactly isobaric
  # and MS cannot tell the two apart. Separate code so a sequence can be
  # recorded correctly even though the mass is identical.
  FANA = list(formula = c(C=5,H=9,F=1,O=4),       name = "2'-deoxy-2'-fluoroarabinose (FANA; isobaric with 'f')"),

  # N3'->P5' phosphoramidate backbones: the 3'-OH is a 3'-NH2, so the
  # change lives in the sugar and the linkage stays an ordinary 'o'/'s'.
  # Deoxyribose with 3'-OH -> 3'-NH2: -O +NH, i.e. -0.9840 Da.
  NP = list(formula = c(C=5,H=11,N=1,O=3),        name = "3'-amino-2',3'-dideoxyribose (N3'->P5' phosphoramidate)")
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
# Every linkage below is the PO bridge with its single non-bridging -OH
# swapped for something else: bridge = HPO3 - OH + X. That is where PS
# (X = SH), methylphosphonate (X = CH3), phosphonoacetate (X = CH2COOH)
# and mesyl phosphoramidate (X = NHSO2CH3) all come from.
LINKAGE_FORMULAS <- list(
  o  = list(formula = c(H=1,P=1,O=3),             name = "phosphodiester (PO)"),
  p  = list(formula = c(H=1,P=1,O=3),             name = "phosphodiester (PO) -- BioPharma Finder notation"),
  s  = list(formula = c(H=1,P=1,S=1,O=2),         name = "phosphorothioate (PS) -- also BioPharma Finder notation"),
  u  = list(formula = c(H=1,P=1,S=1,O=2),         name = "phosphorothioate (PS, stereochemistry variant)"),
  mp = list(formula = c(C=1,H=3,P=1,O=2),         name = "methylphosphonate", verify = TRUE),

  # --- Next-generation backbone chemistries ---
  # Mesyl phosphoramidate ("mu", msPA): -O-P(=O)(NH-SO2-CH3)-O-, the
  # Stetsenko chemistry now used in RNase H ASOs for high nuclease
  # resistance without PS's protein binding. X = NHSO2CH3, so
  # +76.9936 Da vs PO and +61.0164 Da vs PS.
  msp = list(formula = c(C=1,H=4,N=1,O=4,P=1,S=1),
             name = "mesyl phosphoramidate (msPA)", verify = TRUE),

  # Phosphonoacetate (PACE): -O-P(=O)(CH2COOH)-O-. X = CH2COOH,
  # +42.0106 Da vs PO. thioPACE is PACE with the P=O replaced by P=S,
  # a further +15.9772 Da.
  pace  = list(formula = c(C=2,H=3,O=4,P=1),
               name = "phosphonoacetate (PACE)", verify = TRUE),
  tpace = list(formula = c(C=2,H=3,O=3,P=1,S=1),
               name = "thiophosphonoacetate (thioPACE)", verify = TRUE),

  # Phosphoryl guanidine (PGO): X is a guanidine residue in place of the
  # -OH, giving an electrically neutral backbone. Two variants are in
  # common use and they are NOT the same mass -- 'pgo' is the
  # 1,3-dimethylimidazolidin-2-imine (cyclic) group, 'tmg' the
  # tetramethylguanidine (acyclic) group, which is H2 heavier. Check
  # which one your synthesis used.
  pgo = list(formula = c(C=5,H=10,N=3,O=2,P=1),
             name = "phosphoryl guanidine (1,3-dimethylimidazolidin-2-imine variant)",
             verify = TRUE),
  tmg = list(formula = c(C=5,H=12,N=3,O=2,P=1),
             name = "phosphoryl guanidine (tetramethylguanidine, Tmg variant)",
             verify = TRUE),

  # Alkyl phosphonates. Same swap, X = the alkyl group: neutral linkages
  # used at siRNA positions 6-7 from the 5' end to blunt off-target
  # activity, and in ASOs to steer RNase H1 cleavage. 'mp' (methyl,
  # X = CH3) is above; these are the larger members of the same series.
  prp = list(formula = c(C=3,H=7,O=2,P=1),
             name = "propyl phosphonate", verify = TRUE),
  ibu = list(formula = c(C=4,H=9,O=2,P=1),
             name = "isobutyl phosphonate (iBu)", verify = TRUE),
  chx = list(formula = c(C=6,H=11,O=2,P=1),
             name = "cyclohexyl phosphonate (cHex)", verify = TRUE),
  mop = list(formula = c(C=4,H=9,O=3,P=1),
             name = "methoxypropyl phosphonate (MOP)", verify = TRUE)
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
  GalNAc3_triantennary = list(formula = c(C=78,H=140,N=11,O=34,P=1), name = "3'-triantennary GalNAc (BioPharma Finder 'g'; distinct from GalNAc3 above)", attach = "replace_H", verify = TRUE),

  # -- Fatty acid conjugates ------------------------------------------------
  # Lipid conjugation lets a gapmer ride serum albumin, which lengthens
  # circulation time and shifts biodistribution. Chain lengths from C14 to
  # C22 are in use, attached at the 5' end, the 3' end, or internally.
  #
  # Formulas are the FREE fatty acid with attach = "replace_OH", which
  # subtracts H2O -- the ester condensation onto a terminal hydroxyl. The
  # net addition is therefore the acyl group: palmitoyl +238.2297,
  # myristoyl +210.1984, stearoyl +266.2610, docosanoyl +322.3236 Da.
  #
  # If your construct puts a cleavable linker between the lipid and the
  # oligo -- a d(TCA) phosphodiester linker is common -- that linker is
  # not included here. Add its three nucleotides to the sequence itself,
  # or fold its mass into a custom conjugate entry.
  myristoyl  = list(formula = c(C=14,H=28,O=2), name = "myristic acid (C14:0) conjugate",   attach = "replace_OH"),
  palmitoyl  = list(formula = c(C=16,H=32,O=2), name = "palmitic acid (C16:0) conjugate",   attach = "replace_OH"),
  stearoyl   = list(formula = c(C=18,H=36,O=2), name = "stearic acid (C18:0) conjugate",    attach = "replace_OH"),
  docosanoyl = list(formula = c(C=22,H=44,O=2), name = "docosanoic acid (C22:0) conjugate", attach = "replace_OH")
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

# Backward-compatible alias for scripts written against the older name.
REFERENCE_DICT <- STANDARD_DICT

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

  # Linkage codes are not all one character (mp, msp, pace, tpace, pgo),
  # so the prefix is matched against the dictionary's linkage codes
  # longest-first rather than by taking substr(tk, 1, 1). Longest-first
  # matters: "pgoAd" starts with both "p" and "pgo", and only "pgo"
  # leaves a parseable base+sugar behind. A match must leave at least two
  # characters (one base, one sugar) or it is not a prefix at all.
  link_codes <- names(dict)[vapply(dict, function(e)
    identical(e$kind, "linkage"), logical(1))]
  link_codes <- link_codes[order(-nchar(link_codes))]

  for (i in seq_len(n)) {
    tk <- toks[i]
    prefix <- NA_character_
    # The 5' token has no incoming bond, so it never carries a prefix.
    if (i >= 2) {
      for (lc in link_codes) {
        if (startsWith(tk, lc) && nchar(tk) - nchar(lc) >= 2) {
          prefix <- lc
          break
        }
      }
    }
    if (!is.na(prefix)) {
      # prefix = incoming bond (i-1 -> i) = outgoing bond of position i-1
      linkages[i - 1] <- prefix
      tk <- substr(tk, nchar(prefix) + 1, nchar(tk))
    }
    # The sugar is everything after the base letter, not just one
    # character, so multi-character dictionary codes (MOE, cEt, LNA, NMA)
    # work in triplet notation too -- e.g. "sTMOE", "AcEt". Single-letter
    # codes (d, r, m, f, e) are unaffected.
    bases[i] <- substr(tk, 1, 1)
    sugars[i] <- substr(tk, 2, nchar(tk))
  }
  linkages[n] <- NA_character_
  list(bases = bases, sugars = sugars, linkages = linkages)
}


## ---- Published reference examples ------------------------------------------
# Four approved oligonucleotide therapeutics, one per major modality, used as
# worked examples and as the known-mass anchor for the formula engine's
# self-test. They are illustrations of the chemistry classes this pipeline
# handles -- not default targets and not assumed inputs; run the pipeline on
# whatever sequence you actually care about.
#
# Modality classification follows Takakusa et al. (2023), Nucleic Acid
# Therapeutics 33:83-94, Table 1. Sequences and per-position modification
# patterns are from the primary literature and regulatory filings cited on
# each entry (see the Bibliography in README.md).
#
# Duplex drugs (patisiran, givosiran) are given as their SENSE strand: this
# pipeline profiles one strand at a time, so run each strand separately.

# --- Antisense SSO: nusinersen (SPINRAZA) ---
# 18-mer, uniformly 2'-MOE, fully phosphorothioate, 5-methyl on every
# pyrimidine (MeC = 5-methylcytosine -> base code S; MeU = 5-methyluracil,
# i.e. the thymine base -> base code T).
# Sequence 5'-UCACUUUCAUAAUGCUGG-3'.
NUSINERSEN_TRIPLET <- paste0(
  "Te-sSe-sAe-sSe-sTe-sTe-sTe-sSe-sAe-",
  "sTe-sAe-sAe-sTe-sGe-sSe-sTe-sGe-sGe")
# Published free-acid molecular formula (FDA/EMA product labelling; average
# MW 7127.2 Da). One of the formula engine's two regression anchors.
NUSINERSEN_FORMULA <- "C234H340N61O128P17S17"

# --- Antisense gapmer: inotersen (TEGSEDI) ---
# 20-mer 5-10-5 gapmer: 2'-MOE wings (positions 1-5, 16-20), DNA gap
# (6-15), fully phosphorothioate, all cytosines 5-methyl.
# Sequence 5'-TCTTGGTTACATGAAATCCC-3' (TTR 3'UTR, bases 618-637).
INOTERSEN_TRIPLET <- paste0(
  "Te-sSe-sTe-sTe-sGe-",
  "sGd-sTd-sTd-sAd-sSd-sAd-sTd-sGd-sAd-sAd-",
  "sAe-sTe-sSe-sSe-sSe")
# Published free-acid molecular formula (TEGSEDI product labelling; average
# MW 7183.08 Da). The engine's second regression anchor.
INOTERSEN_FORMULA <- "C230H318N69O121P19S19"

# --- siRNA (LNP): patisiran (ONPATTRO), sense strand ---
# 21-mer, mixed ribose / 2'-O-methyl, all-phosphodiester backbone, dTdT
# 3' overhang. Sequence 5'-GUAACCAAGAGUAUUCCAU-dTdT-3'.
# (Antisense strand: 5'-AUGGAAUACUCUUGGUUAC-dTdT-3'.)
PATISIRAN_SENSE_TRIPLET <- paste0(
  "Gr-oUm-oAr-oAr-oCm-oCm-oAr-oAr-oGr-oAr-",
  "oGr-oUm-oAr-oUm-oUm-oCm-oCm-oAr-oUm-oTd-oTd")

# --- siRNA (GalNAc): givosiran (GIVLAARI), sense strand ---
# 21-mer, alternating 2'-O-methyl / 2'-fluoro, phosphorothioate at the two
# 5'-terminal linkages only (remainder phosphodiester), trivalent GalNAc
# (L96) conjugated at the 3' terminus. Structured input is used here because
# triplet notation carries no conjugate field.
# Sequence 5'-CAGAAAGAGUGUCUCAUCUUA-(GalNAc3)-3'.
GIVOSIRAN_SENSE_SPEC <- list(
  bases    = c("C","A","G","A","A","A","G","A","G","U","G",
               "U","C","U","C","A","U","C","U","U","A"),
  sugars   = c("m","m","m","m","m","m","f","m","f","m","f",
               "m","f","m","f","m","m","m","m","m","m"),
  linkages = c("s","s", rep("o", 18), NA),
  conj5    = "none",
  conj3    = "GalNAc3"
)

# Triplet form of the givosiran sense strand, for the notation-based entry
# points (the 3' GalNAc is carried separately in the registry below, since
# triplet notation has no conjugate field).
GIVOSIRAN_SENSE_TRIPLET <- paste0(
  "Cm-sAm-sGm-oAm-oAm-oAm-oGf-oAm-oGf-oUm-",
  "oGf-oUm-oCf-oUm-oCf-oAm-oUm-oCm-oUm-oUm-oAm")

# Registry of the four examples, for the app's "Load example" dropdown and
# for programmatic access. Each entry carries the modality label used by
# Takakusa et al. (2023) Table 1, a triplet string, and the terminal
# conjugates that triplet notation cannot express on its own.
REFERENCE_OLIGOS <- list(
  nusinersen = list(
    name      = "nusinersen",
    brand     = "SPINRAZA",
    modality  = "Antisense (SSO)",
    chemistry = "PS, uniform 2'-MOE, 5-methyl pyrimidines",
    length    = 18,
    triplet   = NUSINERSEN_TRIPLET,
    conj5     = "none",
    conj3     = "none",
    formula   = NUSINERSEN_FORMULA),
  inotersen = list(
    name      = "inotersen",
    brand     = "TEGSEDI",
    modality  = "Antisense (gapmer)",
    chemistry = "PS, 5-10-5 2'-MOE/DNA gapmer, 5-methylcytosine",
    length    = 20,
    triplet   = INOTERSEN_TRIPLET,
    conj5     = "none",
    conj3     = "none",
    formula   = INOTERSEN_FORMULA),
  patisiran = list(
    name      = "patisiran (sense strand)",
    brand     = "ONPATTRO",
    modality  = "siRNA (LNP)",
    chemistry = "2'-OMe / ribose, phosphodiester backbone, dTdT overhang",
    length    = 21,
    triplet   = PATISIRAN_SENSE_TRIPLET,
    conj5     = "none",
    conj3     = "none",
    formula   = NA_character_),
  givosiran = list(
    name      = "givosiran (sense strand)",
    brand     = "GIVLAARI",
    modality  = "siRNA (GalNAc)",
    chemistry = "2'-OMe / 2'-F, partial PS, 3' trivalent GalNAc",
    length    = 21,
    triplet   = GIVOSIRAN_SENSE_TRIPLET,
    conj5     = "none",
    conj3     = "GalNAc3",
    formula   = NA_character_)
)

## ---- Formula engine self-test ----------------------------------------------
# Re-derives the molecular formulas of two approved drugs -- nusinersen (an
# all-MOE SSO) and inotersen (a MOE/DNA gapmer) -- from their sequences and
# chemistry, and compares each against the published free-acid formula on
# its product labelling. Between them they exercise MOE and DNA sugars,
# 5-methylcytosine and 5-methyluracil bases, and a fully phosphorothioate
# backbone.
#
# This is a regression check on the formula engine (run it after editing the
# dictionary) -- it says nothing about whatever sequence you are actually
# analyzing.
.REFERENCE_ANCHORS <- list(
  list(name = "nusinersen", desc = "18-mer PS/2'-MOE SSO",
       triplet = NUSINERSEN_TRIPLET, formula = NUSINERSEN_FORMULA,
       published_avg = 7127.2),
  list(name = "inotersen",  desc = "20-mer PS 5-10-5 MOE/DNA gapmer",
       triplet = INOTERSEN_TRIPLET,  formula = INOTERSEN_FORMULA,
       published_avg = 7183.08)
)

validate_reference <- function(dict = STANDARD_DICT, verbose = TRUE) {
  results <- lapply(.REFERENCE_ANCHORS, function(a) {
    p <- parse_triplet(a$triplet, dict)
    f <- assemble_oligo_formula(p$bases, p$sugars, p$linkages, dict = dict)
    want_vec  <- parse_formula(a$formula)
    want_mass <- formula_mass(want_vec, mono = TRUE)
    got_mass  <- formula_mass(f, mono = TRUE)
    ppm <- if (want_mass > 0) abs(got_mass - want_mass) / want_mass * 1e6 else NA_real_
    # Compare by atom counts (order-independent), not by string.
    ok_f <- all(f == want_vec)
    ok_m <- !is.na(ppm) && ppm < 1
    if (verbose) {
      cat(sprintf("%s (%s)\n", a$name, a$desc),
          "  bases        : ", paste(p$bases, collapse = ""), "\n",
          "  sugars       : ", paste(p$sugars, collapse = "/"), "\n",
          "  formula got  : ", format_formula(f), "\n",
          "  formula want : ", a$formula, "  ", ifelse(ok_f, "OK", "MISMATCH"), "\n",
          "  mono mass    : ", sprintf("%.6f Da", got_mass), "\n",
          "  average mass : ", sprintf("%.2f", formula_mass(f, mono = FALSE)),
          sprintf(" (published %.2f)", a$published_avg), "\n",
          "  delta        : ", sprintf("%.6f Da (%.4f ppm)", got_mass - want_mass, ppm),
          "  ", ifelse(ok_m, "OK (<1 ppm)", "FAIL"), "\n\n", sep = "")
    }
    list(name = a$name, formula = format_formula(f), formula_vec = f,
         mass = got_mass, avg_mass = formula_mass(f, mono = FALSE),
         ppm = ppm, ok = ok_f && ok_m, parsed = p)
  })
  names(results) <- vapply(.REFERENCE_ANCHORS, function(a) a$name, "")
  all_ok <- all(vapply(results, function(r) isTRUE(r$ok), logical(1)))
  if (verbose) {
    cat("Formula-engine self-test: ",
        if (all_ok) "PASS" else "FAIL",
        " (", sum(vapply(results, function(r) isTRUE(r$ok), logical(1))),
        "/", length(results), " published formulas reproduced)\n", sep = "")
  }
  # Top-level fields describe the primary anchor, so callers that expect a
  # single result (workbook validation panel, report) keep working -- but
  # `ok` must report ALL anchors, so drop the primary anchor's own `ok`
  # before appending rather than ending up with two `ok` entries, where
  # `$ok` would silently return the first drug's flag instead of the
  # overall verdict.
  primary <- results[[1]]
  primary$ok <- NULL
  invisible(c(primary, list(ok = all_ok, all = results)))
}
