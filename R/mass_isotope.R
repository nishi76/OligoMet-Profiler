# =============================================================================
# mass_isotope.R
# Mass, charge-envelope, isotope, and biotransformation calculations for
# oligonucleotide metabolite identification.
#
# Conventions:
#   * Negative ESI:  [M - zH]^z-   =>  m/z = (M + h_offset - z*mp) / z
#     h_offset = 0  (standard, matches eluforsen/FMVS papers)
#     h_offset = 3  (reproduces a legacy workbook's +3H neutral convention)
#   * PS->PO oxidation: replace 1 S with 1 O per event  =>  -15.977156 Da
#   * Adducts (negative mode, 1 cation substituting 1 H):
#       +Na:  +21.981944 Da    +K:  +37.955882 Da    +NH4: +17.026549 Da
#   * Depurination: loss of free base + H2O regain  =>  -base + H2O
# =============================================================================

## ---- Constants -------------------------------------------------------------
.PROTON <- 1.007276466          # proton mass (Da)
.H2O_MASS <- 18.010564686       # monoisotopic water
.PS_TO_PO_SHIFT <- .atomic_mass_mono[["S"]] - .atomic_mass_mono[["O"]]  # 15.977157

# Some legacy reference workbooks compute their charge envelope from a
# neutral mass of [M] + 3.0046 Da (not 3*proton or 3*H_atom -- a
# spreadsheet-specific quirk, consistent to +/-0.00002 Da across z=3..12).
# Default h_offset = 0 gives the standard [M-zH]^z- convention
# (eluforsen/FMVS). Set h_offset to this value only when reproducing m/z
# values from a workbook that used that convention.
.LEGACY_WB_OFFSET <- 3.0046

# Adduct mass shifts (cation replaces one H on the backbone)
.ADDUCT_SHIFT <- c(
  Na  = .atomic_mass_mono[["Na"]] - .atomic_mass_mono[["H"]],   # +21.981944
  K   = .atomic_mass_mono[["K"]]  - .atomic_mass_mono[["H"]],   # +37.955882
  NH4 = .atomic_mass_mono[["N"]]  + 4 * .atomic_mass_mono[["H"]] - .atomic_mass_mono[["H"]]  # +17.026549 (N+3H net)
)

## ---- Metabolite mass info --------------------------------------------------
#' Compute a metabolite's mass and molecular formula
#'
#' @param met A metabolite object (see [generate_metabolites()]).
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return A list: `mono_mass`, `avg_mass`, `formula_str`, `formula_vec`
#'   (the full element vector, e.g. for [ps_oxid_formula()]).
#' @seealso [assemble_oligo_formula()], [charge_envelope()],
#'   [ps_oxidation_series()]
#' @export
metabolite_mass_info <- function(met, dict = STANDARD_DICT) {
  f <- assemble_oligo_formula(met$bases, met$sugars, met$linkages,
                              met$conj5, met$conj3, dict = dict)
  list(
    mono_mass   = formula_mass(f, mono = TRUE),
    avg_mass    = formula_mass(f, mono = FALSE),
    formula_str = format_formula(f),
    formula_vec = f
  )
}

## ---- Charge envelope -------------------------------------------------------
#' Compute the negative-ESI charge envelope for a neutral mass
#'
#' `[M - zH]^z-` convention: `m/z = (M + h_offset - z*proton) / z`.
#'
#' @param mono_mass Monoisotopic neutral mass, in Da.
#' @param z_range Charge states to compute.
#' @param h_offset `0` for the standard `[M-zH]^z-` convention
#'   (matches the eluforsen/FMVS papers); non-zero (e.g. `3.0046`) only
#'   to reproduce a legacy workbook's own neutral-mass convention.
#' @return A data.frame(`z`, `mz`): the monoisotopic peak's m/z at each
#'   charge state.
#' @seealso [isotope_mz_cluster()], [match_ms1()]
#' @export
charge_envelope <- function(mono_mass, z_range = 3:12, h_offset = 0) {
  z <- z_range[z_range >= 1]
  mz <- (mono_mass + h_offset - z * .PROTON) / z
  data.frame(z = z, mz = mz, stringsAsFactors = FALSE)
}

## ---- PS->PO oxidation series ----------------------------------------------
#' Mass after k phosphorothioate-to-phosphodiester desulfurizations
#'
#' Each PS -> PO oxidation event replaces one sulfur with one oxygen
#' (-15.977157 Da).
#'
#' @param mono_mass Monoisotopic mass before oxidation, in Da.
#' @param k Number of desulfurization events.
#' @return `mono_mass - k * 15.977157`.
#' @seealso [ps_oxid_formula()], [ps_oxidation_series()]
#' @export
ps_oxid_mass <- function(mono_mass, k) mono_mass - k * .PS_TO_PO_SHIFT

#' Molecular formula after k phosphorothioate-to-phosphodiester oxidations
#'
#' @param formula_vec A full element formula vector (see
#'   [metabolite_mass_info()]).
#' @param k Number of desulfurization events (subtracted from S, added
#'   to O).
#' @return The formula vector after `k` oxidations.
#' @seealso [ps_oxid_mass()], [ps_oxidation_series()]
#' @export
ps_oxid_formula <- function(formula_vec, k) {
  f <- formula_vec
  f[["S"]] <- f[["S"]] - k
  f[["O"]] <- f[["O"]] + k
  f
}

#' Compute a metabolite's full PS -> PO oxidation series
#'
#' @param met A metabolite object (see [generate_metabolites()]).
#' @param max_oxid Maximum oxidation events to model (capped at `met`'s
#'   own phosphorothioate count).
#' @param z_range,h_offset Unused by the current implementation (reserved
#'   for an envelope column); kept for signature compatibility.
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return A data.frame, one row per oxidation level `k = 0` (unoxidized)
#'   through `min(met$n_ps, max_oxid)`: `k`, `mono_mass`, `avg_mass`,
#'   `formula_str`.
#' @seealso [ps_oxid_mass()], [ps_oxid_formula()]
#' @export
ps_oxidation_series <- function(met, max_oxid = 6, z_range = 3:12,
                                h_offset = 0, dict = STANDARD_DICT) {
  info <- metabolite_mass_info(met, dict)
  kmax <- min(met$n_ps, max_oxid)
  rows <- list()
  for (k in 0:kmax) {
    m <- ps_oxid_mass(info$mono_mass, k)
    fv <- ps_oxid_formula(info$formula_vec, k)
    rows[[k + 1]] <- data.frame(
      k = k, mono_mass = m, avg_mass = formula_mass(fv, mono = FALSE),
      formula_str = format_formula(fv), stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

## ---- Isotope pattern calculation ------------------------------------------
# Uses enviPat if available; falls back to built-in convolution.
# Returns data.frame(mass, abundance) sorted by mass, truncated to top-N.
# Held in an environment (not a plain top-level variable) so it can be
# filled lazily even when this file is loaded as a sealed package namespace,
# where `<<-` on a top-level binding would fail with a locked-binding error.
.isotope_data_env <- new.env(parent = emptyenv())

# In-session memoization cache. isotope_mz_cluster() is called once per
# charge state per adduct (see match_ms1()), but the underlying isotope
# pattern for a given formula does not depend on charge or adduct at all --
# only the m/z shift applied afterward does. Without caching, the same
# expensive enviPat::isopattern() convolution was being recomputed up to
# ~10-40x per metabolite (once per z in z_range, per matching adduct),
# which is what made MS1 matching appear to hang on larger oligos.
.envipat_cache <- new.env(parent = emptyenv())

.isotope_pattern_envipat <- function(formula_str, threshold = 1e-6, n_top = 15) {
  if (is.null(.isotope_data_env$isotopes)) {
    if (!requireNamespace("enviPat", quietly = TRUE)) return(NULL)
    e <- new.env()
    data(isotopes, package = "enviPat", envir = e)
    .isotope_data_env$isotopes <- e$isotopes
  }

  cache_key <- paste(formula_str, threshold, n_top, sep = "|")
  cached <- .envipat_cache[[cache_key]]
  if (!is.null(cached)) return(cached)

  # Do not pass emass/algo explicitly: enviPat's isopattern() prints
  # "NOTE: You are sure that is the mass of an electrone?" whenever emass
  # is supplied as an explicit argument at all (regardless of its value,
  # since 0.00054858 is already enviPat's own default) -- so the previous
  # code triggered that NOTE on every single call. Letting enviPat fall
  # back to its own defaults avoids the check entirely. sink() (both
  # output and message streams) is kept as a second layer of defense for
  # the unrelated "done." status print.
  tmp <- tempfile()
  con <- file(tmp, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  on.exit({
    sink(type = "message"); sink(type = "output")
    close(con); unlink(tmp)
  }, add = TRUE)
  res <- try(
    suppressMessages(suppressWarnings(
      enviPat::isopattern(.isotope_data_env$isotopes, formula_str,
                           threshold = threshold, plotit = FALSE,
                           charge = FALSE)
    )),
    silent = TRUE
  )

  if (is.null(res) || inherits(res, "try-error")) return(NULL)
  pat <- res[[1]]
  if (is.null(pat) || nrow(pat) == 0) return(NULL)
  out <- data.frame(mass = pat[, 1], abundance = pat[, 2] / 100,
                    stringsAsFactors = FALSE)
  out <- out[order(out$mass), ]
  top <- head(out[order(-out$abundance), ], n_top)
  top <- top[order(top$mass), ]

  .envipat_cache[[cache_key]] <- top
  top
}

# Built-in isotope table for fallback (monoisotopic + major isotopes).
.BUILTIN_ISOTOPES <- list(
  C  = list(masses = c(12, 13.003354838), abund = c(0.9893, 0.0107)),
  H  = list(masses = c(1.007825032, 2.014101778), abund = c(0.999885, 0.000115)),
  N  = list(masses = c(14.003074005, 15.000108898), abund = c(0.99636, 0.00364)),
  O  = list(masses = c(15.994914619, 16.999131756, 17.999159613),
            abund = c(0.99757, 0.00038, 0.00205)),
  P  = list(masses = c(30.973761998), abund = c(1.0)),
  S  = list(masses = c(31.972071174, 32.971458762, 33.967867004, 35.967080880),
            abund = c(0.9485, 0.00763, 0.04293, 0.00002)),
  F  = list(masses = c(18.998403163), abund = c(1.0)),
  Na = list(masses = c(22.989769281), abund = c(1.0)),
  K  = list(masses = c(38.963706687, 40.961825257), abund = c(0.9326, 0.0673)),
  I  = list(masses = c(126.904473), abund = c(1.0)),
  Cl = list(masses = c(34.968852682, 36.965903260), abund = c(0.7576, 0.2424)),
  Br = list(masses = c(78.9183376, 80.9162906), abund = c(0.5069, 0.4931))
)

# Convolve two isotope patterns (mass-add, abundance-multiply), truncate.
# Keeps only top-N peaks after aggregation to prevent memory blowup for large
# formulas (e.g. C196H270... where naive convolution would produce millions
# of intermediate peaks).
.convolve_patterns <- function(p1, p2, threshold = 1e-7, max_peaks = 200) {
  masses <- outer(p1$mass, p2$mass, "+")
  abund  <- outer(p1$abundance, p2$abundance, "*")
  dim(masses) <- NULL; dim(abund) <- NULL
  # aggregate by rounded mass (0.0001 Da bins)
  bins <- round(masses, 4)
  out <- tapply(abund, bins, sum)
  df <- data.frame(mass = as.numeric(names(out)), abundance = as.numeric(out),
                   stringsAsFactors = FALSE)
  # Truncate to top-N by abundance to keep pattern size bounded
  df <- df[df$abundance >= threshold, ]
  df <- df[order(-df$abundance), ]
  if (nrow(df) > max_peaks) df <- head(df, max_peaks)
  df
}

# Memoization cache for the built-in convolution engine, mirroring
# .envipat_cache above. Without it, the fallback path recomputed the same
# full convolution once per charge state per adduct -- the dominant cost of
# workbook/report builds on systems without enviPat installed.
.builtin_iso_cache <- new.env(parent = emptyenv())

.isotope_pattern_builtin <- function(formula_vec, threshold = 1e-7, n_top = 15) {
  f <- .as_formula(formula_vec)
  cache_key <- paste(format_formula(f), threshold, n_top, sep = "|")
  cached <- .builtin_iso_cache[[cache_key]]
  if (!is.null(cached)) return(cached)
  pat <- data.frame(mass = 0, abundance = 1, stringsAsFactors = FALSE)
  for (el in names(f)) {
    cnt <- f[el]
    if (cnt == 0) next
    iso <- .BUILTIN_ISOTOPES[[el]]
    if (is.null(iso)) next
    # pattern for one atom
    one <- data.frame(mass = iso$masses, abundance = iso$abund,
                      stringsAsFactors = FALSE)
    # Raise the one-atom pattern to the cnt-th convolution power by binary
    # exponentiation: multiply `pat` by elpat = one^(2^bit) at every set bit
    # of cnt, and square elpat once per bit. The exponent must be cnt itself
    # and the accumulate step must happen inside the loop -- an earlier
    # version started from cnt-1 and appended one final unconditional
    # convolution after the loop, which raised `one` to roughly twice the
    # intended power (nusinersen's pattern came out near 15571 Da against a
    # 7122 Da monoisotopic mass).
    elpat <- one
    remaining <- cnt
    while (remaining > 0) {
      if (remaining %% 2 == 1) pat <- .convolve_patterns(pat, elpat, threshold)
      remaining <- remaining %/% 2
      if (remaining > 0) elpat <- .convolve_patterns(elpat, elpat, threshold)
    }
  }
  pat <- pat[order(-pat$abundance), ]
  out <- head(pat, n_top)
  .builtin_iso_cache[[cache_key]] <- out
  out
}

# Unified isotope pattern dispatcher.
# threshold default was 1e-6 (0.0001% relative abundance) -- for large,
# S/P-heavy oligo formulas this forces enviPat to enumerate an enormous
# isotopologue combination space, taking many seconds per call, even though
# only the top n_top (~15-20) peaks by abundance are ever kept afterward.
# 1e-4 (0.01% relative abundance) is well below any real MS noise floor and
# essentially never changes which peaks survive the n_top truncation below,
# but cuts enviPat's enumeration cost dramatically -- this is what was
# making full-library workbook/report builds (which call this once per
# metabolite per PS-oxidation level) take tens of minutes.
#' Compute a formula's isotope pattern
#'
#' Uses enviPat if available and `use_envipat = TRUE`, else falls back to
#' a built-in binary-exponentiation convolution. Both paths are
#' memoized by (formula, threshold, n_top), since the pattern depends
#' only on the formula, never on charge state or adduct -- without that,
#' the same expensive convolution would otherwise be recomputed once per
#' charge state per adduct (see [isotope_mz_cluster()]/[compute_envelope()]).
#'
#' @param formula A formula string or vector (see [parse_formula()]).
#' @param threshold Minimum relative abundance to keep during
#'   enumeration/convolution (not the final per-peak cutoff -- only the
#'   top `n_top` peaks by abundance are kept regardless). The default
#'   (`1e-4`, 0.01% relative abundance) is well below any real MS noise
#'   floor and rarely changes which peaks survive truncation, but keeps
#'   enviPat's enumeration cost from exploding on large, S/P-heavy oligo
#'   formulas.
#' @param n_top Number of isotope peaks to keep, ranked by abundance.
#' @param use_envipat Whether to try enviPat first (`FALSE` uses the
#'   built-in convolution only).
#' @return A data.frame(`mass`, `abundance`), sorted by mass, at most
#'   `n_top` rows. `NULL` if neither path could compute a pattern (e.g.
#'   enviPat unavailable and the built-in path also failed).
#' @seealso [isotope_mz_cluster()], [compute_envelope()]
#' @export
isotope_pattern <- function(formula, threshold = 1e-4, n_top = 15,
                            use_envipat = TRUE) {
  if (is.character(formula)) formula <- parse_formula(formula)
  if (use_envipat) {
    res <- .isotope_pattern_envipat(format_formula(formula), threshold, n_top)
    if (!is.null(res)) return(res)
  }
  .isotope_pattern_builtin(formula, threshold, n_top)
}

## ---- Isotope m/z cluster for a charge state --------------------------------
#' Compute the observable isotope m/z cluster at one charge state
#'
#' Selects the monoisotopic peak (always kept, even when it is not the
#' most abundant) plus the top `n_top - 1` remaining peaks by abundance,
#' then converts each to m/z at charge `z`.
#'
#' @param formula A formula string or vector (see [parse_formula()]).
#' @param z Charge state.
#' @param n_top Number of isotope peaks to include (monoisotopic +
#'   `n_top - 1` most abundant others).
#' @param h_offset Charge-envelope offset -- see [charge_envelope()].
#' @param use_envipat Whether to use enviPat for the underlying pattern
#'   -- see [isotope_pattern()].
#' @return A data.frame(`iso`, `mz`, `abundance`, `mass`), ordered by
#'   mass (`iso = 0` is monoisotopic, increasing thereafter). `NULL` if
#'   no isotope pattern could be computed.
#' @seealso [isotope_pattern()], [match_ms1()] (its `iso_fit` score)
#' @export
isotope_mz_cluster <- function(formula, z, n_top = 8, h_offset = 0,
                               use_envipat = TRUE) {
  if (is.character(formula)) formula <- parse_formula(formula)
  # fetch a generous pattern, then select monoisotopic + most abundant
  pat <- isotope_pattern(formula, n_top = max(n_top * 3, 20),
                         use_envipat = use_envipat)
  if (is.null(pat) || nrow(pat) == 0) return(NULL)
  pat <- pat[order(pat$mass), ]
  mono <- pat[1, , drop = FALSE]              # monoisotopic = lowest mass
  rest <- pat[-1, , drop = FALSE]
  rest <- head(rest[order(-rest$abundance), , drop = FALSE], n_top - 1)
  sel <- rbind(mono, rest)
  sel <- sel[order(sel$mass), ]
  sel$mz <- (sel$mass + h_offset - z * .PROTON) / z
  sel$iso <- seq_len(nrow(sel)) - 1
  sel[, c("iso", "mz", "abundance", "mass")]
}

## ---- Adduct variants -------------------------------------------------------
#' Mass shift for a single cation adduct
#'
#' Negative-mode adduct where one cation replaces one backbone hydrogen:
#' `+Na` +21.981944 Da, `+K` +37.955882 Da, `+NH4` +17.026549 Da.
#'
#' @param adduct One of `"Na"`, `"K"`, `"NH4"`.
#' @return The mass shift, in Da.
#' @seealso [match_ms1()]
#' @export
adduct_shift <- function(adduct = c("Na", "K", "NH4")) {
  adduct <- match.arg(adduct)
  .ADDUCT_SHIFT[[adduct]]
}

## ---- Depurination variants -------------------------------------------------
#' Compute depurination (abasic-site) mass variants
#'
#' For each purine position (A or G), computes the mass/formula after
#' loss of the free base with regain of water (an abasic site).
#'
#' @param met A metabolite object (see [generate_metabolites()]).
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return A data.frame, one row per purine position: `position`, `base`,
#'   `mono_mass`, `formula_str`, `mass_shift`. `NULL` if `met` has no
#'   purine (A/G) positions.
#' @seealso [oxidation_variant()]
#' @export
depurination_variants <- function(met, dict = STANDARD_DICT) {
  purines <- which(met$bases %in% c("A", "G"))
  if (length(purines) == 0) return(NULL)
  info <- metabolite_mass_info(met, dict)
  do.call(rbind, lapply(purines, function(i) {
    base_f <- dict[[met$bases[i]]]$formula
    fv <- info$formula_vec - base_f + .as_formula(.H2O)
    data.frame(
      position = i, base = met$bases[i],
      mono_mass = formula_mass(fv, mono = TRUE),
      formula_str = format_formula(fv),
      mass_shift = -formula_mass(base_f, mono = TRUE) + .H2O_MASS,
      stringsAsFactors = FALSE)
  }))
}

## ---- +O oxidation variant --------------------------------------------------
#' Compute the +O oxidation mass variant
#'
#' One additional oxygen (e.g. on a base or sugar).
#'
#' @param met A metabolite object (see [generate_metabolites()]).
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return A list: `mono_mass`, `formula_str`, `mass_shift` (the mass of
#'   one oxygen atom).
#' @seealso [depurination_variants()]
#' @export
oxidation_variant <- function(met, dict = STANDARD_DICT) {
  info <- metabolite_mass_info(met, dict)
  fv <- info$formula_vec
  fv[["O"]] <- fv[["O"]] + 1
  list(mono_mass = formula_mass(fv, mono = TRUE),
       formula_str = format_formula(fv),
       mass_shift = .atomic_mass_mono[["O"]])
}

## ---- Full envelope table for a metabolite ----------------------------------
#' Compute the full charge/isotope envelope table for a metabolite
#'
#' One row per (PS-oxidation level, charge state, isotope peak) --
#' the table behind the app's Charge Envelopes workbook sheet and report
#' section.
#'
#' @param met A metabolite object (see [generate_metabolites()]).
#' @param z_range Charge states to include.
#' @param n_iso Isotope peaks per (oxidation level, charge state) -- see
#'   [isotope_mz_cluster()].
#' @param max_oxid Maximum PS -> PO oxidation events to model (capped at
#'   `met`'s own phosphorothioate count).
#' @param h_offset Charge-envelope offset -- see [charge_envelope()].
#' @param use_envipat Whether to use enviPat for isotope patterns -- see
#'   [isotope_pattern()].
#' @param dict A chemistry dictionary (see [build_dictionary()]).
#' @return A data.frame: `met_id`, `met_name`, `k_oxid`, `z`, `iso`,
#'   `mz`, `abundance`, `mono_mass`, `formula` -- one row per (oxidation
#'   level, charge state, isotope peak).
#' @seealso [isotope_mz_cluster()], [ps_oxidation_series()]
#' @export
compute_envelope <- function(met, z_range = 3:12, n_iso = 8,
                             max_oxid = 0, h_offset = 0,
                             use_envipat = TRUE, dict = STANDARD_DICT) {
  info <- metabolite_mass_info(met, dict)
  kmax <- min(met$n_ps, max_oxid)
  rows <- list()
  for (k in 0:kmax) {
    fv <- ps_oxid_formula(info$formula_vec, k)
    # The isotope pattern depends only on the formula, never on the charge
    # state -- compute the peak selection once per oxidation level and then
    # derive m/z for every z from it (same selection rule as
    # isotope_mz_cluster: monoisotopic peak always kept + top-(n_iso-1) by
    # abundance, ordered by mass).
    pat <- isotope_pattern(fv, n_top = max(n_iso * 3, 20),
                           use_envipat = use_envipat)
    if (is.null(pat) || nrow(pat) == 0) next
    pat <- pat[order(pat$mass), ]
    mono <- pat[1, , drop = FALSE]
    rest <- pat[-1, , drop = FALSE]
    rest <- head(rest[order(-rest$abundance), , drop = FALSE], n_iso - 1)
    sel <- rbind(mono, rest)
    sel <- sel[order(sel$mass), ]
    n_pk <- nrow(sel)
    k_mono <- formula_mass(fv, mono = TRUE)
    k_formula <- format_formula(fv)
    for (z in z_range) {
      rows[[length(rows) + 1]] <- data.frame(
        met_id = met$id, met_name = met$name, k_oxid = k,
        z = z, iso = seq_len(n_pk) - 1,
        mz = (sel$mass + h_offset - z * .PROTON) / z,
        abundance = sel$abundance,
        mono_mass = k_mono,
        formula = k_formula,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
