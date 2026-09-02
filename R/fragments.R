# =============================================================================
# fragments.R
# McLuckey MS/MS fragment ion generation and confirmation for oligonucleotide
# metabolite identification.
#
# Implements the McLuckey nomenclature for oligonucleotide tandem MS:
#   5'-terminal ions: a, b, c, d  (cleavage at C3'-O3', O3'-P, P-O5', C4'-C3')
#   3'-terminal ions: w, x, y, z  (complementary cleavage sites)
#   Base-loss ions:   a-B, b-B    (loss of 3'-terminal base + H2O)
#   Internal ions:    w-a, w-b    (double cleavage, w-type 5' + a/b-type 3')
#
# Fragment m/z computed in negative ESI mode: [ion - zH]^z-
#   m/z = (fragment_mass - z * proton) / z
#
# Grounded in:
#   - McLuckey et al., JASMS 1992: original oligonucleotide fragment nomenclature
#   - Eluforsen (Kim et al., 2019): PS diagnostic ion at m/z 94.9452
#   - FMVS (Ye et al., 2025): confirmation score + sequence coverage >50%
#     threshold; y/b predominant for MOE-PS; a-B/w for DNA-PS; internal w-a/w-b/w-d
#   - OligoDistiller (Liu et al., 2025): w-a internal fragments, sequence coverage
#     = cleavage sites observed / total possible
# =============================================================================

## ---- Fragment ion mass formulas -------------------------------------------
# For a 5'-terminal fragment of length k (positions 1..k):
#   F_k = intact k-mer formula (3'-OH terminus, conj5, no conj3)
#   L_k = linkage formula at cleavage site k (bond between pos k and k+1)
#
#   a_k   = F_k - O          (C3'-O3' cleavage, 5' retains C3'-H)
#   a-B_k = a_k - base_k + H2O  (additional base loss from 3'-terminal pos)
#   b_k   = F_k              (O3'-P cleavage, retains 3'-OH)
#   b-B_k = b_k - base_k + H2O
#   c_k   = F_k + L_k - H    (P-O5' cleavage, retains P-OH)
#   d_k   = a_k - base_k + H2O  (C4'-C3' cleavage, approx as a with base loss)
#
# For a 3'-terminal fragment of length n-k (positions k+1..n):
#   G_{nk} = intact (n-k)-mer formula (5'-OH terminus, conj3, no conj5)
#
#   w_{nk}   = G_{nk} + L_k + O     (C3'-O3' cleavage, 3' retains O3'-H + linkage)
#   x_{nk}   = G_{nk} + L_k - H     (O3'-P cleavage, retains P-OH)
#   y_{nk}   = G_{nk}               (P-O5' cleavage, retains 5'-OH)
#   z_{nk}   = w_{nk} + base_{k+1} - H2O  (C4'-C3', approx: w + 5'-terminal base)
#
# Complementary check: a_k + w_{nk} = F_k + G_{nk} + L_k = M + H2O
# (F_k + G_{nk} + L_k = M + H2O by the condensation accounting)

## ---- PS diagnostic ions ---------------------------------------------------
# Characteristic fragment ions for phosphorothioate backbones.
# m/z 94.9452 and 192.9746 are PS-specific markers (from eluforsen/FMVS papers).
.PS_DIAGNOSTIC_MZ <- c(94.9452, 192.9746)

## ---- Helper: build a fragment object --------------------------------------
.make_frag <- function(ion_type, direction, cleavage_site, frag_length,
                       formula_vec, dict, z_range, h_offset, met_id,
                       base_loss = NA, internal_5 = NA, internal_3 = NA) {
  mass <- formula_mass(formula_vec, mono = TRUE)
  mz_rows <- lapply(z_range, function(z) {
    data.frame(z = z, mz = (mass - z * .PROTON) / z,
               stringsAsFactors = FALSE)
  })
  mz_df <- do.call(rbind, mz_rows)
  list(
    met_id       = met_id,
    ion_type     = ion_type,
    direction    = direction,
    cleavage_site = cleavage_site,
    frag_length  = frag_length,
    formula      = format_formula(formula_vec),
    formula_vec  = formula_vec,
    mono_mass    = mass,
    mz_table     = mz_df,
    base_loss    = base_loss,
    internal_5   = internal_5,
    internal_3   = internal_3
  )
}

## ---- Generate terminal fragment ions --------------------------------------
# Generate all McLuckey terminal ions for a metabolite.
# opts:
#   ion_types  : character vector of ion types to generate
#   z_range    : charge states for m/z calculation (default 1:2)
#   h_offset   : mass offset for non-standard envelope conventions
#   include_dz : include d/z ions (approximate, default FALSE)
#
# Default ion_types omits c and x: in negative-ion CID of oligonucleotides
# w and (a-B) ions dominate observed spectra by a wide margin, c/y appear
# but weaker, and b/x -- x especially -- are rarely the ions actually seen.
# Pass ion_types = c("a","aB","b","bB","c","w","x","y") for full McLuckey
# coverage (e.g. HCD, or other fragmentation methods that behave differently).
generate_fragments <- function(met, dict = STANDARD_DICT,
                                ion_types = c("a", "aB", "b", "bB", "w", "y"),
                                z_range = 1:2, h_offset = 0,
                                include_dz = FALSE) {
  if (include_dz) ion_types <- union(ion_types, c("d", "z"))
  n <- met$n
  if (n < 3) return(list())  # too short for meaningful fragments
  frags <- list()

  for (k in seq_len(n - 1)) {
    # 5'-terminal fragment: positions 1..k
    F_k <- assemble_oligo_formula(met$bases[1:k], met$sugars[1:k],
                                   met$linkages[1:k], met$conj5, "none", dict)
    lk_code <- met$linkages[k]
    L_k <- if (!is.na(lk_code) && nzchar(lk_code)) dict[[lk_code]]$formula else .empty_formula()
    base_k <- dict[[met$bases[k]]]$formula

    # a-ion: F_k - O
    if ("a" %in% ion_types) {
      fa <- F_k; fa[["O"]] <- fa[["O"]] - 1
      frags <- c(frags, list(.make_frag("a", "5'", k, k, fa, dict, z_range,
                                         h_offset, met$id)))
    }
    # a-B ion: a - base_k + H2O
    if ("aB" %in% ion_types) {
      faB <- F_k; faB[["O"]] <- faB[["O"]] - 1
      faB <- add_formulas(faB, -base_k)
      faB <- add_formulas(faB, .as_formula(.H2O))
      frags <- c(frags, list(.make_frag("a-B", "5'", k, k, faB, dict, z_range,
                                         h_offset, met$id, base_loss = met$bases[k])))
    }
    # b-ion: F_k
    if ("b" %in% ion_types) {
      frags <- c(frags, list(.make_frag("b", "5'", k, k, F_k, dict, z_range,
                                         h_offset, met$id)))
    }
    # b-B ion: b - base_k + H2O
    if ("bB" %in% ion_types) {
      fbB <- add_formulas(add_formulas(F_k, -base_k), .as_formula(.H2O))
      frags <- c(frags, list(.make_frag("b-B", "5'", k, k, fbB, dict, z_range,
                                         h_offset, met$id, base_loss = met$bases[k])))
    }
    # c-ion: F_k + L_k - H
    if ("c" %in% ion_types) {
      fc <- add_formulas(F_k, L_k); fc[["H"]] <- fc[["H"]] - 1
      frags <- c(frags, list(.make_frag("c", "5'", k, k, fc, dict, z_range,
                                         h_offset, met$id)))
    }
    # d-ion (approximate): a - base_k + H2O
    if ("d" %in% ion_types) {
      fd <- F_k; fd[["O"]] <- fd[["O"]] - 1
      fd <- add_formulas(add_formulas(fd, -base_k), .as_formula(.H2O))
      frags <- c(frags, list(.make_frag("d", "5'", k, k, fd, dict, z_range,
                                         h_offset, met$id, base_loss = met$bases[k])))
    }

    # 3'-terminal fragment: positions (k+1)..n
    idx <- (k + 1):n
    G_nk <- assemble_oligo_formula(met$bases[idx], met$sugars[idx],
                                    met$linkages[idx], "none", met$conj3, dict)
    base_kp1 <- dict[[met$bases[k + 1]]]$formula

    # w-ion: G_nk + L_k + O
    if ("w" %in% ion_types) {
      fw <- add_formulas(G_nk, L_k); fw[["O"]] <- fw[["O"]] + 1
      frags <- c(frags, list(.make_frag("w", "3'", k, n - k, fw, dict, z_range,
                                         h_offset, met$id)))
    }
    # x-ion: G_nk + L_k - H
    if ("x" %in% ion_types) {
      fx <- add_formulas(G_nk, L_k); fx[["H"]] <- fx[["H"]] - 1
      frags <- c(frags, list(.make_frag("x", "3'", k, n - k, fx, dict, z_range,
                                         h_offset, met$id)))
    }
    # y-ion: G_nk
    if ("y" %in% ion_types) {
      frags <- c(frags, list(.make_frag("y", "3'", k, n - k, G_nk, dict, z_range,
                                         h_offset, met$id)))
    }
    # z-ion (approximate): w + base_{k+1} - H2O
    if ("z" %in% ion_types) {
      fz <- add_formulas(G_nk, L_k); fz[["O"]] <- fz[["O"]] + 1
      fz <- add_formulas(add_formulas(fz, base_kp1), -(.as_formula(.H2O)))
      frags <- c(frags, list(.make_frag("z", "3'", k, n - k, fz, dict, z_range,
                                         h_offset, met$id, base_loss = met$bases[k + 1])))
    }
  }
  frags
}

## ---- Generate internal fragments (double cleavage) ------------------------
# Internal fragments from two backbone cleavages.
# Types: w-a (w-type 5' end + a-type 3' end), w-b, w-d
# Fragment spans positions (i+1)..j where i < j, i = 5' cleavage, j = 3' cleavage
#
# w-a(i,j) = I_{i+1..j} + L_i - H   (w adds L_i + O - H, a removes O)
# w-b(i,j) = I_{i+1..j} + L_i + O - H  (w adds L_i + O - H, b retains 3'-OH)
# w-d(i,j) = I_{i+1..j} + L_i - H - base_j + H2O  (w + d with base loss)
generate_internal_fragments <- function(met, dict = STANDARD_DICT,
                                         types = c("wa", "wb"),
                                         z_range = 1:2, h_offset = 0,
                                         min_len = 2, max_len = 8) {
  n <- met$n
  if (n < 4) return(list())
  frags <- list()

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):min(n - 1, i + max_len - 1)) {
      if (j - i < min_len - 1) next
      idx <- (i + 1):j
      I_ij <- assemble_oligo_formula(met$bases[idx], met$sugars[idx],
                                      met$linkages[idx], "none", "none", dict)
      lk_i <- met$linkages[i]
      L_i <- if (!is.na(lk_i) && nzchar(lk_i)) dict[[lk_i]]$formula else .empty_formula()
      base_j <- dict[[met$bases[j]]]$formula

      if ("wa" %in% types) {
        fwa <- add_formulas(I_ij, L_i); fwa[["H"]] <- fwa[["H"]] - 1
        frags <- c(frags, list(.make_frag("w-a", "internal", i, j - i, fwa,
                                           dict, z_range, h_offset, met$id,
                                           internal_5 = i, internal_3 = j)))
      }
      if ("wb" %in% types) {
        fwb <- add_formulas(I_ij, L_i); fwb[["O"]] <- fwb[["O"]] + 1
        fwb[["H"]] <- fwb[["H"]] - 1
        frags <- c(frags, list(.make_frag("w-b", "internal", i, j - i, fwb,
                                           dict, z_range, h_offset, met$id,
                                           internal_5 = i, internal_3 = j)))
      }
      if ("wd" %in% types) {
        fwd <- add_formulas(I_ij, L_i); fwd[["H"]] <- fwd[["H"]] - 1
        fwd <- add_formulas(add_formulas(fwd, -base_j), .as_formula(.H2O))
        frags <- c(frags, list(.make_frag("w-d", "internal", i, j - i, fwd,
                                           dict, z_range, h_offset, met$id,
                                           base_loss = met$bases[j],
                                           internal_5 = i, internal_3 = j)))
      }
    }
  }
  frags
}

## ---- Fragment intensity heuristic ------------------------------------------
# Rule-based *relative* intensity weight for a fragment ion. This is NOT a
# calibrated intensity model -- there is no oligonucleotide equivalent of
# the large annotated spectral libraries that peptide tools (Prosit,
# MS2PIP) train on, so nothing here is fit to data. It only encodes
# fragmentation propensities that are already well established in the
# literature already cited at the top of this file:
#
#   - PS (phosphorothioate) linkages are markedly more labile than PO under
#     CID/HCD (this is also why PS_DIAGNOSTIC_MZ exists below), so a
#     cleavage site sitting on a PS linkage is weighted up.
#   - Ye et al. 2025 (FMVS) report y/b ions predominant for MOE-PS
#     chemistry and a-B/w ions predominant for DNA-PS chemistry -- the
#     sugar immediately adjacent to the cleavage site selects which ion
#     type gets the boost.
#   - Purine glycosidic bonds (A/G-type bases) are more labile than
#     pyrimidine ones under CID, so a-B/b-B/w-d base-loss ions are weighted
#     up when the lost base is a purine. Purine-ness is read off the base's
#     formula (>= 4 ring nitrogens) rather than hardcoded letter codes, so
#     it still works for custom/modified purine bases entered via the
#     Custom Chemistry table.
#   - Internal (double-cleavage) ions are consistently minor relative to
#     terminal ions in reported oligo MS2 spectra, so they carry a fixed
#     down-weight rather than the ion-type/sugar logic above.
#
# Treat the result as a coarse "expect this peak to be relatively taller"
# ranking, not a predicted abundance -- confirm every real assignment
# against acquired data (see ../DISCLAIMER.md).
.PS_LINKAGE_BOOST     <- 1.5  # PS vs PO cleavage-site lability
.CHEMISTRY_ION_BOOST  <- 1.6  # dominant ion type for the local sugar chemistry
.PURINE_LOSS_BOOST    <- 1.3  # purine vs pyrimidine base-loss lability
.INTERNAL_ION_WEIGHT  <- 0.35 # internal ions are minor vs. terminal ions
.MINOR_ION_WEIGHT     <- 0.7  # ion types with no specific literature boost

# Purines (adenine/guanine-type bases, including modified ones such as
# 2,6-diaminopurine or hypoxanthine) have a fused 5-/6-membered ring system
# with >= 4 ring nitrogens; pyrimidines (cytosine/thymine/uracil-type,
# including 5-methyl/5-hydroxymethyl variants) have at most 3. Reading this
# off the formula -- rather than matching specific base codes -- keeps the
# heuristic working for any base the user adds via the Custom Chemistry
# table.
.is_purine_base <- function(base_code, dict) {
  entry <- dict[[base_code]]
  !is.null(entry) && isTRUE(unname(entry$formula[["N"]]) >= 4)
}

# PS (phosphorothioate) linkages replace one non-bridging phosphate oxygen
# with sulfur, so they're identified by formula (any S) rather than the "s"
# code alone -- this also catches "u" (the PS stereochemistry variant) and
# any custom PS-like linkage added via the Custom Chemistry table.
.is_ps_linkage <- function(linkage_code, dict) {
  if (is.na(linkage_code) || !nzchar(linkage_code)) return(FALSE)
  entry <- dict[[linkage_code]]
  !is.null(entry) && isTRUE(unname(entry$formula[["S"]]) > 0)
}

# Base ion-type weight for one terminal fragment, before the PS-linkage and
# purine-loss adjustments below. `sugar_code` is the sugar immediately
# adjacent to the cleavage site (see fragment_intensity_weight()).
.terminal_ion_weight <- function(ion_type, sugar_code) {
  is_dna_sugar <- identical(sugar_code, "d")
  is_moe_ome_sugar <- sugar_code %in% c("e", "MOE", "m")
  switch(ion_type,
    "a-B" = if (is_dna_sugar) .CHEMISTRY_ION_BOOST else 1.0,
    "w"   = if (is_dna_sugar) .CHEMISTRY_ION_BOOST else 1.0,
    "y"   = if (is_moe_ome_sugar) .CHEMISTRY_ION_BOOST else 1.0,
    "b"   = if (is_moe_ome_sugar) .CHEMISTRY_ION_BOOST else 1.0,
    .MINOR_ION_WEIGHT)
}

# Relative intensity weight for one fragment object (as produced by
# .make_frag()), independent of fragment charge. Always positive; peaks for
# a given spectrum should be rescaled (e.g. to a 0-100 max) by the caller,
# matching the MS1 library's convention -- see build_ms2_library().
fragment_intensity_weight <- function(f, met, dict = STANDARD_DICT) {
  if (identical(f$direction, "internal")) {
    w <- .INTERNAL_ION_WEIGHT
    if (.is_ps_linkage(met$linkages[f$internal_5], dict)) w <- w * .PS_LINKAGE_BOOST
    return(w)
  }

  site <- f$cleavage_site
  # 5'-ions break the bond just after position `site`; 3'-ions break the
  # same bond from the other side -- either way, "adjacent sugar" means the
  # sugar at `site` for 5'-ions and at `site + 1` for 3'-ions.
  sugar_code <- if (identical(f$direction, "5'")) met$sugars[site] else met$sugars[site + 1]
  w <- .terminal_ion_weight(f$ion_type, sugar_code)

  if (.is_ps_linkage(met$linkages[site], dict)) w <- w * .PS_LINKAGE_BOOST

  if (!is.na(f$base_loss) && .is_purine_base(f$base_loss, dict)) {
    w <- w * .PURINE_LOSS_BOOST
  }

  w
}

## ---- Flatten fragments to a display table ---------------------------------
# Returns data.frame with one row per (fragment, charge state).
fragment_table <- function(frags) {
  if (length(frags) == 0) return(data.frame())
  rows <- lapply(frags, function(f) {
    cbind(data.frame(
      met_id = f$met_id, ion_type = f$ion_type, direction = f$direction,
      cleavage_site = f$cleavage_site, frag_length = f$frag_length,
      formula = f$formula, mono_mass = sprintf("%.4f", f$mono_mass),
      base_loss = ifelse(is.na(f$base_loss), "", f$base_loss),
      stringsAsFactors = FALSE
    ), f$mz_table)
  })
  do.call(rbind, rows)
}

## ---- Match fragments against MS2 peaks ------------------------------------
# Given theoretical fragments and observed MS2 peaks, find matches.
#
# ms2_peaks: data.frame with columns mz, intensity (and optional rt)
# tol_ppm:   mass tolerance in ppm
# z_range:   charge states to try for matching
#
# Returns data.frame with matched fragments and their scores.
match_fragments <- function(frags, ms2_peaks, tol_ppm = 25,
                             z_range = 1:2, h_offset = 0) {
  if (length(frags) == 0 || nrow(ms2_peaks) == 0) return(data.frame())

  matches <- list()
  for (f in frags) {
    for (z in z_range) {
      # theoretical m/z at this charge
      theo_mz <- (f$mono_mass - z * .PROTON) / z
      # find closest observed peak
      dmz <- abs(ms2_peaks$mz - theo_mz)
      best <- which.min(dmz)
      ppm_err <- dmz[best] / theo_mz * 1e6
      if (ppm_err <= tol_ppm) {
        matches[[length(matches) + 1]] <- data.frame(
          met_id = f$met_id, ion_type = f$ion_type,
          direction = f$direction, cleavage_site = f$cleavage_site,
          frag_length = f$frag_length, z = z,
          theo_mz = theo_mz, obs_mz = ms2_peaks$mz[best],
          ppm_error = ppm_err, intensity = ms2_peaks$intensity[best],
          formula = f$formula, base_loss = ifelse(is.na(f$base_loss), "", f$base_loss),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(matches) == 0) return(data.frame())
  do.call(rbind, matches)
}

## ---- Check for PS diagnostic ions in MS2 data ------------------------------
# Returns list of found diagnostic ions with observed m/z and intensity.
check_ps_diagnostic <- function(ms2_peaks, tol_ppm = 50) {
  if (nrow(ms2_peaks) == 0) return(data.frame())
  rows <- list()
  for (dmz in .PS_DIAGNOSTIC_MZ) {
    dmz_obs <- abs(ms2_peaks$mz - dmz)
    best <- which.min(dmz_obs)
    ppm_err <- dmz_obs[best] / dmz * 1e6
    if (ppm_err <= tol_ppm) {
      rows[[length(rows) + 1]] <- data.frame(
        diagnostic_mz = dmz, obs_mz = ms2_peaks$mz[best],
        ppm_error = ppm_err, intensity = ms2_peaks$intensity[best],
        stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

## ---- Sequence coverage ----------------------------------------------------
# Fraction of backbone cleavage sites confirmed by at least one matched ion.
# cleavage sites: 1..(n-1) (between each pair of adjacent nucleotides)
# A site is "covered" if any matched fragment has cleavage_site == that site.
#
# matched: data.frame from match_fragments()
# n:       oligonucleotide length
sequence_coverage <- function(matched, n) {
  total_sites <- n - 1
  if (total_sites < 1) return(list(coverage = 0, covered_sites = integer(0),
                                    total_sites = 0))
  if (nrow(matched) == 0) return(list(coverage = 0, covered_sites = integer(0),
                                       total_sites = total_sites))
  covered <- unique(matched$cleavage_site)
  list(coverage = length(covered) / total_sites,
       covered_sites = covered,
       total_sites = total_sites)
}

## ---- Confirmation score ---------------------------------------------------
# Composite score (0-100) for metabolite identification confidence.
# Based on FMVS paper: total confirmation score = similarity of experimental
# to theoretical isotope envelope (MS1); sequence coverage (MS2) > 50% threshold.
#
# We compute a practical score combining:
#   - sequence coverage (0-50 points)
#   - number of matched fragment ions (0-25 points)
#   - mass accuracy of matches (0-15 points)
#   - PS diagnostic ion presence (0-10 points)
#
# matched:     data.frame from match_fragments()
# n:           oligonucleotide length
# diagnostics: data.frame from check_ps_diagnostic()
confirmation_score <- function(matched, n, diagnostics = NULL) {
  cov <- sequence_coverage(matched, n)
  n_match <- if (is.null(matched) || nrow(matched) == 0) 0 else nrow(matched)
  n_diag <- if (is.null(diagnostics) || nrow(diagnostics) == 0) 0 else nrow(diagnostics)

  # Coverage score (0-50): linear up to 80% coverage, then capped
  cov_score <- min(cov$coverage / 0.8, 1) * 50

  # Match count score (0-25): 10+ matches = full score
  match_score <- min(n_match / 10, 1) * 25

  # Mass accuracy score (0-15): based on median ppm error
  if (n_match > 0) {
    med_ppm <- median(matched$ppm_error)
    acc_score <- max(0, 1 - med_ppm / 25) * 15  # 0 ppm = 15, 25 ppm = 0
  } else {
    acc_score <- 0
  }

  # Diagnostic ion score (0-10): 5 per diagnostic ion found
  diag_score <- min(n_diag * 5, 10)

  total <- cov_score + match_score + acc_score + diag_score
  list(
    total_score = round(total, 1),
    coverage = cov$coverage,
    coverage_score = round(cov_score, 1),
    n_matches = n_match,
    match_score = round(match_score, 1),
    median_ppm = if (n_match > 0) round(median(matched$ppm_error), 2) else NA,
    accuracy_score = round(acc_score, 1),
    n_diagnostics = n_diag,
    diagnostic_score = round(diag_score, 1),
    covered_sites = cov$covered_sites,
    total_sites = cov$total_sites,
    confident = total >= 75 && cov$coverage >= 0.5
  )
}

## ---- Full fragment confirmation for a metabolite ---------------------------
# Convenience: generate fragments, match, score in one call.
# Returns list with fragments, matches, diagnostics, coverage, score.
confirm_metabolite <- function(met, ms2_peaks, dict = STANDARD_DICT,
                                tol_ppm = 25, z_range = 1:2,
                                ion_types = c("a", "aB", "b", "bB", "w", "y"),
                                include_internal = FALSE,
                                include_dz = FALSE, h_offset = 0) {
  frags <- generate_fragments(met, dict, ion_types, z_range, h_offset, include_dz)
  if (include_internal) {
    frags <- c(frags, generate_internal_fragments(met, dict, z_range = z_range,
                                                    h_offset = h_offset))
  }
  matched <- match_fragments(frags, ms2_peaks, tol_ppm, z_range, h_offset)
  diags <- check_ps_diagnostic(ms2_peaks, tol_ppm = max(tol_ppm, 50))
  score <- confirmation_score(matched, met$n, diags)
  list(
    metabolite = met, fragments = frags, matched = matched,
    diagnostics = diags, coverage = score$coverage, score = score
  )
}

## ---- PRM inclusion list export --------------------------------------------
# Generate a targeted inclusion list for PRM acquisition from the metabolite
# library. Returns data.frame with precursor m/z, charge, metabolite name.
# This is a first-class output of the library generator for MS acquisition.
prm_inclusion_list <- function(mets, dict = STANDARD_DICT, z_range = 3:12,
                                h_offset = 0, max_oxid = 0,
                                ppm_window = 10) {
  rows <- list()
  for (met in mets) {
    info <- metabolite_mass_info(met, dict)
    kmax <- min(met$n_ps, max_oxid)
    for (k in 0:kmax) {
      mass <- ps_oxid_mass(info$mono_mass, k)
      for (z in z_range) {
        mz <- (mass + h_offset - z * .PROTON) / z
        if (mz < 100 || mz > 3000) next  # typical Q1 range
        rows[[length(rows) + 1]] <- data.frame(
          met_id = met$id, met_name = met$name,
          kind = met$kind, n = met$n,
          k_oxid = k, z = z,
          precursor_mz = round(mz, 4),
          isolation_ppm = ppm_window,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}
