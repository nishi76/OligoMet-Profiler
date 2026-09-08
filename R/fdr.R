# =============================================================================
# fdr.R
# Heuristic target-decoy FDR estimate for MS1 metabolite matching.
#
# Opt-in only -- this is NOT wired into annotate_metabolites()'s automatic
# summary. It runs the full matching pipeline n_decoys+1 times, which is
# real added compute an ordinary run shouldn't pay, and it's an unvalidated
# heuristic (see the caveat below) that the user should invoke deliberately.
#
#   FDR = mean(n_decoy_matches) / n_target_matches
#
# Decoy strategy: permute base identity across positions, keeping every
# other per-position attribute (sugar, linkage, conjugates) fixed. This is
# the oligo analog of peptide target-decoy sequence reversal -- same
# physical-chemistry mass contributions at each position (so the decoy
# parent mass is identical to the real one), different sequence identity
# (so sequence-dependent fragment masses genuinely differ). Permuting
# sugars/linkages instead would change the parent mass itself (PS vs PO
# linkages differ in mass, as do sugar modifications), which would trivially
# fail to even MS1-match and produce a degenerate always-zero decoy rate --
# not a meaningful decoy.
#
# Caveat (real, not hedging): unlike protein sequence reversal, a
# base-permuted oligo decoy still draws from the same small base alphabet
# as real biology, so decoy fragment masses can coincide with real
# dictionary entries at non-trivial rates. This is a rough global estimate,
# not a validated per-match FDR -- treat it as a sanity check on match
# volume, not a rigorous confidence bound.
# =============================================================================

## ---- Decoy sequence generation ---------------------------------------------
# Permute spec$bases (a random shuffle across positions); sugars, linkages,
# conj5, conj3 all stay identical to the real spec.
.make_decoy_spec <- function(spec, dict = STANDARD_DICT, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  decoy <- spec
  decoy$bases <- sample(spec$bases)
  decoy$raw <- paste0(spec$raw, " [DECOY]")
  parse_input(decoy, dict = dict, notation = "structured")
}

## ---- Target-decoy FDR estimate ---------------------------------------------
# spec: the real oligo_spec (from parse_input()) that `mets` was generated
#   from -- needed here to build decoy specs with the same structure.
# mets: the real metabolite series (generate_metabolites(spec, opts)).
# opts: the same opts used to build `mets`, reused to build each decoy's
#   metabolite series with identical truncation/endo-cleavage depth.
# All other parameters mirror annotate_metabolites()'s exactly, since the
# whole point of a target-decoy comparison is scoring both under identical
# pipeline settings.
#
# FDR is estimated at the MS1 level only (most defensible given the decoy
# caveat above) -- this does not attempt a combined MS1+MS2 FDR.
estimate_fdr <- function(spec, mets, ms1_features, ms2_data = NULL, dict = STANDARD_DICT,
                          n_decoys = 1, opts = list(),
                          ppm_tol = 10, z_range = 3:12,
                          adducts = c("H", "Na", "K", "NH4"),
                          max_oxid = 6, h_offset = 0, n_iso = 5, use_envipat = TRUE,
                          frag_tol_ppm = 25, frag_z_range = 1:2, include_internal = FALSE,
                          seed = NULL) {
  target <- annotate_metabolites(mets, ms1_features, ms2_data, dict, ppm_tol, z_range,
                                  adducts, max_oxid, h_offset, frag_tol_ppm, frag_z_range,
                                  n_iso, use_envipat, include_internal)
  n_target_matches <- nrow(target$ms1_matches)

  decoy_summaries <- list()
  n_decoy_matches_vec <- integer(n_decoys)
  for (i in seq_len(n_decoys)) {
    decoy_spec <- .make_decoy_spec(spec, dict, seed = if (!is.null(seed)) seed + i else NULL)
    decoy_mets <- generate_metabolites(decoy_spec, opts)
    decoy_result <- annotate_metabolites(decoy_mets, ms1_features, ms2_data, dict, ppm_tol,
                                          z_range, adducts, max_oxid, h_offset, frag_tol_ppm,
                                          frag_z_range, n_iso, use_envipat, include_internal)
    n_decoy_matches_vec[i] <- nrow(decoy_result$ms1_matches)
    decoy_summaries[[i]] <- decoy_result$summary
  }

  fdr <- if (n_target_matches > 0) mean(n_decoy_matches_vec) / n_target_matches else NA_real_
  list(fdr = fdr, n_target_matches = n_target_matches,
       n_decoy_matches = mean(n_decoy_matches_vec), n_decoys = n_decoys,
       target_summary = target$summary, decoy_summaries = decoy_summaries)
}
