# =============================================================================
# export_acquisition.R
# Export targeted mass lists for Thermo Orbitrap Exploris acquisition methods.
#
# Column layout matches the Method Editor's "Targeted Mass filter" mass list
# table exactly (Mass List Type = "m/z & z", Time Mode = "Start/End Time"),
# per the Orbitrap Exploris 120 Software Manual, "Targeted Inclusion --
# Targeted Mass filter" topic:
#   Compound, Formula, Adduct, m/z, z, Intensity Threshold,
#   t start (min), t stop (min), HCD Collision Energies (%),
#   Maximum Injection Time (ms)
# Formula/Adduct are left blank throughout -- we supply m/z and z directly,
# which is the other half of that same "m/z & z" input mode and doesn't
# require an adduct string. In the Method Editor, set Time Mode to
# "Start/End Time" (or ignore the t start/t stop columns entirely and pick
# "Unscheduled") when importing.
#
# IMPORTANT: PRM/tMS2 acquisition targets a *precursor* m/z + charge and
# records the full fragment spectrum -- individual MS2 fragment ions are
# not acquisition inputs the instrument accepts. So "MS2 target list" below
# means the precursor list used to trigger targeted MS2/PRM scans (with an
# HCD collision energy attached), not a list of fragment ions. A separate,
# non-Thermo-format fragment reference table is provided for interpreting
# the resulting spectra after acquisition (e.g. building a Skyline
# transition list or manually annotating peaks) -- see
# ms2_fragment_reference() below; do not import that one into the Method
# Editor.
# =============================================================================

## ---- Shared row builder -----------------------------------------------------
.thermo_mass_list_rows <- function(mets, dict, z_range, h_offset, max_oxid,
                                    rt_start, rt_end, mz_min, mz_max,
                                    hcd_nce = NULL, max_inj_ms = NULL) {
  rows <- list()
  for (met in mets) {
    info <- metabolite_mass_info(met, dict)
    kmax <- min(met$n_ps, max_oxid)
    for (k in 0:kmax) {
      mass <- ps_oxid_mass(info$mono_mass, k)
      compound <- if (k == 0) met$name else paste0(met$name, " +", k, "Ox")
      for (z in z_range) {
        mz <- (mass + h_offset - z * .PROTON) / z
        if (mz < mz_min || mz > mz_max) next
        rows[[length(rows) + 1]] <- data.frame(
          Compound = compound,
          Formula = "",
          Adduct = "",
          `m/z` = round(mz, 4),
          z = z,
          `Intensity Threshold` = "",
          `t start (min)` = rt_start,
          `t stop (min)` = rt_end,
          `HCD Collision Energies (%)` = if (is.null(hcd_nce)) "" else hcd_nce,
          `Maximum Injection Time (ms)` = if (is.null(max_inj_ms)) "" else max_inj_ms,
          check.names = FALSE, stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) {
    return(data.frame(Compound = character(0), Formula = character(0),
                      Adduct = character(0), `m/z` = numeric(0), z = integer(0),
                      `Intensity Threshold` = character(0),
                      `t start (min)` = numeric(0), `t stop (min)` = numeric(0),
                      `HCD Collision Energies (%)` = character(0),
                      `Maximum Injection Time (ms)` = character(0),
                      check.names = FALSE, stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

## ---- MS1 targeted inclusion list -------------------------------------------
# One row per (metabolite, PS-oxidation level, charge state) across the
# charge envelope -- for the Full Scan experiment's Targeted Mass
# (inclusion) filter, to prioritize/whitelist these masses during
# acquisition. HCD Collision Energies left blank (not applicable to a
# Full Scan).
thermo_ms1_inclusion_list <- function(mets, dict = STANDARD_DICT, z_range = 3:12,
                                      h_offset = 0, max_oxid = 0,
                                      rt_start = 0, rt_end = 30,
                                      mz_min = 100, mz_max = 6000,
                                      max_targets = 500) {
  df <- .thermo_mass_list_rows(mets, dict, z_range, h_offset, max_oxid,
                                rt_start, rt_end, mz_min, mz_max)
  if (nrow(df) > max_targets) {
    warning("MS1 inclusion list has ", nrow(df), " rows; truncating to the ",
            "highest-priority ", max_targets, " (parent + low truncations ",
            "first). Narrow z_range/max_oxid, or raise max_targets, to keep more.")
    df <- df[seq_len(max_targets), ]
  }
  df
}

## ---- MS2 PRM / targeted-MS2 precursor target list --------------------------
# Same precursor list, but tagged with an HCD collision energy so it can be
# imported as the target list for a Targeted MS2 (tMS2) / PRM scan. Default
# NCE of 20% is a reasonable starting point for HCD of PS-backbone
# oligonucleotides (see McLuckey/FMVS/OligoDistiller references in
# fragments.R) but should be optimized per instrument/method -- this is not
# a validated instrument parameter, just a starting value.
thermo_ms2_prm_target_list <- function(mets, dict = STANDARD_DICT, z_range = 3:8,
                                       h_offset = 0, max_oxid = 0,
                                       rt_start = 0, rt_end = 30,
                                       mz_min = 100, mz_max = 6000,
                                       nce = 20, max_injection_ms = 118,
                                       max_targets = 300) {
  df <- .thermo_mass_list_rows(mets, dict, z_range, h_offset, max_oxid,
                                rt_start, rt_end, mz_min, mz_max,
                                hcd_nce = nce, max_inj_ms = max_injection_ms)
  if (nrow(df) > max_targets) {
    warning("MS2 PRM target list has ", nrow(df), " rows; truncating to ",
            max_targets, ". PRM duty cycle degrades quickly with target ",
            "count -- narrow z_range/max_oxid/metabolite selection, or ",
            "raise max_targets and accept a longer cycle time.")
    df <- df[seq_len(max_targets), ]
  }
  df
}

## ---- MS2 fragment ion reference (NOT a Method Editor import) ---------------
# Theoretical McLuckey fragment ions (a/a-B/b/b-B/c/w/x/y, plus internal
# fragments) per metabolite, for interpreting acquired PRM/tMS2 spectra
# after the run -- manual annotation, or building a Skyline-style
# precursor->product transition list. This is reference data, not an
# acquisition input; the instrument records the full MS2 spectrum for each
# targeted precursor regardless of which fragments are listed here.
ms2_fragment_reference <- function(mets, dict = STANDARD_DICT, z_range = 1:2,
                                   include_internal = FALSE,
                                   ion_types = c("a", "aB", "b", "bB", "w", "y")) {
  rows <- list()
  for (met in mets) {
    if (met$n < 3) next
    info <- metabolite_mass_info(met, dict)
    frags <- generate_fragments(met, dict, ion_types = ion_types, z_range = z_range)
    if (include_internal) {
      frags <- c(frags, generate_internal_fragments(met, dict, z_range = z_range))
    }
    for (f in frags) {
      for (i in seq_len(nrow(f$mz_table))) {
        rows[[length(rows) + 1]] <- data.frame(
          met_id = f$met_id, met_name = met$name,
          precursor_mono_mass = round(info$mono_mass, 4),
          ion_type = f$ion_type, direction = f$direction,
          cleavage_site = f$cleavage_site, frag_length = f$frag_length,
          formula = f$formula, fragment_mono_mass = round(f$mono_mass, 4),
          z = f$mz_table$z[i], fragment_mz = round(f$mz_table$mz[i], 4),
          base_loss = ifelse(is.na(f$base_loss), "", f$base_loss),
          check.names = FALSE, stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}
