# =============================================================================
# default_params.R
# Single source of truth for the pipeline's default parameter values.
#
# Both the Shiny dashboard (inst/app/app.R, sidebar input `value=`s) and the
# command-line driver (run_custom_oligo.R, its PARAMS list) read their
# defaults from DEFAULT_PIPELINE_PARAMS rather than hardcoding numbers
# separately, so the two interfaces cannot drift apart the way they
# previously did (e.g. one exposing method_length/hcd_nce/target caps that
# the other silently ignored).
#
# Changing a default here changes it everywhere. A caller that needs a
# different value for one run should override it locally (e.g.
# `modifyList(DEFAULT_PIPELINE_PARAMS, list(max_oxid = 3))`), not edit this
# file.
# =============================================================================

DEFAULT_PIPELINE_PARAMS <- list(
  # -- Metabolite generation --
  max_3p         = 10,    # max 3' exonuclease truncations
  max_5p         = 10,    # max 5' exonuclease truncations
  endo           = TRUE,  # include endonuclease internal fragments
  endo_sites     = "auto",# "auto" (gap for gapmers, else all), "all", or "gap"
  min_frag_len   = 5,     # minimum internal-fragment length (nt); short endo
                          # fragments are usually minor/analytically irrelevant
  conj_catab     = TRUE,  # GalNAc conjugate catabolism: stepwise sugar loss +
                          # linker cleavage (only fires when a GalNAc conjugate
                          # is present; dominates early DMPK profiles)
  incl_scope     = "all", # inclusion-list scope: "all" or "tier1" (dominant
                          # species: parent + 3'-exo ladder + short 5'-exo)
  incl_max       = 0,     # inclusion-list cap (0 = no cap); keeps the top rows
                          # by tier then charge-state intensity

  # -- Mass & isotope --
  z_min          = 4,     # charge envelope: min charge state (auto-scaled from
                          # the parent mass at library build; see suggest_z_range)
  z_max          = 12,    # charge envelope: max charge state (auto-scaled)
  n_iso          = 7,     # isotope peaks per charge state (HRMS confirmation)
  max_oxid       = 6,     # max PS -> PO oxidation events modeled
  h_offset       = 0,     # 0 = standard [M-zH]^z-; nonzero for legacy conventions
  use_envipat    = TRUE,  # enviPat for isotope patterns (FALSE = built-in convolution)

  # -- Orbitrap Exploris acquisition method --
  method_length    = 30,  # method length / RT window end (min)
  ms2_z_min        = 4,   # PRM target list: min charge state
  ms2_z_max        = 7,   # PRM target list: max charge state
  hcd_nce          = 20,  # HCD normalized collision energy (%) -- starting value
  ms1_target_cap   = 500, # max rows in the MS1 inclusion list
  ms2_target_cap   = 300, # max rows in the MS2 PRM target list

  # -- MS matching (optional) --
  enable_ms      = FALSE,
  ppm_tol        = 10,                       # MS1 matching tolerance (ppm)
  adducts        = c("H", "Na", "K", "NH4"),
  mods           = character(0),             # optional in-place modification
                                             # classes added to the metabolite
                                             # search: "depurination",
                                             # "oxidation", "deamination",
                                             # "cyclic_phosphate" (off by default)
  frag_tol_ppm   = 25,    # fragment matching tolerance (ppm)
  frag_z_max     = 2,     # max fragment charge state considered -- z=1 dominates
                           # oligo MS2, z=2 appears for larger fragments, z=3+ is rare
  min_intensity  = 1e4,   # fixed-value noise mode: hard background/noise floor
                           # (absolute intensity) below which a peak/ROI is never
                           # considered a real feature -- was already the batch/
                           # Python ROI pipeline's own default (roi.py) but wasn't
                           # exposed to either UI path, and the single-file R
                           # path's extract_ms1_features() default of 100 is far
                           # too permissive for real HRMS data (effectively no
                           # filtering), which is exactly what floods matching
                           # with noise
  sn_threshold   = 3,     # S/N noise mode (the default mode): absolute threshold
                           # = (this file's own noise level) * sn_threshold --
                           # 3x is a conventional LOD-adjacent starting point,
                           # not a validated cutoff for any specific instrument

  # -- Batch MS processing (optional) --
  enable_batch      = FALSE,
  batch_run_ms2     = TRUE,
  batch_n_workers   = 4,
  batch_deconv_ppm  = 20
)
