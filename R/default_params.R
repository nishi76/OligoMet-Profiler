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
  endo_sites     = "all", # "all" positions, or "gap" (DNA gap only)
  min_frag_len   = 3,     # minimum fragment length (nt)

  # -- Mass & isotope --
  z_min          = 3,     # charge envelope: min charge state
  z_max          = 12,    # charge envelope: max charge state
  n_iso          = 5,     # isotope peaks per charge state
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
  frag_tol_ppm   = 25,    # fragment matching tolerance (ppm)
  frag_z_max     = 3,     # max fragment charge state considered

  # -- Batch MS processing (optional) --
  enable_batch      = FALSE,
  batch_run_ms2     = TRUE,
  batch_n_workers   = 4,
  batch_deconv_ppm  = 20
)
