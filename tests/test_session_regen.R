# test_session_regen.R -- validates the assumption the two-phase workflow's
# "session resume" design rests on: build_dictionary() + parse_input() +
# generate_metabolites() + prm_inclusion_list() are pure, deterministic
# functions of the saved session parameters (see .session_input_ids and the
# Phase 1 handler in inst/app/app.R). A saved session only stores those
# parameters (not the computed spec/mets/prm objects themselves) -- "reimport
# the previous session" means restore the parameters and click "1. Generate
# Library" again, which only works if doing so reproduces byte-identical
# results. This test builds the library twice from the same fixed parameter
# set (simulating "session saved, then session restored and Phase 1
# re-run") and asserts the two runs are identical.

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (.f in c("about.R", "chemistry_dict.R", "oligo_io.R", "metabolites.R",
             "mass_isotope.R", "fragments.R")) {
  source(file.path(.pkg_root, "R", .f))
}

cat("==== Session-regeneration determinism validation ====\n\n")

# A fixed "session" -- stand-ins for the saved input$... values a real
# session JSON would restore (see .session_input_ids in inst/app/app.R).
session_params <- list(
  seq_str = INOTERSEN_TRIPLET,
  oligo_name = "inotersen",
  max_3p = 4, max_5p = 3, endo = TRUE, endo_sites = "all", min_frag_len = 6,
  z_min = 3, z_max = 8, max_oxid = 2, h_offset = 0
)

# Build the full "Phase 1" library once from the fixed params.
.build_library_from_session <- function(p) {
  dict <- build_dictionary()
  spec <- parse_input(p$seq_str, dict = dict)
  met_opts <- list(
    oligo_name = p$oligo_name, max_3p = p$max_3p, max_5p = p$max_5p,
    endo = p$endo, endo_sites = p$endo_sites, min_frag_len = p$min_frag_len
  )
  mets <- generate_metabolites(spec, opts = met_opts)
  prm <- prm_inclusion_list(mets, dict, z_range = p$z_min:p$z_max,
                             max_oxid = p$max_oxid, h_offset = p$h_offset)
  list(dict = dict, spec = spec, mets = mets, prm = prm)
}

cat("--- Run 1 (\"session saved\") ---\n")
run1 <- .build_library_from_session(session_params)
cat("Metabolites:", length(run1$mets), " PRM entries:", nrow(run1$prm), "\n")

cat("\n--- Run 2 (\"session restored, Phase 1 re-run\") ---\n")
run2 <- .build_library_from_session(session_params)
cat("Metabolites:", length(run2$mets), " PRM entries:", nrow(run2$prm), "\n")

cat("\n--- Identity checks ---\n")
stopifnot(identical(run1$spec, run2$spec))
cat("spec: identical -- PASS\n")
stopifnot(identical(run1$mets, run2$mets))
cat("mets (", length(run1$mets), "metabolites ): identical -- PASS\n")
stopifnot(identical(run1$prm, run2$prm))
cat("prm (", nrow(run1$prm), "rows ): identical -- PASS\n")

# Also check the actual metabolite/PRM TABLES (what a user would see/export)
# are identical, not just the underlying list structure -- a stricter,
# more user-facing form of the same assertion.
stopifnot(identical(metabolite_table(run1$mets), metabolite_table(run2$mets)))
cat("metabolite_table(): identical -- PASS\n")

cat("\n==== Session-regeneration is deterministic: ALL CHECKS PASSED ====\n")
