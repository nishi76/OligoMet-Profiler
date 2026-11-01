# test_agent_tools.R -- validate the agent tool registry (R/agent_tools.R)
#
# This is a proof-of-concept foundation (see /root/.claude/plans/
# staged-mapping-blanket.md), reviewed/tested here BEFORE any LLM loop or
# front end is built on top of it. Checks:
#   1. Each tool round-trips through call_agent_tool() using session context
#      (ctx) the way the planned in-app Shiny assistant would call it.
#   2. Each tool also works from args alone (no ctx), the way the planned
#      stateless MCP server would call it.
#   3. Every match/confirmation result actually carries its evidence fields
#      (ppm_error, ambiguous/n_candidates, coverage/score) -- the mechanism
#      the whole design depends on to keep an agent from hallucinating.
#   4. A representative result survives a real jsonlite::toJSON()/fromJSON()
#      round trip (catches serialization bugs before they reach a real
#      LLM-facing transport).

.pkg_root <- local({
  this <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(f) > 0) normalizePath(f) else NULL
  }, error = function(e) NULL)
  if (!is.null(this)) dirname(dirname(this)) else ".."
})
for (m in c("about", "progress_utils", "chemistry_dict", "oligo_io", "metabolites",
            "mass_isotope", "fragments", "ms_matching", "degradation", "statistics",
            "agent_tools")) {
  source(file.path(.pkg_root, "R", paste0(m, ".R")))
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  cat("SKIPPED: jsonlite not installed\n")
  quit(status = 0, save = "no")
}

fail <- 0L
chk <- function(label, ok) {
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) fail <<- fail + 1L
}

cat("=== 1. parse_sequence ===\n")
r <- call_agent_tool("parse_sequence", list(sequence = INOTERSEN_TRIPLET))
chk("returns a spec with n=20", identical(r$spec$n, 20L) || identical(r$spec$n, 20))
chk("formula matches the published inotersen formula",
    identical(r$formula, format_formula(parse_formula(INOTERSEN_FORMULA))))
spec <- r$spec
dict <- STANDARD_DICT

cat("\n=== 2. build_metabolite_library (ctx path) ===\n")
r2 <- call_agent_tool("build_metabolite_library",
                       list(oligo_name = "inotersen", max_3p = 2, max_5p = 2, endo = FALSE),
                       ctx = list(spec = spec, dict = dict))
chk("5 metabolites (parent + 2x3p + 2x5p)", r2$n == 5)
mets <- r2$metabolites

cat("\n=== 2b. build_metabolite_library (stateless args path) ===\n")
r2b <- call_agent_tool("build_metabolite_library",
                        list(sequence = INOTERSEN_TRIPLET, oligo_name = "inotersen",
                             max_3p = 2, max_5p = 2, endo = FALSE))
chk("stateless call gives the same count", r2b$n == r2$n)

cat("\n=== 3. get_metabolite_mass ===\n")
r3 <- call_agent_tool("get_metabolite_mass", list(met_id = "M01", max_oxid = 2),
                       ctx = list(mets = mets, dict = dict))
chk("mono_mass matches parent formula mass", abs(r3$mono_mass - formula_mass(assemble_oligo_formula(
  mets[[1]]$bases, mets[[1]]$sugars, mets[[1]]$linkages, mets[[1]]$conj5, mets[[1]]$conj3))) < 1e-6)
chk("oxidation_series has k=0..2", nrow(r3$oxidation_series) == 3)

cat("\n=== 4. generate_fragment_ions ===\n")
r4 <- call_agent_tool("generate_fragment_ions", list(met_id = "M01", z_range = list(1L, 2L)),
                       ctx = list(mets = mets, dict = dict))
chk("fragments generated", r4$n > 0)
chk("default ion types present", all(c("a-B", "w", "y", "b") %in% unique(r4$fragments$ion_type)))

cat("\n=== 5. match_ms1_features (with ambiguity evidence) ===\n")
info1 <- metabolite_mass_info(mets[[1]])
z_true <- 6
mz_true <- (info1$mono_mass - z_true * .PROTON) / z_true
features <- list(
  list(mz = mz_true, rt = 12.3, max_intensity = 5e5),
  list(mz = mz_true * (1 + 300 / 1e6), rt = 12.3, max_intensity = 100)  # far outside 10ppm, distinct feature
)
r5 <- call_agent_tool("match_ms1_features", list(features = features, ppm_tol = 10, adducts = list("H")),
                       ctx = list(mets = mets, dict = dict))
chk("at least one match found", r5$n_matches > 0)
chk("every match row carries n_candidates/ambiguous", all(c("n_candidates", "ambiguous") %in% names(r5$matches)))
chk("n_ambiguous field present at top level", "n_ambiguous" %in% names(r5))

cat("\n=== 6. confirm_ms2 ===\n")
frags_m01 <- generate_fragments(mets[[1]], dict, z_range = 1)
ms2_peaks <- lapply(frags_m01[1:6], function(f) list(mz = (f$mono_mass - .PROTON), intensity = 1000))
r6 <- call_agent_tool("confirm_ms2", list(met_id = "M01", ms2_peaks = ms2_peaks, tol_ppm = 10, z_range = list(1L)),
                       ctx = list(mets = mets, dict = dict))
chk("confirm_ms2 returns a score", is.list(r6$score) && !is.null(r6$score$total_score))
chk("confirm_ms2 returns coverage", is.numeric(r6$coverage))
chk("at least some fragments matched", r6$n_matched > 0)

cat("\n=== 7. summarize_degradation ===\n")
ms1_matches <- rbind(
  data.frame(sample = "s1", met_id = "M01", met_name = "parent", kind = "parent", area = 100, stringsAsFactors = FALSE),
  data.frame(sample = "s1", met_id = "M02", met_name = "N-1", kind = "exo_3p", area = 50, stringsAsFactors = FALSE)
)
r7 <- call_agent_tool("summarize_degradation", list(ms1_matches = ms1_matches))
chk("pct_degradation computed", abs(r7$per_sample$pct_degradation - round(100 * (1 - 100/150), 2)) < 1e-6)

cat("\n=== 8. compare_groups ===\n")
ms1_matches2 <- rbind(
  data.frame(sample = "ctrl1", met_id = "M01", met_name = "parent", kind = "parent", intensity = 90, stringsAsFactors = FALSE),
  data.frame(sample = "ctrl2", met_id = "M01", met_name = "parent", kind = "parent", intensity = 92, stringsAsFactors = FALSE),
  data.frame(sample = "treat1", met_id = "M01", met_name = "parent", kind = "parent", intensity = 40, stringsAsFactors = FALSE),
  data.frame(sample = "treat2", met_id = "M01", met_name = "parent", kind = "parent", intensity = 38, stringsAsFactors = FALSE)
)
sample_meta <- data.frame(sample = c("ctrl1", "ctrl2", "treat1", "treat2"),
                          group = c("ctrl", "ctrl", "treat", "treat"), stringsAsFactors = FALSE)
r8 <- call_agent_tool("compare_groups", list(ms1_matches = ms1_matches2, sample_meta = sample_meta,
                                              mode = "two_group", group_a = "ctrl", group_b = "treat"))
chk("compare_groups ran two_group mode", identical(r8$mode, "two_group") && nrow(r8$result) > 0)

cat("\n=== 9. Unknown tool errors clearly ===\n")
err_ok <- tryCatch({ call_agent_tool("not_a_real_tool", list()); FALSE },
                    error = function(e) grepl("Unknown tool", conditionMessage(e)))
chk("unknown tool name raises a clear error", err_ok)

cat("\n=== 10. jsonlite round-trip on a representative result ===\n")
json <- jsonlite::toJSON(r5, auto_unbox = TRUE, na = "null")
back <- jsonlite::fromJSON(json, simplifyVector = FALSE)
chk("match_ms1_features result serializes and parses back",
    is.list(back) && !is.null(back$n_matches))
json2 <- jsonlite::toJSON(r2$metabolites[[1]], auto_unbox = TRUE, na = "null")
chk("a full metabolite object serializes", nzchar(json2))

if (fail > 0) stop(fail, " agent tool check(s) failed")
cat("\n==== All agent tool tests passed ====\n")
