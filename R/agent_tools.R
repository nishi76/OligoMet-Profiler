# =============================================================================
# agent_tools.R
# Proof-of-concept tool registry for an MSAgent-style LLM agent on top of
# OligoMet-Profiler. See the design doc this implements:
#   /root/.claude/plans/staged-mapping-blanket.md
# (or ask for it again if that path doesn't exist in your environment --
# it's the "LLM Agent for OligoMet-Profiler (MSAgent-inspired)" plan).
#
# This file adds NO new analytical logic. Every tool below is a thin,
# JSON-in/JSON-out wrapper around an existing, already-tested function --
# chemistry_dict.R/oligo_io.R/metabolites.R/mass_isotope.R/fragments.R/
# ms_matching.R/degradation.R/statistics.R do the real work, unchanged.
#
# Each entry in AGENT_TOOLS is list(name, description, input_schema, run).
#   run(args, ctx = list()) -> a plain R list/data.frame, JSON-serializable
#   via jsonlite::toJSON() at whatever transport calls it (a plumber API, an
#   MCP server, or an in-process Shiny call -- none of those exist yet; this
#   file is the shared foundation the plan's Phase 1 (agent loop) and both
#   front ends (Track A: in-app chat, Track B: MCP server) are built on).
#
#   args: named list decoded from the LLM's tool-call arguments.
#   ctx:  OPTIONAL already-computed objects from earlier in the same
#         conversation (dict, spec, mets, ms1_matches) -- lets an in-process
#         caller (the planned Shiny chat tab) reuse rv$mets/rv$dict directly
#         instead of round-tripping large objects through JSON on every
#         call. A stateless caller (the planned MCP server) just passes
#         ctx = list() and supplies the same information through args
#         instead -- every tool that can read something from ctx can also
#         reconstruct it from args, so both calling conventions work against
#         the exact same tool.
#
# Every tool that has real evidence attached (ppm error, isotope fit,
# ambiguity, coverage, confirmation score) returns it verbatim and
# unrounded. This is the anti-hallucination mechanism the whole design rests
# on: an agent loop built on top of this file must have the model cite
# these fields for any claim it makes, never its own guess -- see the plan
# doc's discussion of MSAgent Fig 3a (confidence-accuracy correlation
# r = 0.438 with tool grounding vs -0.219 for an LLM alone).
# =============================================================================

## ---- Small helpers shared by several tools ---------------------------------

# A JSON array of objects decodes (via jsonlite::fromJSON) to a list of named
# lists, not a data.frame -- this reconstructs the data.frame the underlying
# analysis functions expect. Already-a-data.frame (the in-process/ctx case)
# passes through unchanged.
.rows_to_df <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  if (!is.list(x) || length(x) == 0) return(data.frame())
  do.call(rbind, lapply(x, function(row) as.data.frame(row, stringsAsFactors = FALSE)))
}

# Resolve the chemistry dictionary a tool call should use: an already-built
# dict handed through ctx (Track A reuses rv$dict this way), else build one
# fresh from any overrides in args (Track B/MCP; empty overrides = STANDARD_DICT).
.resolve_dict <- function(args, ctx) {
  if (!is.null(ctx$dict)) return(ctx$dict)
  build_dictionary(overrides = args$dictionary_overrides %||% list())
}

# Resolve the oligo_spec a tool call should act on: ctx$spec (already parsed
# earlier in this conversation) takes priority, then args$spec (a spec object
# round-tripped from a previous parse_sequence call), then args$sequence
# (parsed fresh). Errors clearly if none of the three is usable.
.resolve_spec <- function(args, ctx, dict) {
  if (!is.null(ctx$spec)) return(ctx$spec)
  if (!is.null(args$spec)) return(args$spec)
  if (!is.null(args$sequence)) {
    return(parse_input(args$sequence, dict = dict,
                       notation = args$notation %||% "auto"))
  }
  stop("No sequence available: pass 'sequence' (or 'spec' from a prior ",
       "parse_sequence call), or call this tool with session context.")
}

# Resolve the metabolite library a tool call should act on: ctx$mets first,
# else args$metabolites (a list round-tripped from a prior
# build_metabolite_library call). If only a single metabolite is needed and
# args$met_id is given, narrows to that one entry.
.resolve_mets <- function(args, ctx) {
  mets <- if (!is.null(ctx$mets)) ctx$mets else args$metabolites
  if (is.null(mets)) {
    stop("No metabolite library available: call build_metabolite_library ",
         "first (or pass 'metabolites' from its result), or call this tool ",
         "with session context.")
  }
  mets
}

.resolve_one_met <- function(args, ctx) {
  if (!is.null(args$metabolite)) return(args$metabolite)
  mets <- .resolve_mets(args, ctx)
  if (!is.null(args$met_id)) {
    idx <- which(vapply(mets, function(m) identical(m$id, args$met_id), logical(1)))
    if (length(idx) == 0) stop("met_id '", args$met_id, "' not found in the metabolite library")
    return(mets[[idx[1]]])
  }
  stop("Pass either 'metabolite' (a single metabolite object) or 'met_id' ",
       "plus a metabolite library (ctx$mets or args$metabolites).")
}

## ---- Tool: parse_sequence ---------------------------------------------------
.tool_parse_sequence <- list(
  name = "parse_sequence",
  description = paste(
    "Parse an oligonucleotide sequence (triplet notation, OligoDistiller",
    "notation, FASTA, or a structured bases/sugars/linkages/conjugate spec)",
    "into the canonical oligo_spec this whole pipeline operates on, and",
    "report its assembled molecular formula and monoisotopic mass.",
    "Call this first for any new sequence."),
  input_schema = list(
    type = "object",
    properties = list(
      sequence = list(type = "string", description = "The oligo sequence in any supported notation."),
      notation = list(type = "string", enum = c("auto", "triplet", "oligodistiller", "fasta", "structured"),
                      description = "Notation to assume; 'auto' detects it."),
      dictionary_overrides = list(type = "object",
                                  description = "Optional custom chemistry overrides, as accepted by build_dictionary(overrides=).")
    ),
    required = list("sequence")
  ),
  run = function(args, ctx = list()) {
    dict <- .resolve_dict(args, ctx)
    spec <- parse_input(args$sequence, dict = dict, notation = args$notation %||% "auto")
    f <- assemble_oligo_formula(spec$bases, spec$sugars, spec$linkages,
                                spec$conj5, spec$conj3, dict = dict)
    list(
      spec = spec,
      length = spec$n,
      formula = format_formula(f),
      mono_mass = formula_mass(f, mono = TRUE),
      avg_mass = formula_mass(f, mono = FALSE),
      description = format_spec(spec)
    )
  }
)

## ---- Tool: build_metabolite_library ------------------------------------------
.tool_build_metabolite_library <- list(
  name = "build_metabolite_library",
  description = paste(
    "Generate the deterministic theoretical metabolite library for a parsed",
    "sequence -- the parent plus 3'/5' exonuclease truncation series and",
    "endonuclease internal fragments (both free-hydroxyl and terminal-",
    "phosphate variants), de-duplicated by structure. This is NOT a guess:",
    "every species is enumerated from known nuclease chemistry, so treat",
    "the result as exhaustive for the requested truncation depth, not as a",
    "candidate short-list to second-guess."),
  input_schema = list(
    type = "object",
    properties = list(
      spec = list(type = "object", description = "A spec from parse_sequence (optional if 'sequence' given)."),
      sequence = list(type = "string", description = "Raw sequence, if 'spec' isn't already available."),
      oligo_name = list(type = "string"),
      max_3p = list(type = "integer", description = "Max 3' exonuclease truncations (default 10)."),
      max_5p = list(type = "integer", description = "Max 5' exonuclease truncations (default 10)."),
      endo = list(type = "boolean", description = "Include endonuclease fragments (default true)."),
      min_frag_len = list(type = "integer", description = "Minimum fragment length to keep (default 3).")
    ),
    required = list()
  ),
  run = function(args, ctx = list()) {
    dict <- .resolve_dict(args, ctx)
    spec <- .resolve_spec(args, ctx, dict)
    opts <- list(
      oligo_name = args$oligo_name %||% "OLIGO",
      max_3p = args$max_3p, max_5p = args$max_5p,
      endo = args$endo, min_frag_len = args$min_frag_len
    )
    opts <- opts[!vapply(opts, is.null, logical(1))]
    mets <- generate_metabolites(spec, opts = opts, dict = dict)
    list(
      metabolites = mets,
      table = metabolite_table(mets),
      n = length(mets),
      note = paste("The full 'metabolites' array is meant to be carried forward as",
                   "session state (ctx$mets) for later tool calls, not necessarily",
                   "re-read in full every turn -- 'table' is the compact summary",
                   "actually worth inspecting turn to turn.")
    )
  }
)

## ---- Tool: get_metabolite_mass ------------------------------------------------
.tool_get_metabolite_mass <- list(
  name = "get_metabolite_mass",
  description = paste(
    "Compute the exact mass, formula, PS->PO oxidation series, and charge",
    "envelope (theoretical m/z per charge state) for one metabolite from an",
    "already-built library."),
  input_schema = list(
    type = "object",
    properties = list(
      met_id = list(type = "string"),
      metabolite = list(type = "object", description = "A single metabolite object, if not resolving by met_id."),
      metabolites = list(type = "array", description = "The library to look met_id up in (optional if session context has it)."),
      max_oxid = list(type = "integer", description = "Max PS->PO oxidation events to enumerate (default 6)."),
      z_range = list(type = "array", items = list(type = "integer"), description = "Charge states (default 3-12)."),
      h_offset = list(type = "number", description = "Envelope offset; leave 0 unless matching a legacy workbook.")
    ),
    required = list()
  ),
  run = function(args, ctx = list()) {
    dict <- .resolve_dict(args, ctx)
    met <- .resolve_one_met(args, ctx)
    z_range <- if (!is.null(args$z_range)) as.integer(unlist(args$z_range)) else 3:12
    h_offset <- args$h_offset %||% 0
    info <- metabolite_mass_info(met, dict)
    list(
      met_id = met$id, name = met$name, kind = met$kind,
      formula = info$formula_str, mono_mass = info$mono_mass, avg_mass = info$avg_mass,
      n_ps = met$n_ps, n_po = met$n_po,
      oxidation_series = ps_oxidation_series(met, max_oxid = args$max_oxid %||% 6,
                                             z_range = z_range, h_offset = h_offset, dict = dict),
      charge_envelope = charge_envelope(info$mono_mass, z_range = z_range, h_offset = h_offset)
    )
  }
)

## ---- Tool: generate_fragment_ions --------------------------------------------
.tool_generate_fragment_ions <- list(
  name = "generate_fragment_ions",
  description = paste(
    "Generate theoretical McLuckey MS/MS fragment ions (a-B, w, y, b by",
    "default) for one metabolite, at the requested fragment charge states.",
    "Use this to predict what an MS2 spectrum SHOULD look like before",
    "comparing it against an acquired one with confirm_ms2."),
  input_schema = list(
    type = "object",
    properties = list(
      met_id = list(type = "string"),
      metabolite = list(type = "object"),
      metabolites = list(type = "array"),
      ion_types = list(type = "array", items = list(type = "string"),
                       description = "Default c('aB','w','y','b') -- the dominant negative-mode CID series."),
      z_range = list(type = "array", items = list(type = "integer"), description = "Fragment charge states (default 1-2)."),
      include_dz = list(type = "boolean", description = "Also include d/z ions (default false).")
    ),
    required = list()
  ),
  run = function(args, ctx = list()) {
    dict <- .resolve_dict(args, ctx)
    met <- .resolve_one_met(args, ctx)
    ion_types <- if (!is.null(args$ion_types)) unlist(args$ion_types) else c("aB", "w", "y", "b")
    z_range <- if (!is.null(args$z_range)) as.integer(unlist(args$z_range)) else 1:2
    frags <- generate_fragments(met, dict, ion_types = ion_types, z_range = z_range,
                               include_dz = isTRUE(args$include_dz))
    tab <- fragment_table(frags)
    list(fragments = tab, n = nrow(tab))
  }
)

## ---- Tool: match_ms1_features -------------------------------------------------
.tool_match_ms1_features <- list(
  name = "match_ms1_features",
  description = paste(
    "Match a metabolite library against observed MS1 features (m/z, rt,",
    "intensity) within a ppm tolerance, across charge states/adducts/PS-",
    "oxidation levels. IMPORTANT: the result carries n_candidates/ambiguous",
    "for every match -- how many DISTINCT metabolites could explain the",
    "same observed feature. Never report a match as confident without",
    "checking this; most matches at typical tolerances have more than one",
    "candidate, and that must be surfaced, not silently dropped."),
  input_schema = list(
    type = "object",
    properties = list(
      metabolites = list(type = "array"),
      features = list(type = "array", description = "List of {mz, rt, max_intensity} (or {mz, rt, intensity}) observed features."),
      ppm_tol = list(type = "number", description = "Default 10."),
      z_range = list(type = "array", items = list(type = "integer"), description = "Default 3-12."),
      adducts = list(type = "array", items = list(type = "string"), description = "Default ['H','Na','K','NH4']."),
      max_oxid = list(type = "integer", description = "Default 6."),
      h_offset = list(type = "number", description = "Default 0.")
    ),
    required = list("features")
  ),
  run = function(args, ctx = list()) {
    dict <- .resolve_dict(args, ctx)
    mets <- .resolve_mets(args, ctx)
    features <- .rows_to_df(if (!is.null(ctx$features)) ctx$features else args$features)
    if ("intensity" %in% names(features) && !("max_intensity" %in% names(features))) {
      features$max_intensity <- features$intensity
    }
    z_range <- if (!is.null(args$z_range)) as.integer(unlist(args$z_range)) else 3:12
    adducts <- if (!is.null(args$adducts)) unlist(args$adducts) else c("H", "Na", "K", "NH4")
    result <- match_ms1(mets, features, dict, ppm_tol = args$ppm_tol %||% 10,
                        z_range = z_range, adducts = adducts,
                        max_oxid = args$max_oxid %||% 6, h_offset = args$h_offset %||% 0)
    list(
      matches = result, n_matches = nrow(result),
      n_ambiguous = if (nrow(result) > 0) sum(result$ambiguous) else 0L,
      note = "n_ambiguous counts MATCH ROWS flagged ambiguous, not distinct metabolites -- see each row's own n_candidates."
    )
  }
)

## ---- Tool: confirm_ms2 --------------------------------------------------------
.tool_confirm_ms2 <- list(
  name = "confirm_ms2",
  description = paste(
    "Confirm one metabolite against an acquired MS2 spectrum: matches",
    "theoretical fragment ions to observed peaks, checks for PS diagnostic",
    "ions, and returns a coverage/confirmation score. This score, not a",
    "self-reported confidence, is what any claim of 'confirmed by MS2'",
    "must be grounded in."),
  input_schema = list(
    type = "object",
    properties = list(
      met_id = list(type = "string"),
      metabolite = list(type = "object"),
      metabolites = list(type = "array"),
      ms2_peaks = list(type = "array", description = "List of {mz, intensity} observed MS2 peaks."),
      tol_ppm = list(type = "number", description = "Default 25."),
      z_range = list(type = "array", items = list(type = "integer"), description = "Default 1-2.")
    ),
    required = list("ms2_peaks")
  ),
  run = function(args, ctx = list()) {
    dict <- .resolve_dict(args, ctx)
    met <- .resolve_one_met(args, ctx)
    ms2 <- .rows_to_df(args$ms2_peaks)
    z_range <- if (!is.null(args$z_range)) as.integer(unlist(args$z_range)) else 1:2
    conf <- confirm_metabolite(met, ms2, dict, tol_ppm = args$tol_ppm %||% 25, z_range = z_range)
    list(
      met_id = met$id, name = met$name,
      score = conf$score, coverage = conf$coverage,
      n_matched = nrow(conf$matched), matched_fragments = conf$matched,
      diagnostics = conf$diagnostics
    )
  }
)

## ---- Tool: summarize_degradation ----------------------------------------------
.tool_summarize_degradation <- list(
  name = "summarize_degradation",
  description = paste(
    "Compute % degradation of the parent oligo per sample from matched MS1",
    "peak areas/intensities, plus the composition breakdown by metabolite",
    "class (exo_3p/exo_5p/endo_*) and the top degradants ranked by signal."),
  input_schema = list(
    type = "object",
    properties = list(
      ms1_matches = list(type = "array", description = "The 'matches' array from match_ms1_features (with a 'sample' column for batch data)."),
      top_n = list(type = "integer", description = "Default 10.")
    ),
    required = list("ms1_matches")
  ),
  run = function(args, ctx = list()) {
    m <- .rows_to_df(if (!is.null(ctx$ms1_matches)) ctx$ms1_matches else args$ms1_matches)
    deg <- degradation_summary(m, top_n = args$top_n %||% 10)
    list(per_sample = deg$per_sample, composition = deg$composition,
        top_degradants = deg$top_degradants, signal_used = deg$signal_used)
  }
)

## ---- Tool: compare_groups ------------------------------------------------------
.tool_compare_groups <- list(
  name = "compare_groups",
  description = paste(
    "Statistically compare matched metabolite signal across experimental",
    "groups (two-group Welch t-test, multi-group ANOVA, or a time-series",
    "trend), with Benjamini-Hochberg correction. Requires per-sample group/",
    "timepoint metadata."),
  input_schema = list(
    type = "object",
    properties = list(
      ms1_matches = list(type = "array"),
      sample_meta = list(type = "array", description = "List of {sample, group} or {sample, timepoint}."),
      mode = list(type = "string", enum = c("two_group", "multi_group", "time_series")),
      group_a = list(type = "string", description = "Required when mode = 'two_group'."),
      group_b = list(type = "string", description = "Required when mode = 'two_group'.")
    ),
    required = list("ms1_matches", "sample_meta", "mode")
  ),
  run = function(args, ctx = list()) {
    m <- .rows_to_df(if (!is.null(ctx$ms1_matches)) ctx$ms1_matches else args$ms1_matches)
    meta <- .rows_to_df(args$sample_meta)
    abund <- build_abundance_matrix(m)
    long <- abundance_long(abund, meta)
    result <- switch(args$mode,
      two_group = compare_two_groups(long, args$group_a, args$group_b),
      multi_group = compare_multi_groups(long),
      time_series = compare_time_series(long),
      stop("mode must be one of two_group/multi_group/time_series"))
    list(mode = args$mode, result = result)
  }
)

## ---- Registry ------------------------------------------------------------------
AGENT_TOOLS <- list(
  .tool_parse_sequence,
  .tool_build_metabolite_library,
  .tool_get_metabolite_mass,
  .tool_generate_fragment_ions,
  .tool_match_ms1_features,
  .tool_confirm_ms2,
  .tool_summarize_degradation,
  .tool_compare_groups
)
names(AGENT_TOOLS) <- vapply(AGENT_TOOLS, function(t) t$name, character(1))

# Dispatch by name -- the one entry point both front ends (plumber/MCP and
# the in-process Shiny caller) are expected to call through, so validation
# and error formatting stay in one place.
call_agent_tool <- function(name, args = list(), ctx = list()) {
  tool <- AGENT_TOOLS[[name]]
  if (is.null(tool)) {
    stop("Unknown tool '", name, "'. Available: ", paste(names(AGENT_TOOLS), collapse = ", "))
  }
  tool$run(args, ctx)
}

# JSON Schema for every tool's parameters, in the shape Anthropic/OpenAI tool-
# use expects (name/description/input_schema per tool) -- Phase 1's agent
# core hands this straight to whichever backend it calls.
agent_tool_specs <- function() {
  lapply(AGENT_TOOLS, function(t) list(name = t$name, description = t$description,
                                       input_schema = t$input_schema))
}
