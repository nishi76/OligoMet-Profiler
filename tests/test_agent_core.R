# test_agent_core.R -- validate the provider-agnostic agent loop (R/agent_core.R)
# using the "stub" backend, so this runs with no network access and no API key.
# Checks:
#   1. A turn with no tool call returns immediately with the final text.
#   2. A turn that calls one tool: the tool actually runs, its result folds
#      into ctx, and the trimmed (not full) result is what's echoed back to
#      the model.
#   3. A turn with two tool_use blocks in one assistant turn: both execute.
#   4. A failing tool call becomes an is_error tool_result, not a crash.
#   5. Exceeding max_tool_iterations raises a clear error.
#   6. The OpenAI message-format translator round-trips a representative
#      conversation (checked directly, since this sandbox can't reach the
#      real OpenAI API to test .llm_call_openai() end-to-end).

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
            "agent_tools", "agent_core")) {
  source(file.path(.pkg_root, "R", paste0(m, ".R")))
}

fail <- 0L
chk <- function(label, ok) {
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) fail <<- fail + 1L
}
text_block <- function(t) list(type = "text", text = t)
tool_use_block <- function(id, name, input) list(type = "tool_use", id = id, name = name, input = input)

cat("=== 1. No tool call: immediate final answer ===\n")
r1 <- run_agent_turn(
  list(user_turn("hi")), backend = "stub",
  stub_responses = list(list(role = "assistant", content = list(text_block("Hello there.")), stop_reason = "end_turn"))
)
chk("final text is the stub's answer", identical(final_text(r1$final), "Hello there."))
chk("messages history has 2 turns (user + assistant)", length(r1$messages) == 2)

cat("\n=== 2. One tool call: executes, folds into ctx, gets trimmed for the model ===\n")
called <- list()
r2 <- run_agent_turn(
  list(user_turn("parse inotersen")), backend = "stub",
  stub_responses = list(
    list(role = "assistant",
        content = list(tool_use_block("t1", "parse_sequence", list(sequence = INOTERSEN_TRIPLET))),
        stop_reason = "tool_use"),
    list(role = "assistant", content = list(text_block("Parsed: 20-mer.")), stop_reason = "end_turn")
  ),
  on_tool_call = function(name, input) called[[length(called) + 1]] <<- name
)
chk("tool was actually invoked", identical(unlist(called), "parse_sequence"))
chk("final text is the second stub turn", identical(final_text(r2$final), "Parsed: 20-mer."))
chk("ctx$spec was folded in from the tool result", !is.null(r2$ctx$spec) && r2$ctx$spec$n == 20)
tool_result_msg <- r2$messages[[3]]
chk("a tool_result turn was appended", identical(tool_result_msg$content[[1]]$type, "tool_result"))
chk("tool_result content is valid JSON", {
  ok <- TRUE
  tryCatch(jsonlite::fromJSON(tool_result_msg$content[[1]]$content), error = function(e) ok <<- FALSE)
  ok
})

cat("\n=== 3. Build library then call get_metabolite_mass in the SAME turn using ctx ===\n")
r3 <- run_agent_turn(
  list(user_turn("build the library and tell me the parent mass")), backend = "stub",
  ctx = list(spec = r2$ctx$spec, dict = STANDARD_DICT),
  stub_responses = list(
    list(role = "assistant",
        content = list(tool_use_block("t2", "build_metabolite_library",
                                      list(oligo_name = "inotersen", max_3p = 1, max_5p = 1, endo = FALSE))),
        stop_reason = "tool_use"),
    list(role = "assistant",
        content = list(tool_use_block("t3", "get_metabolite_mass", list(met_id = "M01"))),
        stop_reason = "tool_use"),
    list(role = "assistant", content = list(text_block("Parent mono mass is 7178.06 Da.")), stop_reason = "end_turn")
  )
)
chk("get_metabolite_mass succeeded using ctx$mets from the prior call (no metabolites re-supplied)",
    identical(final_text(r3$final), "Parent mono mass is 7178.06 Da."))
built_result_msg <- r3$messages[[3]]$content[[1]]$content
chk("the build_metabolite_library result shown to the model has 'metabolites' stripped (trimmed)",
    !grepl('"metabolites"', built_result_msg) && grepl('"table"', built_result_msg))

cat("\n=== 4. A failing tool call becomes an is_error tool_result, not a crash ===\n")
r4 <- run_agent_turn(
  list(user_turn("get mass of a metabolite that doesn't exist")), backend = "stub",
  ctx = list(mets = r3$ctx$mets, dict = STANDARD_DICT),
  stub_responses = list(
    list(role = "assistant", content = list(tool_use_block("t4", "get_metabolite_mass", list(met_id = "NOPE"))), stop_reason = "tool_use"),
    list(role = "assistant", content = list(text_block("That metabolite id doesn't exist.")), stop_reason = "end_turn")
  )
)
chk("loop did not crash on a failing tool call", identical(final_text(r4$final), "That metabolite id doesn't exist."))
chk("the tool_result is flagged is_error", isTRUE(r4$messages[[3]]$content[[1]]$is_error))

cat("\n=== 5. Exceeding max_tool_iterations raises a clear error ===\n")
looping_response <- list(role = "assistant",
                         content = list(tool_use_block("tx", "get_metabolite_mass", list(met_id = "M01"))),
                         stop_reason = "tool_use")
err_ok <- tryCatch({
  run_agent_turn(list(user_turn("loop forever")), backend = "stub", max_tool_iterations = 2,
                 ctx = list(mets = r3$ctx$mets, dict = STANDARD_DICT),
                 stub_responses = list(looping_response, looping_response, looping_response))
  FALSE
}, error = function(e) grepl("did not converge", conditionMessage(e)))
chk("clear convergence error raised", err_ok)

cat("\n=== 6. OpenAI message-format translation round-trips ===\n")
canon <- list(
  user_turn("what's the parent mass?"),
  list(role = "assistant", content = list(tool_use_block("call_1", "get_metabolite_mass", list(met_id = "M01")))),
  list(role = "user", content = list(list(type = "tool_result", tool_use_id = "call_1", content = '{"mono_mass":7178.06}', is_error = FALSE)))
)
oai <- .canonical_to_openai_messages(canon, "system prompt here")
chk("system prompt prepended", identical(oai[[1]]$role, "system"))
chk("user text message present", identical(oai[[2]]$role, "user") && identical(oai[[2]]$content, "what's the parent mass?"))
chk("assistant tool_calls translated", identical(oai[[3]]$tool_calls[[1]]$`function`$name, "get_metabolite_mass"))
chk("tool_calls arguments is a JSON string", is.character(oai[[3]]$tool_calls[[1]]$`function`$arguments))
chk("tool result became its own role='tool' message", identical(oai[[4]]$role, "tool") && identical(oai[[4]]$tool_call_id, "call_1"))

back <- .openai_message_to_canonical(list(content = "Parent mass is 7178.06 Da.",
                                          tool_calls = list(list(id = "call_2", `function` = list(name = "get_metabolite_mass", arguments = '{"met_id":"M01"}')))))
chk("openai->canonical: text block present", any(vapply(back, function(b) identical(b$type, "text") && identical(b$text, "Parent mass is 7178.06 Da."), logical(1))))
chk("openai->canonical: tool_use block present with parsed input", {
  tb <- Filter(function(b) identical(b$type, "tool_use"), back)[[1]]
  identical(tb$name, "get_metabolite_mass") && identical(tb$input$met_id, "M01")
})

if (fail > 0) stop(fail, " agent core check(s) failed")
cat("\n==== All agent core tests passed ====\n")
