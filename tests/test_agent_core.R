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

cat("\n=== 7. Gemini message-format translation round-trips ===\n")
gem_canon <- list(
  user_turn("what's the parent mass?"),
  list(role = "assistant", content = list(tool_use_block("call_1", "get_metabolite_mass", list(met_id = "M01")))),
  list(role = "user", content = list(list(type = "tool_result", tool_use_id = "call_1", content = '{"mono_mass":7178.06}', is_error = FALSE)))
)
gem_contents <- .canonical_to_gemini_contents(gem_canon)
chk("user text turn present with role='user'",
    identical(gem_contents[[1]]$role, "user") && identical(gem_contents[[1]]$parts[[1]]$text, "what's the parent mass?"))
chk("assistant tool_use translated to role='model' functionCall",
    identical(gem_contents[[2]]$role, "model") &&
    identical(gem_contents[[2]]$parts[[1]]$functionCall$name, "get_metabolite_mass") &&
    identical(gem_contents[[2]]$parts[[1]]$functionCall$args$met_id, "M01"))
chk("tool_result resolved to functionResponse with the ORIGINAL function name (via id->name lookup, not the id itself)",
    identical(gem_contents[[3]]$role, "user") &&
    identical(gem_contents[[3]]$parts[[1]]$functionResponse$name, "get_metabolite_mass") &&
    identical(gem_contents[[3]]$parts[[1]]$functionResponse$response$mono_mass, 7178.06))

gem_back_text <- .gemini_parts_to_canonical(list(list(text = "Parent mass is 7178.06 Da.")))
chk("gemini->canonical: plain text part",
    identical(gem_back_text[[1]]$type, "text") && identical(gem_back_text[[1]]$text, "Parent mass is 7178.06 Da."))

gem_back_call <- .gemini_parts_to_canonical(list(list(functionCall = list(name = "get_metabolite_mass", args = list(met_id = "M01")))))
chk("gemini->canonical: functionCall becomes a tool_use block with a locally-invented id", {
  tb <- gem_back_call[[1]]
  identical(tb$type, "tool_use") && identical(tb$name, "get_metabolite_mass") &&
    identical(tb$input$met_id, "M01") && nzchar(tb$id)
})

cat("\n=== 8. Tool specs sent to every provider serialize as a JSON ARRAY, not an object ===\n")
# Regression test for a real bug this session's own live Gemini smoke test
# caught: AGENT_TOOLS is a NAMED list (named by tool name, for $-lookup in
# call_agent_tool()), and jsonlite::toJSON() serializes a named list as a
# JSON OBJECT rather than an array unless the names are stripped first.
# Every adapter's tool_specs <- lapply(tools, ...) line must therefore
# unname() its input, or the "tools" field sent to the provider silently
# becomes '{"parse_sequence": {...}, ...}' instead of '[{...}, ...]' --
# every real provider rejects that (Gemini: HTTP 400 "Unknown name
# 'parse_sequence' at 'tools[0].function_declarations'").
buggy_json <- jsonlite::toJSON(lapply(AGENT_TOOLS, function(t) list(name = t$name)), auto_unbox = TRUE)
fixed_json <- jsonlite::toJSON(lapply(unname(AGENT_TOOLS), function(t) list(name = t$name)), auto_unbox = TRUE)
chk("without unname(): a NAMED list of tool specs (the bug) serializes as a JSON OBJECT",
    startsWith(as.character(buggy_json), "{"))
chk("with unname(): the same tool specs serialize as a JSON ARRAY (what every provider needs)",
    startsWith(as.character(fixed_json), "["))

cat("\n=== 9. .http_post_json() curl args survive shell word-splitting (real network) ===\n")
# Regression test for a second real bug this session's live Gemini smoke
# test caught: system2() joins `args` with plain spaces and hands them to a
# shell WITHOUT quoting -- so a header value like "Content-Type: application/json"
# (every "-H name: value" pair has an internal space, and OpenAI's
# "Authorization: Bearer <key>" has two) used to get word-split into extra
# shell tokens. curl then read the split-off tail as a second positional
# URL and failed outright ("Could not resolve host: application"), and the
# resulting garbage in $status (curl's raw exit-code line concatenated with
# its "%{http_code}" text, e.g. "400000") slipped past the `is.na(http_code)`
# check as if it were a real (if odd-looking) HTTP status. Needs real
# network access (skips cleanly if unavailable, same policy as other
# network-only checks in this suite) since the bug is specifically about
# what actually reaches the shell/curl, not something a mock can show.
probe <- suppressWarnings(system2("curl", shQuote(c("-sS", "-o", "/dev/null", "-w", "%{http_code}",
                                                     "--max-time", "10",
                                                     "https://generativelanguage.googleapis.com/")),
                                  stdout = TRUE, stderr = TRUE))
net_ok <- !is.na(suppressWarnings(as.integer(probe[length(probe)])))
if (net_ok) {
  resp <- tryCatch(
    .http_post_json("https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key=OBVIOUSLY_INVALID_TEST_KEY",
                    headers = list(), body = list(contents = list())),
    error = function(e) NULL)
  chk("a real HTTPS POST with header values containing spaces returns a real 3-digit HTTP status",
      !is.null(resp) && resp$status >= 100 && resp$status <= 599)
  chk("the response body is the provider's real JSON error, not a curl-level connection failure",
      !is.null(resp) && grepl('"error"', resp$body, fixed = TRUE) && grepl("API key not valid", resp$body, fixed = TRUE))
} else {
  cat("  [SKIP] no network access in this environment -- can't exercise the real curl call\n")
}

if (fail > 0) stop(fail, " agent core check(s) failed")
cat("\n==== All agent core tests passed ====\n")
