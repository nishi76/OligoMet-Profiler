# =============================================================================
# agent_core.R
# Phase 1 of the MSAgent-inspired LLM agent (see /root/.claude/plans/
# staged-mapping-blanket.md). Provider-agnostic ReAct-style agent loop --
# think (call an LLM) -> act (invoke a tool from R/agent_tools.R) ->
# observe (feed the result back) -> repeat until the model gives a final
# text answer -- built once here and shared by both front ends (the Shiny
# "Ask OligoMet" tab and the MCP server), same as the tool registry itself.
#
# Deviation from the plan's sketch, flagged explicitly: the plan mentioned
# httr2 for the HTTP calls. httr2 (and httr) are not installable in this
# development sandbox (no CRAN access, not packaged for apt here), so this
# shells out to the `curl` binary instead via .http_post_json() below --
# curl is essentially always present, this adds no new R dependency at all,
# and every backend adapter is written against that one helper, so swapping
# in httr2 later (if a deployment has it) only means replacing that one
# function, not any of the provider translation logic.
#
# ---- Canonical message format --------------------------------------------
# One turn = list(role = "user" | "assistant", content = list(<blocks>)).
# A block is one of:
#   list(type = "text", text = "...")
#   list(type = "tool_use", id = "...", name = "...", input = list(...))
#   list(type = "tool_result", tool_use_id = "...", content = "...", is_error = FALSE)
# This is (deliberately) Anthropic's own Messages API content-block shape --
# the richest of the two, and translating it to/from OpenAI's format (a
# flatter role="tool" message per result, tool_calls on the assistant
# message) is a one-way mechanical conversion either direction.
# =============================================================================

## ---- Dependency-free HTTP POST (JSON in, JSON out) --------------------------
# Header values (notably the x-api-key secret) are passed to curl via a -K
# config file, never as literal argv entries. Two independent reasons:
# argv is visible to other users on the same machine via ps/Task Manager,
# and -- the reason this actually leaked in practice -- system2() itself
# prints the FULL command line, api key included, in an R warning whenever
# curl exits non-zero. A transient network failure (DNS hiccup, timeout,
# offline machine) would otherwise dump a live secret straight to the
# R/RStudio console and any log or bug report that happens to capture it.
# suppressWarnings() below additionally blocks that auto-emitted warning
# outright, since the exit status is already checked explicitly and a
# sanitized stop() message is raised instead.
# Quotes a value for a curl -K config file. Escapes backslash BEFORE
# quote, and escapes backslash at all -- a Windows tempfile() path
# (tmp_body/tmp_out in .http_post_json() below) is backslash-separated,
# and curl's own -K config-file parser treats an unescaped backslash
# inside a quoted value as an escape character (dropping it, per curl's
# documented "\\, \", \t, \n, \r, \v are the only recognized escapes; a
# backslash before any other letter is ignored"), silently corrupting
# the path into a single run-on word (confirmed by reproducing the exact
# reported failure -- "curl: Failed to open ...file.json" with every
# backslash missing -- against a real curl -K config, and confirming the
# fix restores the original path via a real file whose name contains a
# literal backslash). A standalone (not nested-in-closure) function so
# it's directly unit-testable -- see tests/test_agent_core.R.
.cfg_quote <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  paste0("\"", x, "\"")
}

.http_post_json <- function(url, headers = list(), body, timeout_sec = 120) {
  body_json <- jsonlite::toJSON(body, auto_unbox = TRUE, na = "null", digits = NA)
  tmp_body <- tempfile(fileext = ".json")
  tmp_out <- tempfile()
  tmp_cfg <- tempfile(fileext = ".curlcfg")
  on.exit(unlink(c(tmp_body, tmp_out, tmp_cfg)), add = TRUE)
  writeLines(body_json, tmp_body, useBytes = TRUE)

  all_headers <- c(headers, list(`Content-Type` = "application/json"))
  header_lines <- sprintf("header = %s",
    .cfg_quote(paste0(names(all_headers), ": ", unlist(all_headers))))
  cfg_lines <- c(
    header_lines,
    paste("url =", .cfg_quote(url)),
    "silent", "show-error",
    paste("max-time =", as.integer(timeout_sec)),
    "request = \"POST\"",
    paste("data-binary =", .cfg_quote(paste0("@", tmp_body))),
    paste("output =", .cfg_quote(tmp_out)),
    "write-out = \"%{http_code}\""
  )
  writeLines(cfg_lines, tmp_cfg)

  status_lines <- suppressWarnings(
    system2("curl", args = c("-K", tmp_cfg), stdout = TRUE, stderr = TRUE))
  exit_code <- attr(status_lines, "status")
  http_code <- suppressWarnings(as.integer(status_lines[length(status_lines)]))
  resp_text <- if (file.exists(tmp_out)) paste(readLines(tmp_out, warn = FALSE, encoding = "UTF-8"), collapse = "\n") else ""
  # curl prints "000" for --write-out's %{http_code} whenever no HTTP
  # response was ever received at all (DNS failure, connection refused,
  # TLS handshake failure, a corporate proxy/firewall silently dropping the
  # connection, a timeout before headers arrive, ...) -- as.integer("000")
  # is 0, NOT NA, so this must be checked separately from the NA/malformed
  # case below or it silently falls through as if it were a real, empty,
  # HTTP 0 response (confirmed by reproducing this exact path against an
  # unreachable host: it previously returned status=0/body="" instead of
  # raising here, which then surfaced downstream as a confusing
  # "unparseable JSON (HTTP 0)" error with no hint that the real problem was
  # never reaching the server in the first place).
  if (is.na(http_code) || http_code == 0L || (!is.null(exit_code) && exit_code != 0L)) {
    msg_lines <- if (!is.na(http_code)) status_lines[-length(status_lines)] else status_lines
    diag <- paste(Filter(nzchar, msg_lines), collapse = " ")
    stop("Could not reach ", url, " (curl exit status ",
         if (is.null(exit_code)) "unknown" else exit_code, "): ",
         if (nzchar(diag)) diag else "no response received.",
         " Check network access, proxy settings, and any corporate firewall ",
         "that might block outbound HTTPS to this host.")
  }
  list(status = http_code, body = resp_text)
}

.default_api_key <- function(backend) {
  key <- switch(backend,
    anthropic = Sys.getenv("ANTHROPIC_API_KEY", ""),
    openai = Sys.getenv("OPENAI_API_KEY", ""),
    stop("Unknown backend: ", backend))
  if (!nzchar(key)) {
    stop("No API key set for backend '", backend, "'. Set ",
        switch(backend, anthropic = "ANTHROPIC_API_KEY", openai = "OPENAI_API_KEY"),
        " in the environment (never type it into the app UI).")
  }
  key
}

.default_model <- function(backend) {
  # Overridable via env var so a stale/inaccessible model name doesn't need
  # a code change to work around -- e.g. a key on an older plan/org without
  # access to the current default can set OLIGOMET_ANTHROPIC_MODEL to
  # whatever model that key does have.
  switch(backend,
    anthropic = Sys.getenv("OLIGOMET_ANTHROPIC_MODEL", "claude-sonnet-5"),
    openai = Sys.getenv("OLIGOMET_OPENAI_MODEL", "gpt-4o"),
    stop("Unknown backend: ", backend))
}

## ---- Anthropic Messages API adapter ------------------------------------------
# Canonical messages/tools are already in (almost exactly) Anthropic's own
# shape, so this is close to a passthrough plus the HTTP/auth wiring.
.llm_call_anthropic <- function(messages, tools, system_prompt, model, api_key,
                                max_tokens = 4096) {
  tool_specs <- lapply(tools, function(t) list(name = t$name, description = t$description,
                                               input_schema = t$input_schema))
  body <- list(model = model, max_tokens = max_tokens, system = system_prompt,
              messages = messages, tools = tool_specs)
  resp <- .http_post_json("https://api.anthropic.com/v1/messages",
                          headers = list(`x-api-key` = api_key, `anthropic-version` = "2023-06-01"),
                          body = body)
  parsed <- tryCatch(jsonlite::fromJSON(resp$body, simplifyVector = FALSE),
                     error = function(e) stop("Anthropic API returned unparseable JSON (HTTP ",
                                              resp$status, "): ", substr(resp$body, 1, 500)))
  if (resp$status >= 400) {
    msg <- if (!is.null(parsed$error$message)) parsed$error$message else resp$body
    stop("Anthropic API error (HTTP ", resp$status, "): ", msg)
  }
  list(role = "assistant", content = parsed$content,
      stop_reason = if (identical(parsed$stop_reason, "tool_use")) "tool_use" else "end_turn")
}

## ---- OpenAI Chat Completions adapter -----------------------------------------
.canonical_to_openai_messages <- function(messages, system_prompt) {
  out <- list(list(role = "system", content = system_prompt))
  for (m in messages) {
    if (identical(m$role, "assistant")) {
      text_blocks <- Filter(function(b) identical(b$type, "text"), m$content)
      tool_blocks <- Filter(function(b) identical(b$type, "tool_use"), m$content)
      entry <- list(role = "assistant",
                    content = if (length(text_blocks) > 0) paste(vapply(text_blocks, function(b) b$text, ""), collapse = "\n") else NULL)
      if (length(tool_blocks) > 0) {
        entry$tool_calls <- lapply(tool_blocks, function(b) list(
          id = b$id, type = "function",
          `function` = list(name = b$name, arguments = as.character(jsonlite::toJSON(b$input, auto_unbox = TRUE, na = "null")))
        ))
      }
      out[[length(out) + 1]] <- entry
    } else {
      # role == "user": either plain text turns, or a batch of tool_result
      # blocks -- OpenAI wants one role="tool" message PER result, not
      # grouped the way Anthropic's tool_result blocks are.
      tool_results <- Filter(function(b) identical(b$type, "tool_result"), m$content)
      text_blocks <- Filter(function(b) identical(b$type, "text"), m$content)
      if (length(tool_results) > 0) {
        for (tr in tool_results) {
          out[[length(out) + 1]] <- list(role = "tool", tool_call_id = tr$tool_use_id, content = tr$content)
        }
      } else if (length(text_blocks) > 0) {
        out[[length(out) + 1]] <- list(role = "user", content = paste(vapply(text_blocks, function(b) b$text, ""), collapse = "\n"))
      }
    }
  }
  out
}

.openai_message_to_canonical <- function(msg) {
  blocks <- list()
  if (!is.null(msg$content) && nzchar(msg$content %||% "")) {
    blocks[[length(blocks) + 1]] <- list(type = "text", text = msg$content)
  }
  for (tc in msg$tool_calls %||% list()) {
    blocks[[length(blocks) + 1]] <- list(
      type = "tool_use", id = tc$id, name = tc$`function`$name,
      input = jsonlite::fromJSON(tc$`function`$arguments, simplifyVector = FALSE)
    )
  }
  blocks
}

.llm_call_openai <- function(messages, tools, system_prompt, model, api_key, max_tokens = 4096) {
  tool_specs <- lapply(tools, function(t) list(type = "function", `function` = list(
    name = t$name, description = t$description, parameters = t$input_schema)))
  body <- list(model = model, max_tokens = max_tokens,
              messages = .canonical_to_openai_messages(messages, system_prompt),
              tools = tool_specs, tool_choice = "auto")
  resp <- .http_post_json("https://api.openai.com/v1/chat/completions",
                          headers = list(Authorization = paste("Bearer", api_key)), body = body)
  parsed <- tryCatch(jsonlite::fromJSON(resp$body, simplifyVector = FALSE),
                     error = function(e) stop("OpenAI API returned unparseable JSON (HTTP ",
                                              resp$status, "): ", substr(resp$body, 1, 500)))
  if (resp$status >= 400) {
    msg <- if (!is.null(parsed$error$message)) parsed$error$message else resp$body
    stop("OpenAI API error (HTTP ", resp$status, "): ", msg)
  }
  choice <- parsed$choices[[1]]
  list(role = "assistant", content = .openai_message_to_canonical(choice$message),
      stop_reason = if (identical(choice$finish_reason, "tool_calls")) "tool_use" else "end_turn")
}

## ---- Grounding discipline: the system prompt every backend gets -------------
AGENT_SYSTEM_PROMPT <- paste(
  "You are an assistant for OligoMet Profiler, an oligonucleotide metabolite",
  "identification pipeline. Every structural species you could ever discuss",
  "is already deterministically enumerated by this pipeline's own chemistry",
  "(known exo/endonuclease cleavage, PS<->PO oxidation, adducts) -- you are",
  "never guessing at an unknown molecule's structure the way a general",
  "small-molecule tool would. Your job is to orchestrate the available tools",
  "and narrate their output faithfully, not to reason about mass spectra from",
  "first principles or general chemistry knowledge.",
  "",
  "Hard rules:",
  "1. Never state a confidence, identification, or quantitative claim that",
  "   isn't directly backed by a field a tool actually returned (ppm_error,",
  "   iso_fit, n_candidates/ambiguous, coverage, confirmation score,",
  "   pct_degradation, p.value/padj). If you don't have the evidence, say so",
  "   and call the tool that would produce it, or tell the user what's missing.",
  "2. Always name the specific evidence behind a claim in your answer (e.g.",
  "   '<10 ppm, isotope fit 0.97, single candidate' or '3 candidate species",
  "   within 10 ppm -- ambiguous, needs MS2 or RT to resolve'), not just a",
  "   verdict.",
  "3. Never silently pick among ambiguous candidates. If match_ms1_features",
  "   reports ambiguous = TRUE / n_candidates > 1 for a match, say so",
  "   explicitly and name the alternatives if you have them.",
  "4. Prefer tools over memory for anything sequence- or spectrum-specific.",
  "   General background on oligonucleotide modalities (e.g. typical gapmer",
  "   catabolism routes) is fine to mention as context, clearly labeled as",
  "   general background, never as if it were this specific result.",
  sep = "\n"
)

## ---- The agent loop -----------------------------------------------------------
# One user message in, drives think->act->observe until the model returns a
# final text answer (or max_tool_iterations is hit). Returns:
#   list(messages = <full updated history>, final = <last assistant turn>,
#        ctx = <updated session context>)
#
# ctx accumulates automatically: after a successful build_metabolite_library/
# parse_sequence/match_ms1_features call, its heavy result fields (spec,
# metabolites, matches) are folded into ctx so later tool calls in this same
# turn -- and, if the caller persists ctx across turns (the Shiny tab does,
# via a reactiveValue) -- later USER turns too, don't need the model to
# re-supply them. The copy of that same result actually shown to the model
# has those heavy fields stripped first (see .trim_for_llm()), so the
# model's own context window only ever sees the compact summary fields.
#
# backend: "anthropic" | "openai" | "stub". "stub" takes `stub_responses`. a
# list of pre-scripted assistant turns (canonical format) returned in order,
# one per iteration -- for tests, so the loop can be verified with no
# network access and no API key.
run_agent_turn <- function(messages, ctx = list(), backend = "anthropic",
                            model = NULL, api_key = NULL,
                            system_prompt = AGENT_SYSTEM_PROMPT,
                            max_tool_iterations = 8,
                            stub_responses = NULL, on_tool_call = NULL) {
  call_fn <- switch(backend,
    anthropic = function(msgs) .llm_call_anthropic(msgs, AGENT_TOOLS, system_prompt,
                                                    model %||% .default_model("anthropic"),
                                                    api_key %||% .default_api_key("anthropic")),
    openai = function(msgs) .llm_call_openai(msgs, AGENT_TOOLS, system_prompt,
                                             model %||% .default_model("openai"),
                                             api_key %||% .default_api_key("openai")),
    stub = local({
      i <- 0L
      function(msgs) { i <<- i + 1L; stub_responses[[i]] }
    }),
    stop("Unknown backend '", backend, "': use 'anthropic', 'openai', or 'stub'.")
  )

  for (iter in seq_len(max_tool_iterations)) {
    resp <- call_fn(messages)
    messages <- c(messages, list(list(role = "assistant", content = resp$content)))
    tool_use_blocks <- Filter(function(b) identical(b$type, "tool_use"), resp$content)
    if (length(tool_use_blocks) == 0 || !identical(resp$stop_reason, "tool_use")) {
      return(list(messages = messages, final = resp, ctx = ctx))
    }

    result_blocks <- lapply(tool_use_blocks, function(tb) {
      if (!is.null(on_tool_call)) on_tool_call(tb$name, tb$input)
      out <- tryCatch(list(ok = TRUE, value = call_agent_tool(tb$name, tb$input, ctx)),
                      error = function(e) list(ok = FALSE, value = conditionMessage(e)))
      if (out$ok) ctx <<- .fold_into_ctx(ctx, out$value)
      shown <- if (out$ok) .trim_for_llm(tb$name, out$value) else out$value
      content_str <- if (out$ok) as.character(jsonlite::toJSON(shown, auto_unbox = TRUE, na = "null", digits = NA))
                      else paste("ERROR:", out$value)
      list(type = "tool_result", tool_use_id = tb$id, content = content_str, is_error = !out$ok)
    })
    messages <- c(messages, list(list(role = "user", content = result_blocks)))
  }
  stop("Agent did not converge to a final answer within ", max_tool_iterations, " tool-use iterations.")
}

# Accumulate a tool result's reusable fields into the running session
# context, so later tool calls (same turn or, if the caller persists ctx,
# later turns) can use ctx$spec/ctx$mets/ctx$ms1_matches instead of the model
# re-supplying them.
.fold_into_ctx <- function(ctx, value) {
  if (!is.null(value$spec)) ctx$spec <- value$spec
  if (!is.null(value$metabolites)) ctx$mets <- value$metabolites
  if (!is.null(value$matches)) ctx$ms1_matches <- value$matches
  ctx
}

# Strip fields from a tool result before it's echoed back into the model's
# own context window, for tools whose full result is large and already
# folded into ctx by .fold_into_ctx() above -- the model sees the compact
# summary (table/n/note) and works off ctx for anything that needs the rest.
.trim_for_llm <- function(name, value) {
  if (identical(name, "build_metabolite_library")) value$metabolites <- NULL
  value
}

## ---- Convenience: build the first user turn ----------------------------------
user_turn <- function(text) list(role = "user", content = list(list(type = "text", text = text)))

# Extract the plain-text answer from a final assistant turn (for a caller
# that just wants the reply, not the full content-block structure).
final_text <- function(final) {
  blocks <- Filter(function(b) identical(b$type, "text"), final$content)
  paste(vapply(blocks, function(b) b$text, ""), collapse = "\n")
}
