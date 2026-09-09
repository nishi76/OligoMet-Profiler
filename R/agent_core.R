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
.http_post_json <- function(url, headers = list(), body, timeout_sec = 120) {
  body_json <- jsonlite::toJSON(body, auto_unbox = TRUE, na = "null", digits = NA)
  tmp_body <- tempfile(fileext = ".json")
  tmp_out <- tempfile()
  on.exit(unlink(c(tmp_body, tmp_out)), add = TRUE)
  writeLines(body_json, tmp_body, useBytes = TRUE)
  header_args <- unlist(lapply(names(headers), function(h) c("-H", paste0(h, ": ", headers[[h]]))))
  args <- c("-sS", "--max-time", as.character(timeout_sec), "-X", "POST", url,
           header_args, "-H", "Content-Type: application/json",
           "--data-binary", paste0("@", tmp_body), "-o", tmp_out, "-w", "%{http_code}")
  # shQuote() is not optional here: system2() joins `args` with plain spaces
  # and hands the result to a shell -- it does NOT quote elements itself, so
  # any element containing a space (every "-H" value: "name: value" always
  # has one, and OpenAI's "Authorization: Bearer <key>" has two) gets
  # word-split into extra shell tokens instead of staying one argument. A
  # live Gemini API call this session actually hit this: curl read the
  # split-off "application/json" tail of the Content-Type header as a second
  # positional URL and failed with "Could not resolve host: application"
  # while curl's own reported exit code (6) got silently swallowed into
  # $status alongside the real (buggy) "%{http_code}" text -- worth knowing
  # if a future error here ever looks like a status code with extra digits.
  status_lines <- system2("curl", args = shQuote(args), stdout = TRUE, stderr = TRUE)
  http_code <- suppressWarnings(as.integer(status_lines[length(status_lines)]))
  resp_text <- if (file.exists(tmp_out)) paste(readLines(tmp_out, warn = FALSE, encoding = "UTF-8"), collapse = "\n") else ""
  if (is.na(http_code)) {
    stop("curl request to ", url, " failed to complete: ", paste(status_lines, collapse = " "))
  }
  list(status = http_code, body = resp_text)
}

.default_api_key <- function(backend) {
  key <- switch(backend,
    anthropic = Sys.getenv("ANTHROPIC_API_KEY", ""),
    openai = Sys.getenv("OPENAI_API_KEY", ""),
    gemini = Sys.getenv("GEMINI_API_KEY", ""),
    stop("Unknown backend: ", backend))
  if (!nzchar(key)) {
    stop("No API key set for backend '", backend, "'. Set ",
        switch(backend, anthropic = "ANTHROPIC_API_KEY", openai = "OPENAI_API_KEY", gemini = "GEMINI_API_KEY"),
        " in the environment, or enter it in the AI Assistant panel's key field.")
  }
  key
}

.default_model <- function(backend) {
  switch(backend, anthropic = "claude-sonnet-5", openai = "gpt-4o", gemini = "gemini-2.0-flash",
        stop("Unknown backend: ", backend))
}

## ---- Anthropic Messages API adapter ------------------------------------------
# Canonical messages/tools are already in (almost exactly) Anthropic's own
# shape, so this is close to a passthrough plus the HTTP/auth wiring.
.llm_call_anthropic <- function(messages, tools, system_prompt, model, api_key,
                                max_tokens = 4096) {
  # unname() matters: tools (AGENT_TOOLS) is a NAMED list (named by tool
  # name for $-lookup elsewhere), and jsonlite::toJSON serializes a named
  # list as a JSON OBJECT rather than an ARRAY -- silently sending
  # '"tools":{"parse_sequence":{...}, ...}' instead of the array of tool
  # specs every provider's API actually expects. Caught via a real Gemini
  # API round-trip (HTTP 400: "Unknown name 'parse_sequence' at
  # 'tools[0].function_declarations'"), which exists for exactly this class
  # of bug -- an object where an array was expected reads, to jsonlite, as
  # "this key doesn't belong here". Same fix needed in every adapter below.
  tool_specs <- lapply(unname(tools), function(t) list(name = t$name, description = t$description,
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
  # unname() -- see the comment on the Anthropic adapter's own tool_specs
  # line above; the same named-list-serializes-as-object bug applies here.
  tool_specs <- lapply(unname(tools), function(t) list(type = "function", `function` = list(
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

## ---- Gemini generateContent adapter -------------------------------------------
# Gemini's tool-calling shape differs from both of the above: a functionCall
# part carries no id (unlike Anthropic's tool_use.id / OpenAI's tool_calls[].id),
# and a functionResponse part is matched back to it by *name*, not id. Our
# canonical tool_result block only carries tool_use_id, so converting TO
# Gemini needs an id->name lookup built from every tool_use block seen so far
# in the history; converting FROM Gemini just invents a local id (Gemini never
# needs it back) so the rest of run_agent_turn()'s bookkeeping works unchanged.
.canonical_to_gemini_contents <- function(messages) {
  id_to_name <- list()
  for (m in messages) {
    for (b in m$content) if (identical(b$type, "tool_use")) id_to_name[[b$id]] <- b$name
  }
  contents <- list()
  for (m in messages) {
    parts <- list()
    if (identical(m$role, "assistant")) {
      for (b in m$content) {
        if (identical(b$type, "text") && nzchar(b$text %||% "")) {
          parts[[length(parts) + 1]] <- list(text = b$text)
        } else if (identical(b$type, "tool_use")) {
          parts[[length(parts) + 1]] <- list(functionCall = list(name = b$name, args = b$input))
        }
      }
      if (length(parts) > 0) contents[[length(contents) + 1]] <- list(role = "model", parts = parts)
    } else {
      for (b in m$content) {
        if (identical(b$type, "text") && nzchar(b$text %||% "")) {
          parts[[length(parts) + 1]] <- list(text = b$text)
        } else if (identical(b$type, "tool_result")) {
          name <- id_to_name[[b$tool_use_id]] %||% "unknown_tool"
          resp <- tryCatch(jsonlite::fromJSON(b$content, simplifyVector = FALSE),
                           error = function(e) list(result = b$content))
          parts[[length(parts) + 1]] <- list(functionResponse = list(name = name, response = resp))
        }
      }
      if (length(parts) > 0) contents[[length(contents) + 1]] <- list(role = "user", parts = parts)
    }
  }
  contents
}

.gemini_parts_to_canonical <- function(parts) {
  blocks <- list()
  call_i <- 0L
  for (p in parts %||% list()) {
    if (!is.null(p$text) && nzchar(p$text)) {
      blocks[[length(blocks) + 1]] <- list(type = "text", text = p$text)
    } else if (!is.null(p$functionCall)) {
      call_i <- call_i + 1L
      blocks[[length(blocks) + 1]] <- list(type = "tool_use", id = paste0("gemini_call_", call_i),
                                           name = p$functionCall$name, input = p$functionCall$args %||% list())
    }
  }
  blocks
}

.llm_call_gemini <- function(messages, tools, system_prompt, model, api_key, max_tokens = 4096) {
  # unname() -- see the comment on the Anthropic adapter's own tool_specs
  # line; the same named-list-serializes-as-object bug is what this session's
  # own live Gemini smoke test actually caught (HTTP 400, "Unknown name
  # 'parse_sequence' at 'tools[0].function_declarations'").
  tool_specs <- lapply(unname(tools), function(t) list(name = t$name, description = t$description,
                                                        parameters = t$input_schema))
  body <- list(system_instruction = list(parts = list(list(text = system_prompt))),
              contents = .canonical_to_gemini_contents(messages),
              tools = list(list(functionDeclarations = tool_specs)),
              generationConfig = list(maxOutputTokens = max_tokens))
  url <- sprintf("https://generativelanguage.googleapis.com/v1beta/models/%s:generateContent?key=%s",
                model, utils::URLencode(api_key, reserved = TRUE))
  resp <- .http_post_json(url, headers = list(), body = body)
  parsed <- tryCatch(jsonlite::fromJSON(resp$body, simplifyVector = FALSE),
                     error = function(e) stop("Gemini API returned unparseable JSON (HTTP ",
                                              resp$status, "): ", substr(resp$body, 1, 500)))
  if (resp$status >= 400) {
    msg <- if (!is.null(parsed$error$message)) parsed$error$message else resp$body
    stop("Gemini API error (HTTP ", resp$status, "): ", msg)
  }
  cand <- parsed$candidates[[1]]
  content_blocks <- .gemini_parts_to_canonical(cand$content$parts)
  has_call <- any(vapply(content_blocks, function(b) identical(b$type, "tool_use"), logical(1)))
  list(role = "assistant", content = content_blocks,
      stop_reason = if (has_call) "tool_use" else "end_turn")
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
    gemini = function(msgs) .llm_call_gemini(msgs, AGENT_TOOLS, system_prompt,
                                             model %||% .default_model("gemini"),
                                             api_key %||% .default_api_key("gemini")),
    stub = local({
      i <- 0L
      function(msgs) { i <<- i + 1L; stub_responses[[i]] }
    }),
    stop("Unknown backend '", backend, "': use 'anthropic', 'openai', 'gemini', or 'stub'.")
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
