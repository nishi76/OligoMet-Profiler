# =============================================================================
# agent_api.R
# A tiny persistent HTTP JSON API exposing R/agent_tools.R's tool registry,
# so a caller in another process/language (the MCP server in
# inst/mcp_server/) doesn't need to shell out to a fresh Rscript per call.
# Part of the plan at /root/.claude/plans/staged-mapping-blanket.md -- Track B.
#
# Deviation from the plan's sketch, flagged explicitly: the plan mentioned
# plumber for this. plumber isn't installable in this development sandbox
# (no CRAN access, not packaged for apt here). httpuv IS already installed
# (Shiny depends on it), so this is written directly against httpuv's own
# minimal app interface instead -- one persistent process, the same
# warm-process latency win the plan (and the MSAgent paper) call for, no new
# dependency beyond what Shiny already pulls in. Swapping in a plumber-based
# router later, if a deployment has that package, would only mean replacing
# this one file -- inst/mcp_server/ talks to it purely over HTTP/JSON and
# wouldn't need to change.
#
# Endpoints:
#   GET  /health  -> "ok"
#   GET  /tools   -> JSON array of {name, description, input_schema}
#   POST /invoke  -> body {"tool": "...", "args": {...}, "ctx": {...}}
#                    -> the tool's JSON result, or {"error": "..."} (HTTP 400)
#
# ctx is normally {} here -- the MCP server is a STATELESS caller per the
# plan (each call independent; no Shiny-style reactiveValues to hold
# session state between calls), so it supplies everything a tool needs
# through args instead. Every tool in the registry already supports this --
# see the ctx-vs-args fallback in each .resolve_*() helper in agent_tools.R.
# =============================================================================

.agent_api_json_response <- function(status, obj) {
  list(status = status,
      headers = list("Content-Type" = "application/json; charset=utf-8"),
      body = as.character(jsonlite::toJSON(obj, auto_unbox = TRUE, na = "null", digits = NA)))
}

.agent_api_read_body <- function(req) {
  raw <- req$rook.input$read()
  if (is.null(raw) || length(raw) == 0) return(list())
  txt <- rawToChar(raw)
  if (!nzchar(trimws(txt))) return(list())
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

.agent_api_app <- list(
  call = function(req) {
    path <- req$PATH_INFO
    method <- req$REQUEST_METHOD
    tryCatch({
      if (identical(method, "GET") && identical(path, "/health")) {
        return(list(status = 200L, headers = list("Content-Type" = "text/plain"), body = "ok"))
      }
      if (identical(method, "GET") && identical(path, "/tools")) {
        return(.agent_api_json_response(200L, unname(agent_tool_specs())))
      }
      if (identical(method, "POST") && identical(path, "/invoke")) {
        payload <- .agent_api_read_body(req)
        if (is.null(payload$tool)) return(.agent_api_json_response(400L, list(error = "Missing 'tool' field")))
        result <- call_agent_tool(payload$tool, payload$args %||% list(), payload$ctx %||% list())
        return(.agent_api_json_response(200L, result))
      }
      .agent_api_json_response(404L, list(error = paste("No such route:", method, path)))
    }, error = function(e) .agent_api_json_response(400L, list(error = conditionMessage(e))))
  }
)

# Starts the API and blocks the calling process (httpuv::runServer() never
# returns) -- meant to run as its own process, see run_agent_api.R at the
# repository root.
run_agent_api <- function(host = "127.0.0.1", port = 8787) {
  message("OligoMet agent API listening on http://", host, ":", port,
          " (", length(AGENT_TOOLS), " tools registered)")
  httpuv::runServer(host, port, .agent_api_app)
}

# Non-blocking variant, for tests/interactive use -- returns the server
# handle; stop it with httpuv::stopServer(handle) or httpuv::stopAllServers().
start_agent_api_background <- function(host = "127.0.0.1", port = 8787) {
  httpuv::startServer(host, port, .agent_api_app)
}
