# =============================================================================
# run_agent_api.R -- starts the persistent tool API (R/agent_api.R) that the
# MCP server (inst/mcp_server/) talks to.
#
# Run with:
#   Rscript run_agent_api.R           # port 8787
#   Rscript run_agent_api.R 9000      # or a specific port
#   AGENT_API_PORT=9000 Rscript run_agent_api.R
#
# See /root/.claude/plans/staged-mapping-blanket.md ("Track B") for how this
# fits together: this process holds the R analysis code warm; the Python MCP
# server is a thin protocol shim in front of it, so a real deployment only
# needs to keep this one R process running (e.g. under supervisord/systemd)
# rather than a fresh Rscript per tool call.
# =============================================================================

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else "."
}, error = function(e) ".")

if (!file.exists(file.path(script_dir, "R", "chemistry_dict.R"))) {
  stop("Cannot locate pipeline modules (chemistry_dict.R not found). ",
       "Run this script from the package root, e.g. Rscript run_agent_api.R.")
}

for (f in c("about.R", "progress_utils.R", "chemistry_dict.R", "oligo_io.R",
           "metabolites.R", "mass_isotope.R", "fragments.R", "ms_matching.R",
           "degradation.R", "statistics.R", "agent_tools.R", "agent_api.R")) {
  source(file.path(script_dir, "R", f))
}

port <- as.integer(Sys.getenv("AGENT_API_PORT", "8787"))
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) > 0) port <- as.integer(cli_args[1])
if (is.na(port)) stop("Invalid port.")

run_agent_api(host = "127.0.0.1", port = port)
