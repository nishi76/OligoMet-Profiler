# OligoMet Profiler MCP server

Exposes the tool registry in `R/agent_tools.R` (parse a sequence, build the
metabolite library, match MS1 features, confirm MS2, summarize degradation,
compare groups, etc.) to any MCP client -- Claude Code, Claude Desktop,
claude.ai -- as ordinary MCP tools. Part of the plan at
`/root/.claude/plans/staged-mapping-blanket.md` ("Track B").

No analytical logic lives here. This is a thin, two-process shim:

```
MCP client (Claude Code / Desktop / claude.ai)
        | stdio, MCP protocol
        v
inst/mcp_server/server.py   (this directory -- protocol translation only)
        | HTTP/JSON
        v
R/agent_api.R  (via run_agent_api.R at the repo root -- the real R analysis code)
```

## Setup

1. **Start the R API** (from the repository root, in its own terminal/process):
   ```
   Rscript run_agent_api.R 8787
   ```
   Leave it running. It holds every analysis module loaded and warm, so a
   tool call doesn't pay R's startup cost each time.

2. **Install the Python side's one dependency:**
   ```
   pip install -r inst/mcp_server/requirements.txt
   ```
   If `import mcp` fails with a `cryptography`/`_cffi_backend` error, that's
   an outdated `cffi` on the system, not this package -- `pip install
   --upgrade cffi` fixes it (hit and fixed this exact way while building it).

3. **Point an MCP client at `server.py`.** For Claude Code, add to your MCP
   config (`claude mcp add`, or the client's config file):
   ```json
   {
     "mcpServers": {
       "oligomet-profiler": {
         "command": "python3",
         "args": ["/absolute/path/to/inst/mcp_server/server.py"],
         "env": { "OLIGOMET_AGENT_API": "http://127.0.0.1:8787" }
       }
     }
   }
   ```
   `OLIGOMET_AGENT_API` only needs setting if the R API isn't on the default
   `http://127.0.0.1:8787`.

## Stateless by design

Every call is independent -- there's no session between one tool call and
the next the way the in-app Shiny chat assistant (Track A) has. Pass
whatever a tool needs directly in its arguments (a `sequence` string, a
previously-returned `spec`/`metabolites` object, etc.) rather than assuming
an earlier call's result carries forward automatically. Each tool's
description in `R/agent_tools.R` documents what it accepts.
