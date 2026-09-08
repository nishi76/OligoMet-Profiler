"""
server.py -- MCP server for OligoMet Profiler (Track B of the plan at
/root/.claude/plans/staged-mapping-blanket.md).

This is a thin protocol shim: it lists tools and dispatches tool calls by
talking HTTP/JSON to the persistent R API in R/agent_api.R (started
separately via run_agent_api.R at the repo root). NO analytical logic lives
here -- the tool registry and every wrapped function stay in R, unchanged.
Every tool call is stateless (ctx = {}): the R side's tools already fall
back to reconstructing whatever they need from the call's own arguments
when there's no session context to draw on (see the .resolve_*() helpers in
R/agent_tools.R), which is exactly the calling convention this front end
needs.

Run standalone (stdio transport, the standard for Claude Code/Desktop):
    OLIGOMET_AGENT_API=http://127.0.0.1:8787 python3 server.py

Prerequisite: the R API must already be running --
    Rscript run_agent_api.R 8787
(see inst/mcp_server/README.md for a from-scratch setup, including the one
environment quirk hit while building this: an outdated `cffi` breaking the
`mcp` package's `cryptography` import on some systems --
`pip install --upgrade cffi` fixes it.)
"""

import asyncio
import json
import os
import urllib.error
import urllib.request

from mcp.server import Server
from mcp.server.stdio import stdio_server
import mcp.types as types

API_BASE = os.environ.get("OLIGOMET_AGENT_API", "http://127.0.0.1:8787").rstrip("/")


def _http_get(path: str, timeout: float = 30.0):
    with urllib.request.urlopen(API_BASE + path, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _http_post_json(path: str, payload: dict, timeout: float = 120.0):
    """POST JSON to the R agent API.

    A response the R API actually sent back -- 200 with a tool's result, or
    400 with {"error": "..."} for an invalid call (unknown tool, missing
    argument, no metabolite library loaded, etc.) -- is returned as a plain
    dict either way: those are informative, expected outcomes the model
    should read and act on, not a failed call. A transport failure (the R
    API isn't reachable at all, DNS/connection refused, timeout) is a
    genuinely failed call, so that raises instead -- the framework then
    reports it to the MCP client as isError=True, rather than looking like
    an ordinary tool result.
    """
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        API_BASE + path, data=data,
        headers={"Content-Type": "application/json"}, method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8")
        try:
            return json.loads(body)
        except json.JSONDecodeError:
            return {"error": body or f"HTTP {e.code} from the R agent API"}
    except urllib.error.URLError as e:
        raise RuntimeError(
            f"Could not reach the R agent API at {API_BASE}: {e.reason}. "
            "Is `Rscript run_agent_api.R` running?"
        ) from e


server = Server("oligomet-profiler")


@server.list_tools()
async def list_tools() -> list[types.Tool]:
    specs = _http_get("/tools")
    return [
        types.Tool(name=s["name"], description=s["description"], inputSchema=s["input_schema"])
        for s in specs
    ]


@server.call_tool()
async def call_tool(name: str, arguments: dict) -> list[types.TextContent]:
    result = _http_post_json("/invoke", {"tool": name, "args": arguments or {}, "ctx": {}})
    return [types.TextContent(type="text", text=json.dumps(result, ensure_ascii=False))]


async def main() -> None:
    async with stdio_server() as (read_stream, write_stream):
        await server.run(read_stream, write_stream, server.create_initialization_options())


if __name__ == "__main__":
    asyncio.run(main())
