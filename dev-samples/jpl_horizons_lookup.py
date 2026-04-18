#!/usr/bin/env python3
"""
jpl_horizons_lookup.py
Version: 0.1.0

A small standalone command-line tool that queries the JPL Horizons Lookup
API and prints the aliases recognized by Horizons for a given object.

This script focuses on the lookup step only. It does not fetch ephemerides
or compute alt/az positions.

------------------------------------------------------------
OVERVIEW
------------------------------------------------------------

Given a query string, the script calls the official Horizons Lookup API and
prints either:

• a candidate list when the query matches multiple bodies
• a compact object summary when the query resolves to a single body
• a "not found" message when no object matches

The lookup API is useful for correlating names, designations, SPK IDs, and
other aliases before fetching ephemerides from Horizons.

------------------------------------------------------------
API USED
------------------------------------------------------------

Horizons Lookup API:
    https://ssd-api.jpl.nasa.gov/api/horizons_lookup.api

Official documentation:
    https://ssd-api.jpl.nasa.gov/doc/horizons_lookup.html

------------------------------------------------------------
DEFAULT BEHAVIOR
------------------------------------------------------------

By default, the script prints a human-readable summary. Optional JSON
output is available for scripting.

------------------------------------------------------------
OPTIONS
------------------------------------------------------------

--group GROUP
    Restrict the search to a Horizons object group.
    Allowed values: ast, com, pln, sat, sct, mb, sb

--json
    Emit the raw normalized payload as JSON.

--compact
    Reduce output formatting to a terse one-line style.

--version
    Show program version.

------------------------------------------------------------
EXAMPLES
------------------------------------------------------------

Lookup by name:

    python jpl_horizons_lookup.py Ceres

Lookup by designation:

    python jpl_horizons_lookup.py "2004 MN4"

Limit to small bodies:

    python jpl_horizons_lookup.py Juno --group sb

Return JSON:

    python jpl_horizons_lookup.py Ceres --json

------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from typing import Any


VERSION = "0.1.0"
BASE_URL = "https://ssd.jpl.nasa.gov/api/horizons_lookup.api"
DEFAULT_USER_AGENT = f"jpl_horizons_lookup/{VERSION}"
ALLOWED_GROUPS = {"ast", "com", "pln", "sat", "sct", "mb", "sb"}


def build_url(query: str, group: str | None) -> str:
    params: dict[str, str] = {"sstr": query}
    if group:
        params["group"] = group
    return BASE_URL + "?" + urllib.parse.urlencode(params)


def fetch_payload(url: str, user_agent: str) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"User-Agent": user_agent})
    with urllib.request.urlopen(req, timeout=20) as resp:
        charset = resp.headers.get_content_charset() or "utf-8"
        body = resp.read().decode(charset)
    payload = json.loads(body)
    if not isinstance(payload, dict):
        raise RuntimeError("Horizons lookup returned an unexpected payload type")
    return payload


def format_candidate(entry: dict[str, Any]) -> str:
    name = str(entry.get("name", "")).strip()
    pdes = str(entry.get("pdes", "")).strip()
    spkid = str(entry.get("spkid", "")).strip()
    type_name = str(entry.get("type", "")).strip()
    parts = [part for part in (name, pdes, f"spkid={spkid}" if spkid else "", type_name) if part]
    return " | ".join(parts)


def format_object_summary(payload: dict[str, Any]) -> list[str]:
    result = payload.get("result") or []
    if not isinstance(result, list) or not result:
        return []
    entry = result[0]
    if not isinstance(entry, dict):
        return []

    lines = []
    name = str(entry.get("name", "")).strip()
    if name:
        lines.append(name)

    fields = [
        ("type", entry.get("type")),
        ("designation", entry.get("pdes")),
        ("spkid", entry.get("spkid")),
    ]
    for label, value in fields:
        text = str(value).strip() if value is not None else ""
        if text and text.lower() != "none":
            lines.append(f"{label}: {text}")

    aliases = entry.get("alias") or []
    if isinstance(aliases, list) and aliases:
        alias_text = ", ".join(str(item).strip() for item in aliases if str(item).strip())
        if alias_text:
            lines.append(f"aliases: {alias_text}")
    return lines


def print_payload(payload: dict[str, Any], *, compact: bool) -> None:
    count = payload.get("count")
    result = payload.get("result")

    if count == 0 or count == "0":
        print("No matching object found.")
        return

    if isinstance(result, list) and len(result) == 1:
        if compact:
            print(format_candidate(result[0]) if isinstance(result[0], dict) else json.dumps(payload, ensure_ascii=True, sort_keys=True))
            return
        for line in format_object_summary(payload):
            print(line)
        return

    if isinstance(result, list) and len(result) > 1:
        if compact:
            for entry in result:
                if isinstance(entry, dict):
                    print(format_candidate(entry))
            return
        print(f"{count} match(es):")
        for entry in result:
            if isinstance(entry, dict):
                print(f"- {format_candidate(entry)}")
        return

    if compact:
        print(json.dumps(payload, ensure_ascii=True, sort_keys=True))
        return

    print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="jpl_horizons_lookup.py",
        description="Query the JPL Horizons Lookup API for object aliases.",
    )
    parser.add_argument("query", help="Search string passed to the Horizons sstr parameter.")
    parser.add_argument(
        "--group",
        choices=sorted(ALLOWED_GROUPS),
        help="Restrict the search to a Horizons object group.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit the raw Horizons lookup payload as JSON.",
    )
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Use terse one-line output for candidates and single-object results.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    query = str(args.query).strip()
    if not query:
        print("Query must not be empty.", file=sys.stderr)
        return 2

    url = build_url(query, str(args.group).strip() if args.group else None)
    try:
        payload = fetch_payload(url, DEFAULT_USER_AGENT)
    except urllib.error.URLError as exc:
        print(f"Horizons lookup request failed: {exc}", file=sys.stderr)
        return 1
    except json.JSONDecodeError as exc:
        print(f"Horizons lookup response was not valid JSON: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        print(f"Horizons lookup failed: {exc}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))
    else:
        print_payload(payload, compact=bool(args.compact))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
