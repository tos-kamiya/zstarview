#!/usr/bin/env python3
"""
jpl_sbdb_search.py
Version: 0.1.0

A small standalone command-line tool that queries the JPL Small-Body
Database (SBDB) search API and prints matching small-body candidates.

This script focuses on the search step only. It does not fetch observer
ephemerides or compute alt/az positions.

------------------------------------------------------------
OVERVIEW
------------------------------------------------------------

Given a query string, the script calls the official SBDB API and prints
either:

• a candidate list when the query is ambiguous
• a compact object summary when the query resolves to a single body
• a "not found" message when no object matches

The SBDB API accepts names, designations, MPC packed designations, and
numeric identifiers. It is suitable for finding asteroids, comets, and
other small bodies including dwarf planets.

------------------------------------------------------------
API USED
------------------------------------------------------------

SBDB API:
    https://ssd-api.jpl.nasa.gov/sbdb.api

Official documentation:
    https://ssd-api.jpl.nasa.gov/doc/sbdb.html

------------------------------------------------------------
DEFAULT BEHAVIOR
------------------------------------------------------------

By default, the script prints a human-readable summary. Optional JSON
output is available for scripting.

------------------------------------------------------------
OPTIONS
------------------------------------------------------------

--json
    Emit the raw normalized payload as JSON.

--compact
    Reduce output formatting to a terse one-line style.

--version
    Show program version.

------------------------------------------------------------
EXAMPLES
------------------------------------------------------------

Search by name:

    python jpl_sbdb_search.py Ceres

Search by partial designation:

    python jpl_sbdb_search.py 141P

Search with JSON output:

    python jpl_sbdb_search.py Eris --json

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
BASE_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"
DEFAULT_USER_AGENT = f"jpl_sbdb_search/{VERSION}"


def build_url(query: str) -> str:
    params = {"sstr": query}
    return BASE_URL + "?" + urllib.parse.urlencode(params)


def fetch_payload(url: str, user_agent: str) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"User-Agent": user_agent})
    with urllib.request.urlopen(req, timeout=20) as resp:
        charset = resp.headers.get_content_charset() or "utf-8"
        body = resp.read().decode(charset)
    payload = json.loads(body)
    if not isinstance(payload, dict):
        raise RuntimeError("SBDB API returned an unexpected payload type")
    return payload


def format_candidate(entry: dict[str, Any]) -> str:
    pdes = str(entry.get("pdes", "")).strip()
    name = str(entry.get("name", "")).strip()
    if pdes and name:
        return f"{pdes}  {name}"
    if pdes:
        return pdes
    if name:
        return name
    return json.dumps(entry, ensure_ascii=True, sort_keys=True)


def format_object_summary(payload: dict[str, Any]) -> list[str]:
    obj = payload.get("object") or {}
    orbit = payload.get("orbit") or {}
    if not isinstance(obj, dict):
        obj = {}
    if not isinstance(orbit, dict):
        orbit = {}

    lines = []
    shortname = str(obj.get("shortname", "")).strip()
    fullname = str(obj.get("fullname", "")).strip()
    if shortname:
        lines.append(shortname)
    elif fullname:
        lines.append(fullname)

    fields = [
        ("designation", obj.get("des")),
        ("name", obj.get("name")),
        ("kind", obj.get("kind")),
        ("prefix", obj.get("prefix")),
        ("orbit_class", (obj.get("orbit_class") or {}).get("name") if isinstance(obj.get("orbit_class"), dict) else None),
        ("spkid", obj.get("spkid")),
        ("orbit_id", obj.get("orbit_id")),
        ("soln_date", orbit.get("soln_date")),
    ]
    for label, value in fields:
        text = str(value).strip() if value is not None else ""
        if text and text.lower() != "none":
            lines.append(f"{label}: {text}")
    return lines


def print_payload(payload: dict[str, Any], *, compact: bool) -> None:
    code = payload.get("code")
    message = str(payload.get("message", "")).strip()

    if code == 200 and "object" in payload:
        if compact:
            obj = payload.get("object") or {}
            if isinstance(obj, dict):
                print(format_candidate({
                    "pdes": obj.get("des", ""),
                    "name": obj.get("name", ""),
                }))
            else:
                print(json.dumps(payload, ensure_ascii=True, sort_keys=True))
            return
        for line in format_object_summary(payload):
            print(line)
        return

    if code == 300 and "list" in payload:
        candidates = payload.get("list") or []
        if not isinstance(candidates, list):
            raise RuntimeError("SBDB API returned an invalid candidate list")
        if compact:
            for entry in candidates:
                if isinstance(entry, dict):
                    print(format_candidate(entry))
            return
        print(message or "Multiple objects matched:")
        for entry in candidates:
            if isinstance(entry, dict):
                print(f"- {format_candidate(entry)}")
        return

    if code == 404:
        print(message or "No matching object found.")
        return

    if compact:
        print(json.dumps(payload, ensure_ascii=True, sort_keys=True))
        return

    if message:
        print(message)
    print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="jpl_sbdb_search.py",
        description="Query the JPL SBDB search API for small-body candidates.",
    )
    parser.add_argument("query", help="Search string passed to the SBDB sstr parameter.")
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit the raw SBDB payload as JSON.",
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

    url = build_url(query)
    try:
        payload = fetch_payload(url, DEFAULT_USER_AGENT)
    except urllib.error.URLError as exc:
        print(f"SBDB request failed: {exc}", file=sys.stderr)
        return 1
    except json.JSONDecodeError as exc:
        print(f"SBDB response was not valid JSON: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        print(f"SBDB search failed: {exc}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))
    else:
        print_payload(payload, compact=bool(args.compact))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
