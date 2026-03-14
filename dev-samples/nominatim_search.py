#!/usr/bin/env python3
"""
nominatim_search.py
Version: 0.2.0

A small standalone command-line tool that queries the OpenStreetMap
Nominatim geocoding service and returns candidate locations with
latitude and longitude.

This script intentionally uses only the Python standard library so it
can be distributed as a single file without dependencies.

------------------------------------------------------------
OVERVIEW
------------------------------------------------------------

Given a place name or free-form location string, the script queries
the Nominatim Search API and returns matching locations.

Results are normalized, sorted by Nominatim's "importance" score
in descending order, and then optionally sliced using a simple
Python-like slice syntax subset.

By default, the script behaves as if:

    --slice 0

were specified, meaning that it returns only the single most
important candidate.

Alternative output modes allow returning sliced results as JSON
or JSONL, or printing only the number of matched candidates.

------------------------------------------------------------
FEATURES
------------------------------------------------------------

• Uses the official Nominatim search endpoint
• Returns normalized output fields
• Importance-based sorting
• Python-like slice selection
• Optional full result output (JSON / JSONL)
• Candidate count output
• No external dependencies
• Suitable for embedding in scripts or tools

------------------------------------------------------------
DEFAULT BEHAVIOR
------------------------------------------------------------

By default:

• Results are sorted by importance, descending
• The default slice is "0"
• A single JSON object is returned
• If no location is found, the program exits with an error

------------------------------------------------------------
SLICE SYNTAX
------------------------------------------------------------

The --slice option accepts a small subset of Python slicing syntax.

Examples:

    0       first result only
    3       fourth result only
    :       all results
    :3      first three results
    1:4     second through fourth results
    2:      third and later results

Notes:

• Step values are not supported
• Negative indices are not supported
• The default is "0"

------------------------------------------------------------
OPTIONS
------------------------------------------------------------

-n, --limit N
    Maximum number of candidates requested from Nominatim
    (range 1–40, default: 5)

-c, --countrycode CODE
    Restrict search to a specific country (ISO 3166-1 alpha-2)

    Example:
        -c jp

-l, --lang LANG
    Language preference via the HTTP Accept-Language header

    Examples:
        -l ja
        -l en
        -l ja,en

-s, --slice SPEC
    Slice results after sorting by importance.
    Default: 0

    Examples:
        -s 0
        -s :
        -s :5
        -s 1:4

--json
    Return sliced results as a JSON array.
    If no candidates are found, return [].

--jsonl
    Return sliced results as JSON Lines.
    If no candidates are found, print nothing.

--count
    Print the number of normalized candidates after sorting.
    This ignores --slice and prints only the total count.

--version
    Show program version.

------------------------------------------------------------
EXAMPLES
------------------------------------------------------------

Default (best match only):

    python nominatim_search.py "Matsue Castle"

Return all candidates:

    python nominatim_search.py "Matsue" -s : --json

Return the top three candidates:

    python nominatim_search.py "Matsue" -s :3 --json

Return only the second candidate:

    python nominatim_search.py "Matsue" -s 1

Restrict to Japan:

    python nominatim_search.py "Matsue" -c jp -s : --json

Print only the number of candidates:

    python nominatim_search.py "Matsue" --count

------------------------------------------------------------
OUTPUT FORMAT
------------------------------------------------------------

Each normalized result contains:

    name        Full display name from Nominatim
    lat         Latitude (float)
    lon         Longitude (float)
    category    OSM category (e.g. place, tourism, historic)
    type        OSM type (e.g. city, castle, station)
    importance  Nominatim importance score

------------------------------------------------------------
NOMINATIM USAGE POLICY
------------------------------------------------------------

This script uses the public Nominatim API:

    https://nominatim.openstreetmap.org

Users should respect the usage policy:

• Send a valid User-Agent
• Avoid high-frequency automated queries

See:
https://operations.osmfoundation.org/policies/nominatim/

------------------------------------------------------------
LICENSE
------------------------------------------------------------

This script is provided as-is for convenience. Users should verify
compliance with the Nominatim usage policy when integrating it into
their own applications.

------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.parse
import urllib.request


VERSION = "0.2.0"

BASE_URL = "https://nominatim.openstreetmap.org/search"
DEFAULT_USER_AGENT = f"nominatim_search/{VERSION}"


def build_url(query: str, limit: int, countrycode: str | None) -> str:
    params = {
        "q": query,
        "format": "jsonv2",
        "limit": str(limit),
        "addressdetails": "0",
        "dedupe": "1",
    }

    if countrycode:
        params["countrycodes"] = countrycode

    return BASE_URL + "?" + urllib.parse.urlencode(params)


def fetch(url: str, language: str, user_agent: str) -> list[dict]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": user_agent,
            "Accept-Language": language,
        },
    )

    with urllib.request.urlopen(req, timeout=15) as resp:
        charset = resp.headers.get_content_charset() or "utf-8"
        body = resp.read().decode(charset)

    data = json.loads(body)

    if not isinstance(data, list):
        raise ValueError("unexpected response format")

    return data


def normalize(items: list[dict]) -> list[dict]:
    results = []

    for item in items:
        try:
            lat = float(item["lat"])
            lon = float(item["lon"])
        except (KeyError, TypeError, ValueError):
            continue

        importance = item.get("importance") or 0.0

        results.append(
            {
                "name": item.get("display_name"),
                "lat": lat,
                "lon": lon,
                "category": item.get("category") or item.get("class"),
                "type": item.get("type"),
                "importance": float(importance),
            }
        )

    results.sort(key=lambda r: r["importance"], reverse=True)
    return results


def parse_slice_spec(spec: str | None) -> slice:
    if spec is None or spec == "":
        return slice(0, 1)

    if ":" not in spec:
        index = int(spec)
        if index < 0:
            raise ValueError("negative indices are not supported")
        return slice(index, index + 1)

    parts = spec.split(":")
    if len(parts) != 2:
        raise ValueError("slice step is not supported")

    start_text, stop_text = parts

    start = None
    stop = None

    if start_text != "":
        start = int(start_text)
        if start < 0:
            raise ValueError("negative indices are not supported")

    if stop_text != "":
        stop = int(stop_text)
        if stop < 0:
            raise ValueError("negative indices are not supported")

    return slice(start, stop)


def apply_slice(results: list[dict], spec: str | None) -> list[dict]:
    sl = parse_slice_spec(spec)
    return results[sl]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Search place names using Nominatim."
    )

    parser.add_argument("query", help="place name or free-form query")

    parser.add_argument(
        "-n",
        "--limit",
        type=int,
        default=5,
        help="number of candidates requested from Nominatim (1–40)",
    )

    parser.add_argument(
        "-c",
        "--countrycode",
        help="restrict search to country code (e.g. jp)",
    )

    parser.add_argument(
        "-l",
        "--lang",
        default="ja",
        help="Accept-Language header (default: ja)",
    )

    parser.add_argument(
        "-s",
        "--slice",
        dest="slice_spec",
        default="0",
        help='result slice after sorting, e.g. "0", ":", ":3", "1:4" (default: 0)',
    )

    parser.add_argument(
        "--user-agent",
        default=DEFAULT_USER_AGENT,
        help="User-Agent header",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}",
    )

    parser.add_argument(
        "--count",
        action="store_true",
        help="print the number of normalized candidates",
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--json",
        action="store_true",
        help="output sliced results as a JSON array",
    )
    group.add_argument(
        "--jsonl",
        action="store_true",
        help="output sliced results as JSONL",
    )

    args = parser.parse_args()

    if not (1 <= args.limit <= 40):
        print("error: limit must be between 1 and 40", file=sys.stderr)
        return 2

    try:
        parse_slice_spec(args.slice_spec)
    except ValueError as e:
        print(f"error: invalid slice: {e}", file=sys.stderr)
        return 2
    except Exception:
        print("error: invalid slice syntax", file=sys.stderr)
        return 2

    url = build_url(args.query, args.limit, args.countrycode)

    try:
        raw = fetch(url, args.lang, args.user_agent)
        results = normalize(raw)
    except urllib.error.HTTPError as e:
        print(f"HTTP error: {e.code} {e.reason}", file=sys.stderr)
        return 1
    except urllib.error.URLError as e:
        print(f"network error: {e.reason}", file=sys.stderr)
        return 1
    except json.JSONDecodeError as e:
        print(f"JSON parse error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    if args.count:
        print(len(results))
        return 0

    selected = apply_slice(results, args.slice_spec)

    if args.json:
        json.dump(selected, sys.stdout, ensure_ascii=False, indent=2)
        print()
        return 0

    if args.jsonl:
        for item in selected:
            print(json.dumps(item, ensure_ascii=False))
        return 0

    if not selected:
        print(f"error: no result for '{args.query}'", file=sys.stderr)
        return 1

    if len(selected) == 1:
        json.dump(selected[0], sys.stdout, ensure_ascii=False, indent=2)
    else:
        json.dump(selected, sys.stdout, ensure_ascii=False, indent=2)
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
