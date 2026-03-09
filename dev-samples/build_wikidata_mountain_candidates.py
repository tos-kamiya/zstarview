#!/usr/bin/env python3
"""Collapse raw WDQS mountain query results into reviewable candidates.

This is the pre-normalization step for mountain viewpoints.

Input: raw JSON result saved from a manual WDQS query, typically the
`highest point of each country` query.

Output: a review-oriented JSON file keyed conceptually by `(country, item)`.
It keeps duplicate coordinates/elevations visible so a later curated seed can
be created deliberately, rather than silently flattening data too early.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


QID_RE = re.compile(r"/(Q\d+)$")
POINT_RE = re.compile(r"^Point\(([-0-9.]+) ([-0-9.]+)\)$")
PREFERRED_LABEL_LANGS = ("en", "ja")
WIKIDATA_API_URL = "https://www.wikidata.org/w/api.php"
USER_AGENT = "zstarview-dev-samples/1.0 (mountain candidate builder)"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collapse raw Wikidata mountain query rows into reviewable candidates."
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Raw JSON result path saved from WDQS.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Review-oriented candidate JSON output path.",
    )
    parser.add_argument(
        "--skip-wikidata-label-fetch",
        action="store_true",
        help="Do not fetch en/ja labels from Wikidata API.",
    )
    return parser.parse_args()


def parse_qid(entity_url: str) -> str:
    match = QID_RE.search((entity_url or "").strip())
    if not match:
        raise ValueError(f"Could not parse QID from {entity_url!r}")
    return match.group(1)


def parse_point(point_text: str) -> tuple[float, float]:
    match = POINT_RE.match((point_text or "").strip())
    if not match:
        raise ValueError(f"Could not parse coordinate point {point_text!r}")
    lon_text, lat_text = match.groups()
    return float(lat_text), float(lon_text)


@dataclass
class CandidateRecord:
    country_qid: str
    country_label: str
    item_qid: str
    item_label: str
    row_count: int = 0
    coordinates: set[tuple[float, float]] = field(default_factory=set)
    elevations: set[float] = field(default_factory=set)
    item_labels: dict[str, str] = field(default_factory=dict)
    country_labels: dict[str, str] = field(default_factory=dict)

    def add_row(self, row: dict[str, Any]) -> None:
        self.row_count += 1
        self.coordinates.add(parse_point(str(row["coord"])))
        self.elevations.add(float(row["elevation"]))

    def as_dict(self, shared_with_country_qids: list[str]) -> dict[str, object]:
        sorted_coords = sorted(self.coordinates)
        sorted_elevations = sorted(self.elevations)
        representative_lat, representative_lon = sorted_coords[0]
        max_elevation_m = max(sorted_elevations) if sorted_elevations else None
        return {
            "country_qid": self.country_qid,
            "country_label": self.country_label,
            "country_labels": dict(sorted(self.country_labels.items())),
            "item_qid": self.item_qid,
            "item_label": self.item_label,
            "item_labels": dict(sorted(self.item_labels.items())),
            "row_count": self.row_count,
            "coordinates": [
                {"latitude_deg": lat, "longitude_deg": lon}
                for lat, lon in sorted_coords
            ],
            "elevations_m": sorted_elevations,
            "latitude_deg": representative_lat,
            "longitude_deg": representative_lon,
            "location_arg": f"{representative_lat:.6f};{representative_lon:.6f}",
            "max_elevation_m": max_elevation_m,
            "coord_count": len(sorted_coords),
            "elevation_count": len(sorted_elevations),
            "needs_coord_review": len(sorted_coords) > 1,
            "needs_elevation_review": len(sorted_elevations) > 1,
            "is_shared_peak": bool(shared_with_country_qids),
            "shared_with_country_qids": shared_with_country_qids,
            "wikidata_url": f"https://www.wikidata.org/wiki/{self.item_qid}",
        }


def load_raw_rows(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict):
        rows = payload.get("results", {}).get("bindings", [])
    else:
        raise ValueError("Raw payload must be a list or a WDQS JSON results object.")
    if not isinstance(rows, list):
        raise ValueError("Raw payload does not contain a row list.")
    return [row for row in rows if isinstance(row, dict)]


def normalize_rows(rows: list[dict[str, Any]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], CandidateRecord] = {}
    countries_by_item: dict[str, set[str]] = {}

    for row in rows:
        country_qid = parse_qid(str(row["country"]))
        item_qid = parse_qid(str(row["item"]))
        key = (country_qid, item_qid)
        record = grouped.get(key)
        if record is None:
            record = CandidateRecord(
                country_qid=country_qid,
                country_label=str(row["countryLabel"]).strip(),
                item_qid=item_qid,
                item_label=str(row["itemLabel"]).strip(),
            )
            grouped[key] = record
        record.add_row(row)
        countries_by_item.setdefault(item_qid, set()).add(country_qid)

    normalized: list[dict[str, object]] = []
    for (country_qid, item_qid), record in grouped.items():
        shared_with = sorted(countries_by_item.get(item_qid, set()) - {country_qid})
        normalized.append(record.as_dict(shared_with))

    return sorted(
        normalized,
        key=lambda item: (
            str(item["country_label"]),
            -float(item["max_elevation_m"]) if item["max_elevation_m"] is not None else 0.0,
            str(item["item_label"]),
        ),
    )


def fetch_entity_labels(
    qids: list[str],
    languages: tuple[str, ...] = PREFERRED_LABEL_LANGS,
    batch_size: int = 50,
) -> dict[str, dict[str, str]]:
    labels_by_qid: dict[str, dict[str, str]] = {}
    language_param = "|".join(languages)

    for start_index in range(0, len(qids), batch_size):
        batch = qids[start_index : start_index + batch_size]
        params = urlencode(
            {
                "action": "wbgetentities",
                "ids": "|".join(batch),
                "props": "labels",
                "languages": language_param,
                "format": "json",
            }
        )
        request = Request(
            f"{WIKIDATA_API_URL}?{params}",
            headers={"User-Agent": USER_AGENT},
        )
        with urlopen(request, timeout=30) as response:
            payload = json.load(response)

        entities = payload.get("entities", {})
        if not isinstance(entities, dict):
            continue
        for qid, entity in entities.items():
            if not isinstance(entity, dict):
                continue
            raw_labels = entity.get("labels", {})
            if not isinstance(raw_labels, dict):
                continue
            labels: dict[str, str] = {}
            for lang in languages:
                lang_value = raw_labels.get(lang)
                if isinstance(lang_value, dict):
                    value = lang_value.get("value")
                    if isinstance(value, str) and value.strip():
                        labels[lang] = value.strip()
            if labels:
                labels_by_qid[qid] = labels

    return labels_by_qid


def merge_entity_labels(
    items: list[dict[str, object]],
    *,
    item_labels_by_qid: dict[str, dict[str, str]],
    country_labels_by_qid: dict[str, dict[str, str]],
) -> list[dict[str, object]]:
    enriched: list[dict[str, object]] = []
    for item in items:
        updated = dict(item)
        updated["item_labels"] = dict(sorted(item_labels_by_qid.get(str(item["item_qid"]), {}).items()))
        updated["country_labels"] = dict(sorted(country_labels_by_qid.get(str(item["country_qid"]), {}).items()))
        enriched.append(updated)
    return enriched


def main() -> None:
    args = parse_args()
    raw_payload = json.loads(args.input.read_text(encoding="utf-8"))
    rows = load_raw_rows(raw_payload)
    items = normalize_rows(rows)
    label_fetch_status = "skipped"
    if not args.skip_wikidata_label_fetch and items:
        item_qids = sorted({str(item["item_qid"]) for item in items})
        country_qids = sorted({str(item["country_qid"]) for item in items})
        try:
            item_labels_by_qid = fetch_entity_labels(item_qids)
            country_labels_by_qid = fetch_entity_labels(country_qids)
        except URLError as exc:
            print(f"Warning: failed to fetch Wikidata labels: {exc}")
            label_fetch_status = "failed"
        else:
            items = merge_entity_labels(
                items,
                item_labels_by_qid=item_labels_by_qid,
                country_labels_by_qid=country_labels_by_qid,
            )
            label_fetch_status = "fetched"
    output = {
        "source": "wikidata-highest-point-query",
        "source_query_result": str(args.input),
        "normalization": {
            "group_by": ["country_qid", "item_qid"],
            "representative_coord_rule": "first sorted distinct coordinate",
            "elevation_rule": "keep all distinct elevations and expose max_elevation_m",
            "wikidata_label_fetch": label_fetch_status,
            "review_flags": [
                "needs_coord_review",
                "needs_elevation_review",
                "is_shared_peak",
            ],
        },
        "items": items,
    }
    args.output.write_text(
        json.dumps(output, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(items)} mountain candidates to {args.output}")


if __name__ == "__main__":
    main()
