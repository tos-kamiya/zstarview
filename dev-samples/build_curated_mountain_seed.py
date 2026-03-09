#!/usr/bin/env python3
"""Create a curated mountain seed JSON from a review candidate JSON."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


EXCLUDED_NAME_TOKENS = ("hill", "ridge", "chemin")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a curated mountain candidate seed from a review JSON using "
            "the currently approved filtering rules."
        )
    )
    parser.add_argument("--input", type=Path, required=True, help="Review JSON input path.")
    parser.add_argument("--output", type=Path, required=True, help="Curated seed JSON output path.")
    return parser.parse_args()


def should_exclude(item: dict[str, Any]) -> bool:
    if bool(item.get("needs_coord_review")):
        return True
    if bool(item.get("needs_elevation_review")):
        return True
    if bool(item.get("is_shared_peak")):
        return True

    names_to_check: list[str] = [str(item.get("item_label", ""))]
    item_labels = item.get("item_labels", {})
    if isinstance(item_labels, dict):
        names_to_check.extend(str(value) for value in item_labels.values())
    lowered = " ".join(name.casefold() for name in names_to_check if name.strip())
    return any(token in lowered for token in EXCLUDED_NAME_TOKENS)


def to_seed_candidate(item: dict[str, Any]) -> dict[str, Any]:
    item_labels = item.get("item_labels", {})
    labels = {
        str(lang): str(value).strip()
        for lang, value in item_labels.items()
        if isinstance(item_labels, dict) and isinstance(lang, str) and isinstance(value, str) and value.strip()
    }
    primary_name = labels.get("en") or str(item.get("item_label", "")).strip()
    names = {primary_name, str(item.get("item_label", "")).strip(), *labels.values()}
    region_hint = ""
    country_labels = item.get("country_labels", {})
    if isinstance(country_labels, dict):
        region_hint = str(country_labels.get("en") or country_labels.get("ja") or "").strip()
    if not region_hint:
        region_hint = str(item.get("country_label", "")).strip()

    candidate = {
        "qid": str(item["item_qid"]),
        "name": primary_name,
        "names": sorted(name for name in names if name),
        "labels": dict(sorted(labels.items())),
        "latitude_deg": float(item["latitude_deg"]),
        "longitude_deg": float(item["longitude_deg"]),
        "elevation_m": float(item["max_elevation_m"]),
        "region_hint": region_hint,
    }
    return candidate


def build_seed(payload: Any) -> dict[str, object]:
    if not isinstance(payload, dict):
        raise ValueError("Review payload must be a dict.")
    raw_items = payload.get("items", [])
    if not isinstance(raw_items, list):
        raise ValueError("Review payload must contain an 'items' list.")

    candidates = [
        to_seed_candidate(item)
        for item in raw_items
        if isinstance(item, dict) and not should_exclude(item)
    ]
    candidates.sort(key=lambda item: (str(item["region_hint"]), str(item["name"])))
    return {
        "source": "wikidata-mountain-review",
        "source_review_result": str(payload.get("source_query_result", "")),
        "curation_rules": {
            "exclude_needs_coord_review": True,
            "exclude_needs_elevation_review": True,
            "exclude_shared_peaks": True,
            "exclude_name_tokens": list(EXCLUDED_NAME_TOKENS),
            "name_rule": "prefer labels.en, else item_label",
            "region_hint_rule": "prefer country_labels.en, else country_label",
        },
        "candidates": candidates,
    }


def main() -> None:
    args = parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    output = build_seed(payload)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(output['candidates'])} curated mountain candidates to {args.output}")


if __name__ == "__main__":
    main()
