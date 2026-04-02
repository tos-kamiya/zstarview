#!/usr/bin/env python3
"""Build review-first tower viewpoint candidates from WDQS raw JSON files."""

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


POINT_RE = re.compile(r"^Point\(([-0-9.]+) ([-0-9.]+)\)$")
QID_RE = re.compile(r"/(Q\d+)$")
PREFERRED_LABEL_LANGS = ("en", "ja")
WIKIDATA_API_URL = "https://www.wikidata.org/w/api.php"
USER_AGENT = "zstarview-dev-samples/1.0 (tower review label fetch)"

DEFAULT_RULES_BY_BASENAME = {
    "query-1.json": (
        "instance_of_observation_tower",
        "instance_of_observation_deck",
    ),
    "query-2.json": ("has_use_observation_deck",),
    "query-3.json": (
        "instance_of_observation_tower",
        "instance_of_observation_deck",
    ),
}


@dataclass
class ReviewCandidate:
    qid: str
    name: str
    latitude_deg: float
    longitude_deg: float
    height_m: float
    viewpoint_height_m: float
    names: set[str] = field(default_factory=set)
    labels: dict[str, str] = field(default_factory=dict)
    matched_rules: set[str] = field(default_factory=set)
    source_files: set[str] = field(default_factory=set)

    def as_dict(self) -> dict[str, object]:
        slug = self.name.replace(" ", "_")
        return {
            "id": f"wikidata:{self.qid}",
            "name": self.name,
            "names": sorted(self.names or {self.name}),
            "labels": dict(sorted(self.labels.items())),
            "qid": self.qid,
            "latitude_deg": self.latitude_deg,
            "longitude_deg": self.longitude_deg,
            "height_m": self.height_m,
            "viewpoint_height_m": self.viewpoint_height_m,
            "matched_rules": sorted(self.matched_rules),
            "source_files": sorted(self.source_files),
            "viewpoint_types": ["tower_or_high_observation"],
            "wikidata_url": f"https://www.wikidata.org/wiki/{self.qid}",
            "location_arg": f"{self.latitude_deg:.6f};{self.longitude_deg:.6f}",
            "slug": slug,
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge multiple WDQS raw JSON files into a review-first tower "
            "candidate list with matched-rule provenance."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        action="append",
        default=[],
        help=(
            "Raw WDQS JSON array file. May be repeated. Known basenames such as "
            "query-1.json and query-2.json get default matched rules."
        ),
    )
    parser.add_argument(
        "--input-rule",
        action="append",
        default=[],
        metavar="PATH=RULE1,RULE2",
        help=(
            "Override matched rules for a specific input file. Example: "
            "raw-data/query-custom.json=has_part_observation_deck"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("dev-samples/wikidata_tower_review_candidates_2026-04-02.json"),
        help="Review-candidate JSON output path.",
    )
    parser.add_argument(
        "--min-height-m",
        type=float,
        default=100.0,
        help="Minimum height to keep after normalization.",
    )
    parser.add_argument(
        "--skip-wikidata-label-fetch",
        action="store_true",
        help="Do not fetch missing en/ja labels from Wikidata API.",
    )
    return parser.parse_args()


def parse_qid(entity_url: str) -> str:
    match = QID_RE.search(entity_url)
    if not match:
        raise ValueError(f"Could not parse QID from {entity_url!r}")
    return match.group(1)


def parse_point(point_text: str) -> tuple[float, float]:
    match = POINT_RE.match(point_text)
    if not match:
        raise ValueError(f"Could not parse coordinate point {point_text!r}")
    lon_text, lat_text = match.groups()
    return float(lat_text), float(lon_text)


def collect_labels(row: dict[str, Any]) -> tuple[dict[str, str], set[str]]:
    labels: dict[str, str] = {}
    names: set[str] = set()

    item_label = row.get("itemLabel")
    if isinstance(item_label, str) and item_label.strip():
        names.add(item_label.strip())

    explicit_item_labels = row.get("itemLabels")
    if isinstance(explicit_item_labels, dict):
        for lang, label in explicit_item_labels.items():
            if (
                isinstance(lang, str)
                and lang in PREFERRED_LABEL_LANGS
                and isinstance(label, str)
                and label.strip()
            ):
                labels[lang] = label.strip()
                names.add(label.strip())

    for key, value in row.items():
        if not isinstance(value, str) or not value.strip():
            continue
        if key.startswith("itemLabel_"):
            lang = key.removeprefix("itemLabel_")
            if lang in PREFERRED_LABEL_LANGS:
                labels[lang] = value.strip()
                names.add(value.strip())

    return labels, names


def choose_primary_name(labels: dict[str, str], names: set[str], fallback: str) -> str:
    for lang in PREFERRED_LABEL_LANGS:
        label = labels.get(lang)
        if label:
            return label
    if names:
        return sorted(names)[0]
    return fallback


def is_qid_like_name(name: str) -> bool:
    return bool(re.fullmatch(r"Q\d+", name))


def parse_input_rule_overrides(values: list[str]) -> dict[str, tuple[str, ...]]:
    overrides: dict[str, tuple[str, ...]] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Invalid --input-rule {value!r}; expected PATH=RULE1,RULE2")
        raw_path, raw_rules = value.split("=", 1)
        rules = tuple(rule.strip() for rule in raw_rules.split(",") if rule.strip())
        if not rules:
            raise ValueError(f"Invalid --input-rule {value!r}; no rules given")
        overrides[str(Path(raw_path))] = rules
    return overrides


def infer_rules_for_path(
    path: Path,
    overrides: dict[str, tuple[str, ...]],
) -> tuple[str, ...]:
    explicit = overrides.get(str(path))
    if explicit is not None:
        return explicit
    basename = path.name
    if basename in DEFAULT_RULES_BY_BASENAME:
        return DEFAULT_RULES_BY_BASENAME[basename]
    return ("candidate_from_query",)


def normalize_candidate_rows(
    rows_by_input: list[tuple[Path, list[dict[str, Any]]]],
    min_height_m: float,
    input_rule_overrides: dict[str, tuple[str, ...]] | None = None,
) -> list[dict[str, object]]:
    overrides = input_rule_overrides or {}
    records: dict[str, ReviewCandidate] = {}

    for input_path, rows in rows_by_input:
        matched_rules = infer_rules_for_path(input_path, overrides)
        for row in rows:
            item_url = row.get("item")
            coord = row.get("coord")
            if not isinstance(item_url, str) or not isinstance(coord, str):
                continue
            qid = parse_qid(item_url)
            latitude_deg, longitude_deg = parse_point(coord)

            height_text = row.get("height")
            try:
                height_m = float(height_text) if height_text not in (None, "") else 0.0
            except (TypeError, ValueError):
                height_m = 0.0
            if height_m < min_height_m:
                continue

            labels, names = collect_labels(row)
            name = choose_primary_name(labels, names, fallback=qid)

            record = records.get(qid)
            if record is None:
                record = ReviewCandidate(
                    qid=qid,
                    name=name,
                    latitude_deg=latitude_deg,
                    longitude_deg=longitude_deg,
                    height_m=height_m,
                    viewpoint_height_m=height_m,
                )
                records[qid] = record
            else:
                if height_m > record.height_m:
                    record.height_m = height_m
                    record.viewpoint_height_m = height_m
                if name and (
                    is_qid_like_name(record.name)
                    or (not is_qid_like_name(name) and len(name) > len(record.name))
                ):
                    record.name = name

            if names:
                record.names.update(names)
            elif not is_qid_like_name(name):
                record.names.add(name)
            for lang, label in labels.items():
                record.labels.setdefault(lang, label)
            record.matched_rules.update(matched_rules)
            record.source_files.add(str(input_path))

    return sorted(
        (record.as_dict() for record in records.values()),
        key=lambda item: (-float(item["height_m"]), str(item["name"])),
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
    labels_by_qid: dict[str, dict[str, str]],
) -> list[dict[str, object]]:
    enriched: list[dict[str, object]] = []

    for item in items:
        labels = {
            str(lang): str(label)
            for lang, label in dict(item.get("labels", {})).items()
            if isinstance(lang, str) and isinstance(label, str) and label.strip()
        }
        qid = str(item["qid"])
        fetched_labels = labels_by_qid.get(qid, {})
        for lang in PREFERRED_LABEL_LANGS:
            label = fetched_labels.get(lang)
            if label:
                labels[lang] = label

        names = {
            str(name).strip()
            for name in item.get("names", [])
            if isinstance(name, str) and str(name).strip()
        }
        names.update(labels.values())
        names.add(str(item["name"]).strip())

        updated = dict(item)
        updated["labels"] = dict(sorted(labels.items()))
        updated["names"] = sorted(names)
        updated["name"] = choose_primary_name(labels, names, fallback=str(item["name"]))
        updated["slug"] = str(updated["name"]).replace(" ", "_")
        enriched.append(updated)

    return enriched


def main() -> None:
    args = parse_args()
    if not args.input:
        raise SystemExit("At least one --input is required")

    input_rule_overrides = parse_input_rule_overrides(args.input_rule)
    rows_by_input = [
        (
            path,
            json.loads(path.read_text(encoding="utf-8")),
        )
        for path in args.input
    ]

    items = normalize_candidate_rows(
        rows_by_input,
        min_height_m=args.min_height_m,
        input_rule_overrides=input_rule_overrides,
    )
    label_fetch_status = "skipped"
    if not args.skip_wikidata_label_fetch and items:
        qids = [str(item["qid"]) for item in items]
        try:
            labels_by_qid = fetch_entity_labels(qids)
        except URLError as exc:
            print(f"Warning: failed to fetch Wikidata labels: {exc}")
            label_fetch_status = "failed"
        else:
            items = merge_entity_labels(items, labels_by_qid)
            label_fetch_status = "fetched"

    output = {
        "source": "wikidata",
        "source_query_results": [str(path) for path in args.input],
        "normalization": {
            "deduplicate_by": "qid",
            "height_rule": "max",
            "coord_rule": "first",
            "min_height_m": args.min_height_m,
            "matched_rules_rule": "distinct_sorted_union",
            "source_files_rule": "distinct_sorted_union",
            "preferred_label_languages": list(PREFERRED_LABEL_LANGS),
            "wikidata_label_fetch": label_fetch_status,
        },
        "items": items,
    }
    args.output.write_text(
        json.dumps(output, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(items)} tower review candidates to {args.output}")


if __name__ == "__main__":
    main()
