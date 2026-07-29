#!/usr/bin/env python3
"""Merge existing bundled viewpoints with extra WDQS high-observation candidates."""

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


@dataclass
class CandidateRecord:
    qid: str
    name: str
    latitude_deg: float
    longitude_deg: float
    height_m: float
    names: set[str] = field(default_factory=set)

    def as_item(self) -> dict[str, object]:
        names = sorted(self.names or {self.name})
        return {
            "id": f"wikidata:{self.qid}",
            "name": self.name,
            "names": names,
            "labels": {},
            "qid": self.qid,
            "latitude_deg": self.latitude_deg,
            "longitude_deg": self.longitude_deg,
            "height_m": self.height_m,
            "viewpoint_height_m": self.height_m,
            "viewpoint_types": ["tower_or_high_observation"],
            "wikidata_url": f"https://www.wikidata.org/wiki/{self.qid}",
            "sources": [
                {
                    "field": "name, names, labels, latitude_deg, longitude_deg, height_m",
                    "name": "Wikidata",
                    "url": f"https://www.wikidata.org/wiki/{self.qid}",
                }
            ],
            "location_arg": f"{self.latitude_deg:.6f};{self.longitude_deg:.6f}",
            "slug": self.name.replace(" ", "_"),
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge existing bundled tower viewpoints with extra Wikidata "
            "high-observation candidates."
        )
    )
    parser.add_argument(
        "--base",
        type=Path,
        default=Path("src/zstarview/data/viewpoints/tower_viewpoints.json"),
        help="Existing bundled viewpoints JSON to use as the base dataset.",
    )
    parser.add_argument(
        "--extra-query",
        action="append",
        default=[],
        help="Extra WDQS JSON array file such as query-2.json. May be repeated.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("src/zstarview/data/viewpoints/tower_viewpoints.json"),
        help="Where to write the merged bundled viewpoints JSON.",
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


def normalize_extra_rows(rows: list[dict[str, Any]]) -> list[dict[str, object]]:
    records: dict[str, CandidateRecord] = {}
    for row in rows:
        item_url = row.get("item")
        coord = row.get("coord")
        if not isinstance(item_url, str) or not isinstance(coord, str):
            continue
        qid = parse_qid(item_url)
        latitude_deg, longitude_deg = parse_point(coord)
        name = str(row.get("itemLabel") or qid).strip()
        height_text = row.get("height")
        try:
            height_m = float(height_text) if height_text not in (None, "") else 0.0
        except (TypeError, ValueError):
            height_m = 0.0
        if height_m <= 0.0:
            continue

        record = records.get(qid)
        if record is None:
            record = CandidateRecord(
                qid=qid,
                name=name,
                latitude_deg=latitude_deg,
                longitude_deg=longitude_deg,
                height_m=height_m,
            )
            records[qid] = record
        else:
            record.height_m = max(record.height_m, height_m)
            if name and len(name) > len(record.name):
                record.name = name
        if name:
            record.names.add(name)

    return sorted(
        (record.as_item() for record in records.values()),
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


def choose_primary_name(labels: dict[str, str], names: set[str]) -> str:
    for lang in PREFERRED_LABEL_LANGS:
        label = labels.get(lang)
        if label:
            return label
    if names:
        return sorted(names)[0]
    raise ValueError("Could not determine a primary label for item")


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
        if item.get("name"):
            names.add(str(item["name"]).strip())

        updated = dict(item)
        updated["labels"] = dict(sorted(labels.items()))
        updated["names"] = sorted(names)
        updated["name"] = choose_primary_name(labels, names)
        updated["slug"] = str(updated["name"]).replace(" ", "_")
        enriched.append(updated)

    return enriched


def merge_items(
    base_items: list[dict[str, object]],
    extra_items: list[dict[str, object]],
) -> list[dict[str, object]]:
    merged: dict[str, dict[str, object]] = {
        str(item["qid"]): dict(item)
        for item in base_items
        if isinstance(item, dict) and "qid" in item
    }

    for item in merged.values():
        if "viewpoint_height_m" not in item:
            if "observer_height_m" in item:
                item["viewpoint_height_m"] = float(item.get("observer_height_m", 0.0))
            else:
                item["viewpoint_height_m"] = float(item.get("height_m", 0.0))
        item.pop("observer_height_m", None)

    for extra in extra_items:
        qid = str(extra["qid"])
        current = merged.get(qid)
        if current is None:
            merged[qid] = dict(extra)
            continue

        current_names = {
            str(name).strip()
            for name in current.get("names", [])
            if isinstance(name, str) and str(name).strip()
        }
        extra_names = {
            str(name).strip()
            for name in extra.get("names", [])
            if isinstance(name, str) and str(name).strip()
        }
        current_names.update(extra_names)
        if current.get("name"):
            current_names.add(str(current["name"]).strip())
        if extra.get("name"):
            current_names.add(str(extra["name"]).strip())
        current["names"] = sorted(name for name in current_names if name)

        current_height = float(current.get("height_m", 0.0))
        extra_height = float(extra.get("height_m", 0.0))
        if extra_height > current_height:
            current["height_m"] = extra_height
        if "viewpoint_height_m" not in current:
            if "observer_height_m" in current:
                current["viewpoint_height_m"] = current.get("observer_height_m", 0.0)
            else:
                current["viewpoint_height_m"] = current.get("height_m", 0.0)
        current.pop("observer_height_m", None)

        current_types = {
            str(value)
            for value in current.get("viewpoint_types", [])
            if isinstance(value, str) and value.strip()
        }
        extra_types = {
            str(value)
            for value in extra.get("viewpoint_types", [])
            if isinstance(value, str) and value.strip()
        }
        if extra_types:
            current["viewpoint_types"] = sorted(current_types | extra_types)

    return sorted(
        merged.values(),
        key=lambda item: (-float(item.get("viewpoint_height_m", 0.0)), str(item.get("name", ""))),
    )


def load_rows(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise ValueError(f"{path} must contain a top-level JSON array")
    return [row for row in payload if isinstance(row, dict)]


def main() -> int:
    args = parse_args()
    base_payload = json.loads(args.base.read_text(encoding="utf-8"))
    base_items = base_payload.get("items", [])
    if not isinstance(base_items, list):
        raise ValueError("Base viewpoints file does not contain an items list.")

    extra_items: list[dict[str, object]] = []
    for extra_path_text in args.extra_query:
        extra_path = Path(extra_path_text)
        extra_items.extend(normalize_extra_rows(load_rows(extra_path)))

    merged_items = merge_items(
        [item for item in base_items if isinstance(item, dict)],
        extra_items,
    )
    if not args.skip_wikidata_label_fetch and merged_items:
        qids = [str(item["qid"]) for item in merged_items]
        try:
            labels_by_qid = fetch_entity_labels(qids)
        except URLError as exc:
            print(f"warning: failed to fetch Wikidata labels: {exc}")
        else:
            merged_items = merge_entity_labels(merged_items, labels_by_qid)

    output_payload = dict(base_payload)
    output_payload["items"] = merged_items
    output_payload["source_query_results"] = [
        str(Path(path_text))
        for path_text in args.extra_query
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(output_payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"[ok] wrote {args.output} items={len(merged_items)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
