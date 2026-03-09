#!/usr/bin/env python3
"""Normalize curated mountain viewpoint candidates into zstarview JSON.

This script is the normalization step, not the discovery step.

Expected workflow:

1. Investigate candidate mountains in a `wikidata_mountain_query_investigation`
   memo and save any raw WDQS / Wikidata results under `dev-samples/`.
2. Hand-curate a seed JSON file from those results and supporting Wikipedia
   pages.
3. Run this script to normalize that seed into a `mountain_viewpoints` JSON.

The script enriches the curated seed with Wikidata labels, aliases, coordinate
location, and elevation when available.

Expected curated seed JSON shape:

{
  "candidates": [
    {
      "qid": "Q39261",
      "name": "Mount Fuji",
      "names": ["Fuji", "富士山"],
      "labels": {"ja": "富士山"},
      "latitude_deg": 35.360776,
      "longitude_deg": 138.727299,
      "elevation_m": 3776,
      "region_hint": "Japan",
      "country_codes": ["JP"],
      "wikipedia_titles": {"en": "Mount Fuji", "ja": "富士山"}
    }
  ]
}

Seed coordinates/elevation are treated as overrides. If omitted, the script will
fill them from Wikidata when available.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


PREFERRED_LABEL_LANGS = ("en", "ja")
PREFERRED_WIKIPEDIA_LANGS = ("en", "ja")
WIKIDATA_API_URL = "https://www.wikidata.org/w/api.php"
USER_AGENT = "zstarview-dev-samples/1.0 (mountain viewpoint builder)"
QID_RE = re.compile(r"^(?:wikidata:)?(Q\d+)$")
WIKIPEDIA_URL_TEMPLATE = "https://{lang}.wikipedia.org/wiki/{title}"
COORDINATE_PROPERTY = "P625"
ELEVATION_PROPERTY = "P2044"
LONG_SENTENCE_WORDS = 12
DROP_EXACT_NAMES = {
    "3003",
    "cold peak",
    "🗻",
    "オリンポス山 (ギリシャ)",
    'Desde la ruta 39, al km 47,5 pueden entrar en los senderos de "Las Cañas", y llegar hasta el Cerro Catedral. Sino entran en el km 60 y listo. Una maravilla!',
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a curated mountain-viewpoint JSON from a Wikipedia-seeded "
            "candidate list with optional Wikidata enrichment."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Curated candidate JSON path.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Normalized JSON output path.",
    )
    parser.add_argument(
        "--skip-wikidata-fetch",
        action="store_true",
        help="Do not fetch labels/aliases/coords/elevation from Wikidata API.",
    )
    parser.add_argument(
        "--allow-missing-elevation",
        action="store_true",
        help="Keep entries even if neither seed nor Wikidata provides elevation.",
    )
    return parser.parse_args()


def parse_qid(text: str) -> str:
    match = QID_RE.match((text or "").strip())
    if not match:
        raise ValueError(f"Could not parse QID from {text!r}")
    return match.group(1)


def _dedupe_sorted(values: set[str]) -> list[str]:
    return sorted(value for value in values if value.strip())


def _clean_alias_text(text: str) -> str:
    cleaned = (text or "").strip()
    cleaned = re.sub(r"^Mt\.?\s+", "Mount ", cleaned)
    cleaned = re.sub(r"^Mtn\.?\s+", "Mountain ", cleaned)
    cleaned = re.sub(r"\s+Mt\.?$", " Mountain", cleaned)
    cleaned = re.sub(r"\s+Mtn\.?$", " Mountain", cleaned)
    cleaned = re.sub(r"^(.*),\s*mount$", r"Mount \1", cleaned, flags=re.IGNORECASE)
    return " ".join(cleaned.split()).strip()


def _should_keep_alias(text: str) -> bool:
    cleaned = (text or "").strip()
    if not cleaned:
        return False
    if cleaned in DROP_EXACT_NAMES:
        return False
    if not any(ch.isalpha() for ch in cleaned):
        return False
    if " / " in cleaned:
        return False
    if len(cleaned.split()) >= LONG_SENTENCE_WORDS and re.search(r'[0-9,"!?.]', cleaned):
        return False
    return True


def _clean_aliases(names: set[str]) -> set[str]:
    cleaned_names: set[str] = set()
    for raw_name in names:
        cleaned = _clean_alias_text(raw_name)
        if not _should_keep_alias(cleaned):
            continue
        cleaned_names.add(cleaned)
    return cleaned_names


def _normalize_candidate_seed(candidate: dict[str, Any]) -> dict[str, Any]:
    qid = parse_qid(str(candidate["qid"]))
    names = {
        str(name).strip()
        for name in candidate.get("names", [])
        if isinstance(name, str) and str(name).strip()
    }
    labels_raw = candidate.get("labels", {})
    labels = {
        str(lang): str(label).strip()
        for lang, label in labels_raw.items()
        if isinstance(labels_raw, dict)
        and isinstance(lang, str)
        and isinstance(label, str)
        and str(label).strip()
    }
    names.update(labels.values())
    seed_name = str(candidate.get("name", "")).strip()
    if seed_name:
        names.add(seed_name)

    wikipedia_titles_raw = candidate.get("wikipedia_titles", {})
    wikipedia_titles = {
        str(lang): str(title).strip()
        for lang, title in wikipedia_titles_raw.items()
        if isinstance(wikipedia_titles_raw, dict)
        and isinstance(lang, str)
        and isinstance(title, str)
        and str(title).strip()
    }

    normalized: dict[str, Any] = {
        "qid": qid,
        "seed_name": seed_name,
        "seed_names": _dedupe_sorted(names),
        "seed_labels": dict(sorted(labels.items())),
        "wikipedia_titles": dict(sorted(wikipedia_titles.items())),
    }

    if candidate.get("latitude_deg") is not None:
        normalized["latitude_deg"] = float(candidate["latitude_deg"])
    if candidate.get("longitude_deg") is not None:
        normalized["longitude_deg"] = float(candidate["longitude_deg"])
    if candidate.get("elevation_m") is not None:
        normalized["elevation_m"] = float(candidate["elevation_m"])
    region_hint = str(candidate.get("region_hint", "")).strip()
    if region_hint:
        normalized["region_hint"] = region_hint
    country_codes = [
        str(code).strip().upper()
        for code in candidate.get("country_codes", [])
        if isinstance(code, str) and str(code).strip()
    ]
    if country_codes:
        normalized["country_codes"] = sorted(set(country_codes))
    return normalized


def normalize_seed_candidates(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, dict):
        raw_candidates = payload.get("candidates", [])
    else:
        raw_candidates = payload
    if not isinstance(raw_candidates, list):
        raise ValueError("Candidate payload must be a list or a dict containing 'candidates'.")

    deduped: dict[str, dict[str, Any]] = {}
    for raw_candidate in raw_candidates:
        if not isinstance(raw_candidate, dict):
            continue
        normalized = _normalize_candidate_seed(raw_candidate)
        existing = deduped.get(normalized["qid"])
        if existing is None:
            deduped[normalized["qid"]] = normalized
            continue
        existing["seed_names"] = _dedupe_sorted(set(existing["seed_names"]) | set(normalized["seed_names"]))
        merged_labels = dict(existing.get("seed_labels", {}))
        merged_labels.update(normalized.get("seed_labels", {}))
        existing["seed_labels"] = dict(sorted(merged_labels.items()))
        merged_titles = dict(existing.get("wikipedia_titles", {}))
        merged_titles.update(normalized.get("wikipedia_titles", {}))
        existing["wikipedia_titles"] = dict(sorted(merged_titles.items()))
        for key in ("seed_name", "region_hint"):
            if not existing.get(key) and normalized.get(key):
                existing[key] = normalized[key]
        if not existing.get("country_codes") and normalized.get("country_codes"):
            existing["country_codes"] = normalized["country_codes"]
        for key in ("latitude_deg", "longitude_deg", "elevation_m"):
            if existing.get(key) is None and normalized.get(key) is not None:
                existing[key] = normalized[key]

    return sorted(deduped.values(), key=lambda item: (str(item.get("seed_name", "")), str(item["qid"])))


def _extract_labels(entity: dict[str, Any]) -> dict[str, str]:
    raw_labels = entity.get("labels", {})
    if not isinstance(raw_labels, dict):
        return {}
    labels: dict[str, str] = {}
    for lang in PREFERRED_LABEL_LANGS:
        lang_value = raw_labels.get(lang)
        if isinstance(lang_value, dict):
            value = lang_value.get("value")
            if isinstance(value, str) and value.strip():
                labels[lang] = value.strip()
    return labels


def _extract_aliases(entity: dict[str, Any]) -> set[str]:
    raw_aliases = entity.get("aliases", {})
    if not isinstance(raw_aliases, dict):
        return set()
    names: set[str] = set()
    for lang_list in raw_aliases.values():
        if not isinstance(lang_list, list):
            continue
        for alias in lang_list:
            if not isinstance(alias, dict):
                continue
            value = alias.get("value")
            if isinstance(value, str) and value.strip():
                names.add(value.strip())
    return names


def _extract_coordinate(entity: dict[str, Any]) -> tuple[float, float] | None:
    claims = entity.get("claims", {})
    if not isinstance(claims, dict):
        return None
    claim_list = claims.get(COORDINATE_PROPERTY, [])
    if not isinstance(claim_list, list):
        return None
    for claim in claim_list:
        mainsnak = claim.get("mainsnak")
        if not isinstance(mainsnak, dict):
            continue
        datavalue = mainsnak.get("datavalue")
        if not isinstance(datavalue, dict):
            continue
        value = datavalue.get("value")
        if not isinstance(value, dict):
            continue
        latitude = value.get("latitude")
        longitude = value.get("longitude")
        if isinstance(latitude, (int, float)) and isinstance(longitude, (int, float)):
            return float(latitude), float(longitude)
    return None


def _extract_elevation(entity: dict[str, Any]) -> float | None:
    claims = entity.get("claims", {})
    if not isinstance(claims, dict):
        return None
    claim_list = claims.get(ELEVATION_PROPERTY, [])
    if not isinstance(claim_list, list):
        return None
    for claim in claim_list:
        mainsnak = claim.get("mainsnak")
        if not isinstance(mainsnak, dict):
            continue
        datavalue = mainsnak.get("datavalue")
        if not isinstance(datavalue, dict):
            continue
        value = datavalue.get("value")
        if not isinstance(value, dict):
            continue
        amount = value.get("amount")
        if isinstance(amount, str):
            try:
                return float(amount)
            except ValueError:
                continue
        if isinstance(amount, (int, float)):
            return float(amount)
    return None


def _extract_wikipedia_urls(entity: dict[str, Any]) -> dict[str, str]:
    sitelinks = entity.get("sitelinks", {})
    if not isinstance(sitelinks, dict):
        return {}
    urls: dict[str, str] = {}
    for lang in PREFERRED_WIKIPEDIA_LANGS:
        site = sitelinks.get(f"{lang}wiki")
        if not isinstance(site, dict):
            continue
        title = site.get("title")
        if isinstance(title, str) and title.strip():
            urls[lang] = WIKIPEDIA_URL_TEMPLATE.format(lang=lang, title=title.replace(" ", "_"))
    return urls


def fetch_entity_data(
    qids: list[str],
    batch_size: int = 50,
) -> dict[str, dict[str, Any]]:
    entities_by_qid: dict[str, dict[str, Any]] = {}

    for start_index in range(0, len(qids), batch_size):
        batch = qids[start_index : start_index + batch_size]
        params = urlencode(
            {
                "action": "wbgetentities",
                "ids": "|".join(batch),
                "props": "labels|aliases|sitelinks|claims",
                "languages": "|".join(PREFERRED_LABEL_LANGS),
                "sitefilter": "|".join(f"{lang}wiki" for lang in PREFERRED_WIKIPEDIA_LANGS),
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
            labels = _extract_labels(entity)
            aliases = _extract_aliases(entity)
            coordinate = _extract_coordinate(entity)
            elevation_m = _extract_elevation(entity)
            wikipedia_urls = _extract_wikipedia_urls(entity)
            entities_by_qid[qid] = {
                "labels": labels,
                "aliases": aliases,
                "coordinate": coordinate,
                "elevation_m": elevation_m,
                "wikipedia_urls": wikipedia_urls,
            }

    return entities_by_qid


def choose_primary_name(
    *,
    seed_name: str,
    labels: dict[str, str],
    names: set[str],
) -> str:
    seed_name = _clean_alias_text(seed_name)
    for lang in PREFERRED_LABEL_LANGS:
        label = labels.get(lang)
        if label:
            return _clean_alias_text(label)
    if seed_name:
        return seed_name
    if names:
        return sorted(names)[0]
    raise ValueError("Could not determine a primary name for mountain candidate")


def merge_entity_data(
    candidates: list[dict[str, Any]],
    entities_by_qid: dict[str, dict[str, Any]],
    *,
    allow_missing_elevation: bool,
) -> list[dict[str, Any]]:
    mountains: list[dict[str, Any]] = []

    for candidate in candidates:
        qid = str(candidate["qid"])
        fetched = entities_by_qid.get(qid, {})

        labels = dict(candidate.get("seed_labels", {}))
        labels.update(fetched.get("labels", {}))
        labels = {
            lang: cleaned
            for lang, label in labels.items()
            if (cleaned := _clean_alias_text(str(label))) and _should_keep_alias(cleaned)
        }

        names = set(candidate.get("seed_names", []))
        names.update(labels.values())
        names.update(fetched.get("aliases", set()))
        names = _clean_aliases(names)

        primary_name = choose_primary_name(
            seed_name=str(candidate.get("seed_name", "")),
            labels=labels,
            names=names,
        )

        latitude_deg = candidate.get("latitude_deg")
        longitude_deg = candidate.get("longitude_deg")
        fetched_coordinate = fetched.get("coordinate")
        if latitude_deg is None or longitude_deg is None:
            if isinstance(fetched_coordinate, tuple):
                latitude_deg, longitude_deg = fetched_coordinate
        if latitude_deg is None or longitude_deg is None:
            raise ValueError(f"Missing coordinate for {qid}")

        elevation_m = candidate.get("elevation_m")
        if elevation_m is None:
            elevation_m = fetched.get("elevation_m")
        if elevation_m is None and not allow_missing_elevation:
            raise ValueError(f"Missing elevation for {qid}")

        wikipedia_urls: dict[str, str] = {}
        for lang, title in candidate.get("wikipedia_titles", {}).items():
            wikipedia_urls[lang] = WIKIPEDIA_URL_TEMPLATE.format(lang=lang, title=str(title).replace(" ", "_"))
        wikipedia_urls.update(fetched.get("wikipedia_urls", {}))

        mountain = {
            "id": f"wikidata:{qid}",
            "qid": qid,
            "name": primary_name,
            "names": _dedupe_sorted(names),
            "labels": dict(sorted(labels.items())),
            "latitude_deg": float(latitude_deg),
            "longitude_deg": float(longitude_deg),
            "wikidata_url": f"https://www.wikidata.org/wiki/{qid}",
            "wikipedia_urls": dict(sorted(wikipedia_urls.items())),
            "location_arg": f"{float(latitude_deg):.6f};{float(longitude_deg):.6f}",
            "slug": primary_name.replace(" ", "_"),
        }
        if elevation_m is not None:
            mountain["elevation_m"] = float(elevation_m)
        if candidate.get("region_hint"):
            mountain["region_hint"] = str(candidate["region_hint"])
        if candidate.get("country_codes"):
            mountain["country_codes"] = list(candidate["country_codes"])
        mountains.append(mountain)

    return sorted(mountains, key=lambda item: (str(item["name"]), str(item["qid"])))


def main() -> None:
    args = parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    candidates = normalize_seed_candidates(payload)
    fetch_status = "skipped"
    entities_by_qid: dict[str, dict[str, Any]] = {}

    if not args.skip_wikidata_fetch and candidates:
        try:
            entities_by_qid = fetch_entity_data([str(candidate["qid"]) for candidate in candidates])
        except URLError as exc:
            print(f"Warning: failed to fetch Wikidata entity data: {exc}")
            fetch_status = "failed"
        else:
            fetch_status = "fetched"

    mountains = merge_entity_data(
        candidates,
        entities_by_qid,
        allow_missing_elevation=args.allow_missing_elevation,
    )

    output = {
        "source": "wikipedia+wikidata-curated",
        "source_seed": str(args.input),
        "normalization": {
            "deduplicate_by": "qid",
            "preferred_label_languages": list(PREFERRED_LABEL_LANGS),
            "preferred_wikipedia_languages": list(PREFERRED_WIKIPEDIA_LANGS),
            "coordinate_rule": "seed override else Wikidata P625",
            "elevation_rule": "seed override else Wikidata P2044",
            "wikidata_fetch": fetch_status,
        },
        "items": mountains,
    }
    args.output.write_text(
        json.dumps(output, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(mountains)} mountain viewpoints to {args.output}")


if __name__ == "__main__":
    main()
