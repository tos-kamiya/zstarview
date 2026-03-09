from __future__ import annotations

import json
import unicodedata
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from .paths import MOUNTAIN_VIEWPOINTS_FILE


def _normalize_name(text: str) -> str:
    folded = unicodedata.normalize("NFKC", text).casefold().strip()
    folded = folded.replace("’", "").replace("'", "")
    folded = folded.replace("–", "-").replace("—", "-")
    return " ".join(folded.split())


def _looks_like_qid_placeholder(text: str) -> bool:
    candidate = (text or "").strip()
    return candidate.startswith("Q") and candidate[1:].isdigit()


def _ascii_fallback_name(text: str) -> str | None:
    source = (text or "").strip()
    if not source or _looks_like_qid_placeholder(source):
        return None
    replaced = (
        source.replace("’", "'")
        .replace("‘", "'")
        .replace("“", '"')
        .replace("”", '"')
        .replace("–", "-")
        .replace("—", "-")
    )
    decomposed = unicodedata.normalize("NFKD", replaced)
    ascii_text = "".join(
        ch for ch in decomposed
        if not unicodedata.combining(ch) and ord(ch) < 128
    )
    ascii_text = " ".join(ascii_text.split()).strip()
    if not ascii_text or ascii_text == source:
        return None
    return ascii_text


@dataclass(frozen=True)
class MountainViewpoint:
    id: str
    qid: str
    name: str
    labels: dict[str, str]
    names: tuple[str, ...]
    latitude_deg: float
    longitude_deg: float
    elevation_m: float
    wikidata_url: str
    wikipedia_urls: dict[str, str]
    location_arg: str
    slug: str
    region_hint: str

    @property
    def persistent_key(self) -> str:
        return self.id

    @property
    def ascii_name(self) -> str | None:
        return _ascii_fallback_name(self.name)


def _as_mountain(item: dict[str, Any]) -> MountainViewpoint:
    labels_raw = item.get("labels", {})
    labels = {
        str(lang): str(label)
        for lang, label in labels_raw.items()
        if isinstance(lang, str) and isinstance(label, str) and label.strip()
    } if isinstance(labels_raw, dict) else {}
    names_raw = item.get("names", [])
    names = tuple(
        str(name).strip()
        for name in names_raw
        if isinstance(name, str) and str(name).strip()
    )
    wikipedia_urls_raw = item.get("wikipedia_urls", {})
    wikipedia_urls = {
        str(lang): str(url)
        for lang, url in wikipedia_urls_raw.items()
        if isinstance(lang, str) and isinstance(url, str) and url.strip()
    } if isinstance(wikipedia_urls_raw, dict) else {}
    return MountainViewpoint(
        id=str(item["id"]),
        qid=str(item["qid"]),
        name=str(item["name"]).strip(),
        labels=labels,
        names=names,
        latitude_deg=float(item["latitude_deg"]),
        longitude_deg=float(item["longitude_deg"]),
        elevation_m=float(item["elevation_m"]),
        wikidata_url=str(item.get("wikidata_url", "")).strip(),
        wikipedia_urls=wikipedia_urls,
        location_arg=str(item.get("location_arg", "")).strip(),
        slug=str(item.get("slug", "")).strip(),
        region_hint=str(item.get("region_hint", "")).strip(),
    )


@lru_cache(maxsize=1)
def load_mountain_viewpoints(
    path: str | Path = MOUNTAIN_VIEWPOINTS_FILE,
) -> tuple[MountainViewpoint, ...]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    items = payload.get("items", [])
    if not isinstance(items, list):
        raise ValueError("mountain viewpoints file does not contain an items list")
    mountains = tuple(_as_mountain(item) for item in items if isinstance(item, dict))
    return mountains


def list_mountain_primary_names(
    mountains: tuple[MountainViewpoint, ...] | None = None,
) -> tuple[str, ...]:
    mountains = load_mountain_viewpoints() if mountains is None else mountains
    names = {
        (mountain.ascii_name or mountain.name).strip()
        for mountain in mountains
        if mountain.name.strip() and not _looks_like_qid_placeholder(mountain.name)
    }
    return tuple(sorted(names, key=lambda name: (_normalize_name(name), name)))


def list_mountain_all_names(
    mountains: tuple[MountainViewpoint, ...] | None = None,
) -> tuple[str, ...]:
    mountains = load_mountain_viewpoints() if mountains is None else mountains
    names_by_key: dict[str, str] = {}
    for mountain in mountains:
        ascii_candidates = tuple(
            ascii_candidate
            for ascii_candidate in (
                mountain.ascii_name,
                *(_ascii_fallback_name(candidate) for candidate in mountain.names),
                *(_ascii_fallback_name(candidate) for candidate in mountain.labels.values()),
            )
            if ascii_candidate
        )
        for candidate in (
            mountain.name,
            *mountain.names,
            *mountain.labels.values(),
            *ascii_candidates,
        ):
            text = candidate.strip()
            if not text:
                continue
            if _looks_like_qid_placeholder(text):
                continue
            normalized = _normalize_name(text)
            existing = names_by_key.get(normalized)
            if existing is None or text < existing:
                names_by_key[normalized] = text
    return tuple(sorted(names_by_key.values(), key=lambda name: (_normalize_name(name), name)))


def mountain_viewpoint_to_dict(mountain: MountainViewpoint) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "id": mountain.id,
        "qid": mountain.qid,
        "name": mountain.name,
        "labels": dict(mountain.labels),
        "names": list(mountain.names),
        "latitude_deg": mountain.latitude_deg,
        "longitude_deg": mountain.longitude_deg,
        "elevation_m": mountain.elevation_m,
        "wikidata_url": mountain.wikidata_url,
        "wikipedia_urls": dict(mountain.wikipedia_urls),
        "location_arg": mountain.location_arg,
        "slug": mountain.slug,
        "region_hint": mountain.region_hint,
    }
    if mountain.ascii_name is not None:
        payload["ascii_name"] = mountain.ascii_name
    return payload


def resolve_mountain_viewpoint(
    query: str,
    mountains: tuple[MountainViewpoint, ...] | None = None,
) -> MountainViewpoint | None:
    mountains = load_mountain_viewpoints() if mountains is None else mountains
    text = (query or "").strip()
    if not text:
        return None

    text_qid = text.removeprefix("wikidata:").strip()
    if text_qid.startswith("Q") and text_qid[1:].isdigit():
        for mountain in mountains:
            if mountain.qid == text_qid:
                return mountain
        return None

    normalized = _normalize_name(text)
    exact_matches: list[MountainViewpoint] = []
    partial_matches: list[MountainViewpoint] = []
    for mountain in mountains:
        ascii_candidates = {
            ascii_candidate
            for ascii_candidate in (
                mountain.ascii_name,
                *(_ascii_fallback_name(candidate) for candidate in mountain.names),
                *(_ascii_fallback_name(candidate) for candidate in mountain.labels.values()),
            )
            if ascii_candidate
        }
        candidates = {mountain.name, *mountain.names, *mountain.labels.values(), *ascii_candidates}
        normalized_candidates = {_normalize_name(candidate) for candidate in candidates if candidate.strip()}
        if normalized in normalized_candidates:
            exact_matches.append(mountain)
            continue
        if any(normalized in candidate for candidate in normalized_candidates):
            partial_matches.append(mountain)

    if exact_matches:
        return max(exact_matches, key=lambda mountain: mountain.elevation_m)
    if len(partial_matches) == 1:
        return partial_matches[0]
    if partial_matches:
        return max(partial_matches, key=lambda mountain: mountain.elevation_m)
    return None
