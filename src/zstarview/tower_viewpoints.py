from __future__ import annotations

import json
import unicodedata
from dataclasses import asdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from .paths import TOWER_VIEWPOINTS_FILE


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
class TowerViewpoint:
    qid: str
    name: str
    labels: dict[str, str]
    names: tuple[str, ...]
    latitude_deg: float
    longitude_deg: float
    height_m: float
    classes: tuple[str, ...]

    @property
    def persistent_key(self) -> str:
        return f"wikidata:{self.qid}"

    @property
    def ascii_name(self) -> str | None:
        return _ascii_fallback_name(self.name)


def _as_tower(item: dict[str, Any]) -> TowerViewpoint:
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
    classes_raw = item.get("classes", [])
    classes = tuple(
        str(name).strip()
        for name in classes_raw
        if isinstance(name, str) and str(name).strip()
    )
    return TowerViewpoint(
        qid=str(item["qid"]),
        name=str(item["name"]).strip(),
        labels=labels,
        names=names,
        latitude_deg=float(item["latitude_deg"]),
        longitude_deg=float(item["longitude_deg"]),
        height_m=float(item["height_m"]),
        classes=classes,
    )


@lru_cache(maxsize=1)
def load_tower_viewpoints(path: str | Path = TOWER_VIEWPOINTS_FILE) -> tuple[TowerViewpoint, ...]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    items = payload.get("items", [])
    if not isinstance(items, list):
        raise ValueError("tower viewpoints file does not contain an items list")
    towers = tuple(_as_tower(item) for item in items if isinstance(item, dict))
    return towers


def list_tower_primary_names(towers: tuple[TowerViewpoint, ...] | None = None) -> tuple[str, ...]:
    towers = load_tower_viewpoints() if towers is None else towers
    names = {
        (tower.ascii_name or tower.name).strip()
        for tower in towers
        if tower.name.strip() and not _looks_like_qid_placeholder(tower.name)
    }
    return tuple(sorted(names, key=lambda name: (_normalize_name(name), name)))


def list_tower_all_names(towers: tuple[TowerViewpoint, ...] | None = None) -> tuple[str, ...]:
    towers = load_tower_viewpoints() if towers is None else towers
    names_by_key: dict[str, str] = {}
    for tower in towers:
        ascii_candidates = tuple(
            ascii_candidate
            for ascii_candidate in (
                tower.ascii_name,
                *(_ascii_fallback_name(candidate) for candidate in tower.names),
                *(_ascii_fallback_name(candidate) for candidate in tower.labels.values()),
            )
            if ascii_candidate
        )
        for candidate in (tower.name, *tower.names, *tower.labels.values(), *ascii_candidates):
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


def tower_viewpoint_to_dict(tower: TowerViewpoint) -> dict[str, Any]:
    payload = asdict(tower)
    payload["names"] = list(tower.names)
    payload["classes"] = list(tower.classes)
    if tower.ascii_name is not None:
        payload["ascii_name"] = tower.ascii_name
    return payload


def resolve_tower_viewpoint(query: str, towers: tuple[TowerViewpoint, ...] | None = None) -> TowerViewpoint | None:
    towers = load_tower_viewpoints() if towers is None else towers
    text = (query or "").strip()
    if not text:
        return None

    text_qid = text.removeprefix("wikidata:").strip()
    if text_qid.startswith("Q") and text_qid[1:].isdigit():
        for tower in towers:
            if tower.qid == text_qid:
                return tower
        return None

    normalized = _normalize_name(text)
    exact_matches: list[TowerViewpoint] = []
    partial_matches: list[TowerViewpoint] = []
    for tower in towers:
        ascii_candidates = {
            ascii_candidate
            for ascii_candidate in (
                tower.ascii_name,
                *(_ascii_fallback_name(candidate) for candidate in tower.names),
                *(_ascii_fallback_name(candidate) for candidate in tower.labels.values()),
            )
            if ascii_candidate
        }
        candidates = {tower.name, *tower.names, *tower.labels.values(), *ascii_candidates}
        normalized_candidates = {_normalize_name(candidate) for candidate in candidates if candidate.strip()}
        if normalized in normalized_candidates:
            exact_matches.append(tower)
            continue
        if any(normalized in candidate for candidate in normalized_candidates):
            partial_matches.append(tower)

    if exact_matches:
        return max(exact_matches, key=lambda tower: tower.height_m)
    if len(partial_matches) == 1:
        return partial_matches[0]
    if partial_matches:
        return max(partial_matches, key=lambda tower: tower.height_m)
    return None
