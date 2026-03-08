from __future__ import annotations

import json
import unicodedata
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
        candidates = {tower.name, *tower.names, *tower.labels.values()}
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
