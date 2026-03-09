from __future__ import annotations

import json
import unicodedata
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Iterable


def normalize_viewpoint_name(text: str) -> str:
    folded = unicodedata.normalize("NFKC", text).casefold().strip()
    folded = folded.replace("’", "").replace("'", "")
    folded = folded.replace("–", "-").replace("—", "-")
    return " ".join(folded.split())


def looks_like_qid_placeholder(text: str) -> bool:
    candidate = (text or "").strip()
    return candidate.startswith("Q") and candidate[1:].isdigit()


def ascii_fallback_name(text: str) -> str | None:
    source = (text or "").strip()
    if not source or looks_like_qid_placeholder(source):
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
    if not any(ch.isalpha() for ch in ascii_text):
        return None
    if not ascii_text or ascii_text == source:
        return None
    return ascii_text


def prefixed_viewpoint_name(kind: str, name: str) -> str:
    if kind == "tower":
        return f"t/{name}"
    if kind == "mountain":
        return f"m/{name}"
    return name


def split_prefixed_viewpoint(text: str) -> tuple[str, str] | None:
    source = (text or "").strip()
    if len(source) >= 3 and source[1] == "/" and source[0] in {"t", "T", "m", "M"}:
        kind = "tower" if source[0] in {"t", "T"} else "mountain"
        return kind, source[2:].strip()
    return None


@dataclass(frozen=True)
class Viewpoint:
    id: str
    qid: str
    kind: str
    name: str
    labels: dict[str, str]
    names: tuple[str, ...]
    latitude_deg: float
    longitude_deg: float
    height_m: float
    meta: dict[str, Any]

    @property
    def persistent_key(self) -> str:
        return self.id

    @property
    def ascii_name(self) -> str | None:
        return ascii_fallback_name(self.name)


def build_viewpoint(
    item: dict[str, Any],
    *,
    kind: str,
    height_key: str,
    meta_keys: Iterable[str],
) -> Viewpoint:
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
    meta = {
        key: item[key]
        for key in meta_keys
        if key in item
    }
    return Viewpoint(
        id=str(item["id"]),
        qid=str(item["qid"]),
        kind=kind,
        name=str(item["name"]).strip(),
        labels=labels,
        names=names,
        latitude_deg=float(item["latitude_deg"]),
        longitude_deg=float(item["longitude_deg"]),
        height_m=float(item.get(height_key, 0.0)),
        meta=meta,
    )


def load_viewpoints(
    path: str | Path,
    *,
    builder: Callable[[dict[str, Any]], Viewpoint],
    error_label: str,
) -> tuple[Viewpoint, ...]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    items = payload.get("items", [])
    if not isinstance(items, list):
        raise ValueError(f"{error_label} file does not contain an items list")
    return tuple(builder(item) for item in items if isinstance(item, dict))


def list_viewpoint_primary_names(
    viewpoints: tuple[Viewpoint, ...],
) -> tuple[str, ...]:
    names = {
        (viewpoint.ascii_name or viewpoint.name).strip()
        for viewpoint in viewpoints
        if viewpoint.name.strip() and not looks_like_qid_placeholder(viewpoint.name)
    }
    return tuple(sorted(names, key=lambda name: (normalize_viewpoint_name(name), name)))


def list_viewpoint_all_names(
    viewpoints: tuple[Viewpoint, ...],
) -> tuple[str, ...]:
    names_by_key: dict[str, str] = {}
    for viewpoint in viewpoints:
        ascii_candidates = tuple(
            ascii_candidate
            for ascii_candidate in (
                viewpoint.ascii_name,
                *(ascii_fallback_name(candidate) for candidate in viewpoint.names),
                *(ascii_fallback_name(candidate) for candidate in viewpoint.labels.values()),
            )
            if ascii_candidate
        )
        for candidate in (
            viewpoint.name,
            *viewpoint.names,
            *viewpoint.labels.values(),
            *ascii_candidates,
        ):
            text = candidate.strip()
            if not text or looks_like_qid_placeholder(text):
                continue
            if " / " in text:
                continue
            normalized = normalize_viewpoint_name(text)
            existing = names_by_key.get(normalized)
            if existing is None or text < existing:
                names_by_key[normalized] = text
    return tuple(
        sorted(
            names_by_key.values(),
            key=lambda name: (normalize_viewpoint_name(name), name),
        )
    )


def viewpoint_to_dict(viewpoint: Viewpoint) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "id": viewpoint.id,
        "qid": viewpoint.qid,
        "kind": viewpoint.kind,
        "name": viewpoint.name,
        "labels": dict(viewpoint.labels),
        "names": list(viewpoint.names),
        "latitude_deg": viewpoint.latitude_deg,
        "longitude_deg": viewpoint.longitude_deg,
        "height_m": viewpoint.height_m,
        "meta": dict(viewpoint.meta),
    }
    if viewpoint.ascii_name is not None:
        payload["ascii_name"] = viewpoint.ascii_name
    return payload


def resolve_viewpoint(
    query: str,
    viewpoints: tuple[Viewpoint, ...],
    *,
    rank_key: Callable[[Viewpoint], float],
) -> Viewpoint | None:
    text = (query or "").strip()
    if not text:
        return None

    text_qid = text.removeprefix("wikidata:").strip()
    if text_qid.startswith("Q") and text_qid[1:].isdigit():
        for viewpoint in viewpoints:
            if viewpoint.qid == text_qid:
                return viewpoint
        return None

    normalized = normalize_viewpoint_name(text)
    allow_partial_match = not normalized.isdigit()
    exact_matches: list[Viewpoint] = []
    partial_matches: list[Viewpoint] = []
    for viewpoint in viewpoints:
        ascii_candidates = {
            ascii_candidate
            for ascii_candidate in (
                viewpoint.ascii_name,
                *(ascii_fallback_name(candidate) for candidate in viewpoint.names),
                *(ascii_fallback_name(candidate) for candidate in viewpoint.labels.values()),
            )
            if ascii_candidate
        }
        candidates = {
            viewpoint.name,
            *viewpoint.names,
            *viewpoint.labels.values(),
            *ascii_candidates,
        }
        normalized_candidates = {
            normalize_viewpoint_name(candidate)
            for candidate in candidates
            if candidate.strip()
        }
        if normalized in normalized_candidates:
            exact_matches.append(viewpoint)
            continue
        if allow_partial_match and any(normalized in candidate for candidate in normalized_candidates):
            partial_matches.append(viewpoint)

    if exact_matches:
        return max(exact_matches, key=rank_key)
    if len(partial_matches) == 1:
        return partial_matches[0]
    if partial_matches:
        return max(partial_matches, key=rank_key)
    return None


def find_exact_viewpoint_matches(
    query: str,
    viewpoints: tuple[Viewpoint, ...],
) -> tuple[Viewpoint, ...]:
    text = (query or "").strip()
    if not text:
        return ()
    normalized = normalize_viewpoint_name(text)
    matches: list[Viewpoint] = []
    for viewpoint in viewpoints:
        ascii_candidates = {
            ascii_candidate
            for ascii_candidate in (
                viewpoint.ascii_name,
                *(ascii_fallback_name(candidate) for candidate in viewpoint.names),
                *(ascii_fallback_name(candidate) for candidate in viewpoint.labels.values()),
            )
            if ascii_candidate
        }
        candidates = {
            viewpoint.name,
            *viewpoint.names,
            *viewpoint.labels.values(),
            *ascii_candidates,
        }
        normalized_candidates = {
            normalize_viewpoint_name(candidate)
            for candidate in candidates
            if candidate.strip()
        }
        if normalized in normalized_candidates:
            matches.append(viewpoint)
    return tuple(matches)
