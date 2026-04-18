from __future__ import annotations

from dataclasses import dataclass

from .models import SearchJumpTarget


@dataclass(frozen=True)
class SearchQuerySpec:
    raw: str
    normalized: str
    selector: str | None = None
    value: str = ""


def parse_search_query(query: str) -> SearchQuerySpec:
    text = str(query or "").strip()
    if not text:
        return SearchQuerySpec(raw="", normalized="")
    if "=" in text:
        key, value = text.split("=", 1)
        selector = key.strip().casefold()
        value_text = value.strip()
        if selector in {"label", "id"} and value_text:
            return SearchQuerySpec(
                raw=text,
                normalized=value_text.casefold(),
                selector=selector,
                value=value_text,
            )
    return SearchQuerySpec(raw=text, normalized=text.casefold(), value=text)


def _target_search_haystack(target: SearchJumpTarget) -> str:
    parts = [
        target.label,
        target.subtitle,
        target.object_key,
        target.command,
        target.kind,
    ]
    return " ".join(part.strip() for part in parts if part).casefold()


def search_target_matches_query(
    target: SearchJumpTarget,
    spec: SearchQuerySpec,
) -> bool:
    if not spec.normalized:
        return False
    if spec.selector == "label":
        return spec.normalized in target.label.casefold()
    if spec.selector == "id":
        haystack = " ".join(
            part.strip()
            for part in (target.object_key, target.command, target.label)
            if part
        ).casefold()
        return spec.normalized in haystack
    return spec.normalized in _target_search_haystack(target)

