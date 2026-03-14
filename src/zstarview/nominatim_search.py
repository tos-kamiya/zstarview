from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request


_VERSION = "0.2.0"
_BASE_URL = "https://nominatim.openstreetmap.org/search"
_DEFAULT_USER_AGENT = f"zstarview/{_VERSION}"


def _build_url(query: str, limit: int, countrycode: str | None) -> str:
    params = {
        "q": query,
        "format": "jsonv2",
        "limit": str(limit),
        "addressdetails": "0",
        "dedupe": "1",
    }
    if countrycode:
        params["countrycodes"] = countrycode
    return _BASE_URL + "?" + urllib.parse.urlencode(params)


def _fetch(url: str, language: str, user_agent: str = _DEFAULT_USER_AGENT) -> list[dict]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": user_agent,
            "Accept-Language": language,
        },
    )
    with urllib.request.urlopen(req, timeout=15) as resp:
        charset = resp.headers.get_content_charset() or "utf-8"
        body = resp.read().decode(charset)

    data = json.loads(body)
    if not isinstance(data, list):
        raise ValueError("unexpected response format")
    return data


def _normalize(items: list[dict]) -> list[dict]:
    results = []
    for item in items:
        try:
            lat = float(item["lat"])
            lon = float(item["lon"])
        except (KeyError, TypeError, ValueError):
            continue

        importance = item.get("importance") or 0.0
        results.append(
            {
                "name": item.get("display_name"),
                "lat": lat,
                "lon": lon,
                "category": item.get("category") or item.get("class"),
                "type": item.get("type"),
                "importance": float(importance),
            }
        )

    results.sort(key=lambda result: result["importance"], reverse=True)
    return results


def search(
    query: str,
    *,
    limit: int = 5,
    countrycode: str | None = None,
    language: str = "en",
    user_agent: str = _DEFAULT_USER_AGENT,
) -> list[dict]:
    url = _build_url(query, limit=limit, countrycode=countrycode)
    return _normalize(_fetch(url, language=language, user_agent=user_agent))


__all__ = [
    "search",
]
