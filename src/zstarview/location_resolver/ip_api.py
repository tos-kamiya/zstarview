from __future__ import annotations

import json
import logging
import urllib.request
from typing import Any

logger = logging.getLogger(__name__)

_API_URL = "http://ip-api.com/json"


def fetch_location_by_ip() -> dict[str, Any]:
    """Fetch current location data from ip-api.com."""
    logger.info("Fetching current location by IP from %s...", _API_URL)
    req = urllib.request.Request(_API_URL)
    with urllib.request.urlopen(req, timeout=10) as resp:
        charset = resp.headers.get_content_charset() or "utf-8"
        body = resp.read().decode(charset)

    data = json.loads(body)
    if data.get("status") != "success":
        raise ValueError(f"IP-API error: {data.get('message', 'unknown error')}")
    return data
