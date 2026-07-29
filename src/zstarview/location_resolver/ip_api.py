from __future__ import annotations

import json
import logging
import time
import urllib.request
from pathlib import Path
from typing import Any

from ..paths import CACHE_PATH
from ..user_agent import build_user_agent

logger = logging.getLogger(__name__)

_IP_API_URL = "http://ip-api.com/json"
_IP_API_LAST_REQUEST_FILE = Path(CACHE_PATH) / "ip_api_last_request.timestamp"


def fetch_location_by_ip() -> dict[str, Any]:
    """Fetch current location data from ip-api.com with a 3-second rate limit."""
    now = time.time()
    if _IP_API_LAST_REQUEST_FILE.exists():
        try:
            last_request_time = float(_IP_API_LAST_REQUEST_FILE.read_text())
            if now - last_request_time < 3.0:
                raise RuntimeError("Rate limit exceeded: too many requests to ip-api.com. Please wait at least 3 seconds.")
        except (ValueError, OSError):
            _IP_API_LAST_REQUEST_FILE.unlink(missing_ok=True)

    logger.info("Fetching current location by IP from %s...", _IP_API_URL)
    req = urllib.request.Request(
        _IP_API_URL,
        headers={"User-Agent": build_user_agent("ip-api")},
    )
    with urllib.request.urlopen(req, timeout=10) as resp:
        charset = resp.headers.get_content_charset() or "utf-8"
        body = resp.read().decode(charset)

    _IP_API_LAST_REQUEST_FILE.write_text(str(now))

    data = json.loads(body)
    if data.get("status") != "success":
        raise ValueError(f"IP-API error: {data.get('message', 'unknown error')}")
    return data
