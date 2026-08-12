"""HTTP access helpers for published GMN trajectory summaries."""

from __future__ import annotations

from urllib.parse import urljoin
from urllib.request import Request, urlopen

from ..user_agent import build_user_agent
from .constants import GMN_DAILY_INDEX_URL


def fetch_gmn_text(url: str, *, timeout_s: float = 30.0) -> str:
    request = Request(
        url,
        headers={
            "Accept": "text/plain,text/html;q=0.9,*/*;q=0.1",
            "User-Agent": build_user_agent("gmn-meteors"),
        },
    )
    with urlopen(request, timeout=float(timeout_s)) as response:
        payload = response.read()
        charset = response.headers.get_content_charset() or "utf-8"
    return payload.decode(charset, errors="replace")


def gmn_daily_file_url(
    filename: str,
    *,
    index_url: str = GMN_DAILY_INDEX_URL,
) -> str:
    return urljoin(index_url, filename)
