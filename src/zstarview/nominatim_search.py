from .location_resolver.nominatim import search_nominatim


def search(
    query: str,
    *,
    limit: int = 5,
    countrycode: str | None = None,
    language: str = "en",
    user_agent: str | None = None,
) -> list[dict]:
    kwargs = {
        "limit": limit,
        "countrycode": countrycode,
        "language": language,
    }
    if user_agent is not None:
        kwargs["user_agent"] = user_agent
    return search_nominatim(query, **kwargs)


__all__ = [
    "search",
]
