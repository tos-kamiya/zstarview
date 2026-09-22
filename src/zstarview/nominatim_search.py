from .location_resolver.nominatim import search_nominatim


def search(
    query: str,
    *,
    limit: int = 5,
    countrycode: str | None = None,
    language: str = "en",
    user_agent: str | None = None,
) -> list[dict]:
    kwargs: dict[str, object] = {
        "limit": limit, "countrycode": countrycode, "language": language
    }
    if user_agent is not None:
        kwargs["user_agent"] = user_agent
    if user_agent is None:
        return search_nominatim(
            query, limit=limit, countrycode=countrycode, language=language
        )
    return search_nominatim(
        query,
        limit=limit,
        countrycode=countrycode,
        language=language,
        user_agent=user_agent,
    )


__all__ = [
    "search",
]
