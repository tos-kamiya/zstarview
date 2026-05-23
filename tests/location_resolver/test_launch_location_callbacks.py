from __future__ import annotations

from zstarview.location_resolver.resolve import resolve_launch_location


def test_resolve_launch_location_uses_injected_callbacks() -> None:
    saved_values: list[str | dict[str, object]] = []

    def save_last_city(value: str | dict[str, object]) -> None:
        saved_values.append(value)

    location = resolve_launch_location(
        None,
        load_last_city_func=lambda: "Tokyo",
        save_last_city_func=save_last_city,
    )

    assert location.display_name
    assert saved_values
    assert saved_values[0]
