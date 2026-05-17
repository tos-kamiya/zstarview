from zstarview.location_resolver import ResolvedLocation, format_splash_location


def test_format_splash_location_includes_lat_lon_and_height() -> None:
    location = ResolvedLocation(
        display_name="t/Tokyo Skytree",
        lat=35.71006,
        lon=139.8107,
        tz="Asia/Tokyo",
        persistence_key="tower:t/Tokyo Skytree",
        observer_height_m=1.7,
        kind="tower",
        ground_elevation_m=12.3,
        location_height_m=634.0,
    )

    assert (
        format_splash_location(location)
        == (
            "Location: t/Tokyo Skytree | Lat: 35.71006, Lon: 139.81070 | "
            "Height: ground 12.3 m, building 634 m, add 1.7 m"
        )
    )


def test_format_splash_location_does_not_repeat_lat_lon_name() -> None:
    location = ResolvedLocation(
        display_name="Lat: 35.71006, Lon: 139.81070",
        lat=35.71006,
        lon=139.8107,
        tz="Asia/Tokyo",
        persistence_key="35.71006;139.81070",
        observer_height_m=1.7,
        kind="coords",
        ground_elevation_m=12.3,
        location_height_m=0.0,
    )

    assert (
        format_splash_location(location)
        == (
            "Location: Lat: 35.71006, Lon: 139.81070 | "
            "Height: ground 12.3 m, add 1.7 m"
        )
    )


def test_format_splash_location_uses_additional_height_for_tower_base() -> None:
    location = ResolvedLocation(
        display_name="t/Tokyo Skytree",
        lat=35.710055555,
        lon=139.810722222,
        tz="Asia/Tokyo",
        persistence_key="t/Tokyo Skytree",
        observer_height_m=635.7,
        kind="tower",
        ground_elevation_m=5.4,
        location_height_m=634.0,
    )

    assert (
        format_splash_location(location)
        == (
            "Location: t/Tokyo Skytree | Lat: 35.71006, Lon: 139.81072 | "
            "Height: ground 5.4 m, building 634 m, add 1.7 m"
        )
    )
