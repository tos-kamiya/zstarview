
from tests._planet_marker_support import *


def test_aircraft_label_uses_black_theme_style_in_day_theme(monkeypatch) -> None:
    class _Painter:
        def save(self) -> None:
            pass

        def restore(self) -> None:
            pass

        def font(self) -> QFont:
            return QFont()

        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawPolyline(self, *_args, **_kwargs) -> None:
            pass

        def drawPolygon(self, *_args, **_kwargs) -> None:
            pass

    label_candidates: list[dict[str, object]] = []
    aircraft_points = [
        AircraftOverlayPoint(
            icao24="abc123",
            callsign="TEST123",
            alt_deg=10.0,
            az_deg=151.0,
            trail_alt_az_points=((10.0, 151.0), (10.2, 151.3)),
            distance_km=5.0,
            age_seconds=10.0,
            alpha_scale=1.0,
            trail_geodetic_points=((35.0, 139.0, 1000.0), (35.0, 139.01, 1000.0)),
        )
    ]
    monkeypatch.setattr(
        render_aircraft,
        "project_aircraft_snapshots",
        lambda *_args, **_kwargs: aircraft_points,
    )
    render_aircraft.draw_aircraft_overlay(
        painter=_Painter(),
        geometry=ScreenGeometry(center=(20, 20), radius=20),
        aircraft_snapshots={"dummy": True},
        viewer_data=ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(0.0, 151.0),
            edge_fov_deg=110.0,
            content_fov_deg=180.0,
            observer_height_m=1.7,
        ),
        time_obj=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
        opacity=1.0,
        label_candidates=label_candidates,
        theme=THEME_STYLES_BY_PRESET["day"],
    )

    assert len(label_candidates) == 1
    style = label_candidates[0]["style"]
    expected_rgb = PALETTE_AIRCRAFT_AND_SATELLITE_RGB
    assert (
        style.text_color.red(),
        style.text_color.green(),
        style.text_color.blue(),
    ) == expected_rgb
    assert (
        style.outline_width
        == render_text.resolve_text_style(
            THEME_STYLES_BY_PRESET["night"], QFont()
        ).outline_width
    )


def test_hover_ignores_unnamed_stars() -> None:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    mouse_pos = QPoint(120, 90)

    stars, meta = _star_table(names=["", "Sirius"])
    highlighted = render_stars.find_highlighted_object(
        _celestial_data_with_stars(stars, meta),
        viewer,
        mouse_pos,
        geometry,
    )

    assert highlighted is not None
    obj, _ = highlighted
    assert obj.get("name") == "Sirius"


def test_hover_returns_none_without_named_star() -> None:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    mouse_pos = QPoint(120, 90)

    stars, meta = _star_table(names=["", ""])
    highlighted = render_stars.find_highlighted_object(
        _celestial_data_with_stars(stars, meta),
        viewer,
        mouse_pos,
        geometry,
    )

    assert highlighted is None


def test_overlay_skips_label_for_planet(monkeypatch) -> None:
    class DummyPainter:
        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawEllipse(self, *_args, **_kwargs) -> None:
            pass

    painter = DummyPainter()
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    planet = PlanetBody(name="mars", alt=45.0, az=180.0, symbol="♂", is_visible=True)
    highlighted_object = (planet, QPointF(120.0, 90.0))

    label_calls: list[str] = []

    def fake_draw_outlined_text(_painter, text, *_args, **_kwargs) -> None:
        label_calls.append(text)

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=highlighted_object,
        text_font=QFont(),
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert "mars" not in label_calls


def test_overlay_info_includes_location_height_and_explicit_observer_height(
    monkeypatch,
) -> None:
    class DummyPainter:
        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawEllipse(self, *_args, **_kwargs) -> None:
            pass

    painter = DummyPainter()
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="t/Tokyo Skytree",
        view_center=(45.0, 180.0),
        observer_height_m=12.0,
        height_add_m=12.0,
        ground_elevation_m=35.0,
        location_height_label="Tower height",
        location_height_m=634.0,
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)

    label_calls: list[str] = []

    def fake_draw_outlined_text(_painter, text, *_args, **_kwargs) -> None:
        label_calls.append(text)

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=None,
        text_font=QFont(),
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert label_calls == [
        "t/Tokyo Skytree",
        "Lat: 35.00000, Lon: 139.00000 | Height: ground 35 m, building 634 m, add 12 m",
        "2026-02-27 00:00:00 UTC",
        "Alt 45°  Az 180° (S)",
    ]


def test_overlay_info_first_line_top_margin_matches_left_margin_when_cursor_is_lower_half(
    monkeypatch,
) -> None:
    class DummyPainter:
        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawEllipse(self, *_args, **_kwargs) -> None:
            pass

    painter = DummyPainter()
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="t/Tokyo Skytree",
        view_center=(45.0, 180.0),
        observer_height_m=12.0,
        height_add_m=12.0,
        ground_elevation_m=35.0,
        location_height_label="Tower height",
        location_height_m=634.0,
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    text_font = QFont()

    first_label_pos = None

    def fake_draw_outlined_text(_painter, text, pos, *_args, **_kwargs) -> None:
        nonlocal first_label_pos
        if text == "t/Tokyo Skytree" and first_label_pos is None:
            first_label_pos = pos

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=None,
        text_font=text_font,
        viewport_rect=QRectF(0.0, 0.0, 240.0, 180.0),
        mouse_pos=QPoint(10, 170),
        bottom_left=False,
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert first_label_pos is not None
    fm = QFontMetrics(text_font)
    bounds = fm.tightBoundingRect("Ag")
    top_margin = float(first_label_pos.y()) + float(bounds.top())
    left_margin = float(fm.lineSpacing())
    assert abs(top_margin - left_margin) <= 1.0


def test_overlay_info_moves_to_bottom_when_cursor_is_in_upper_half(monkeypatch) -> None:
    class DummyPainter:
        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawEllipse(self, *_args, **_kwargs) -> None:
            pass

    painter = DummyPainter()
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="t/Tokyo Skytree",
        view_center=(45.0, 180.0),
        observer_height_m=12.0,
        height_add_m=12.0,
        ground_elevation_m=35.0,
        location_height_label="Tower height",
        location_height_m=634.0,
    )

    first_label_pos = None

    def fake_draw_outlined_text(_painter, text, pos, *_args, **_kwargs) -> None:
        nonlocal first_label_pos
        if text == "t/Tokyo Skytree" and first_label_pos is None:
            first_label_pos = pos

    render_overlay_info.draw_overlay_info(
        painter,
        ScreenGeometry(center=(120, 90), radius=70),
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        enlarge_moon=False,
        highlighted_dso=None,
        highlighted_object=None,
        text_font=QFont(),
        viewport_rect=QRectF(0.0, 0.0, 240.0, 180.0),
        mouse_pos=QPoint(10, 20),
        bottom_left=True,
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert first_label_pos is not None
    fm = QFontMetrics(QFont())
    bounds = fm.tightBoundingRect("Ag")
    line_height = int(fm.lineSpacing() * 1.2)
    bottom_margin = 180.0 - (
        float(first_label_pos.y()) + float(bounds.bottom()) + 3.0 * line_height
    )
    left_margin = float(fm.lineSpacing())
    assert abs(bottom_margin - left_margin) <= 2.0


def test_format_overlay_info_lines_matches_static_overlay_order() -> None:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="t/Tokyo Skytree",
        view_center=(45.0, 180.0),
        observer_height_m=12.0,
        height_add_m=12.0,
        ground_elevation_m=35.0,
        location_height_label="Tower height",
        location_height_m=634.0,
    )

    assert render_background.format_overlay_info_lines(
        _empty_celestial_data([]), viewer, 6.0
    ) == [
        "t/Tokyo Skytree",
        "Lat: 35.00000, Lon: 139.00000 | Height: ground 35 m, building 634 m, add 12 m",
        "2026-02-27 00:00:00 UTC",
        "Alt 45°  Az 180° (S)",
    ]


def test_format_overlay_info_lines_uses_structure_as_base_for_tower_viewpoint() -> None:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="t/Tokyo Skytree",
        view_center=(45.0, 180.0),
        observer_height_m=635.7,
        height_add_m=1.7,
        ground_elevation_m=5.4,
        location_height_label="Tower height",
        location_height_m=634.0,
    )

    assert render_background.format_overlay_info_lines(
        _empty_celestial_data([]), viewer, 6.0
    ) == [
        "t/Tokyo Skytree",
        "Lat: 35.00000, Lon: 139.00000 | Height: ground 5.4 m, building 634 m, add 1.7 m",
        "2026-02-27 00:00:00 UTC",
        "Alt 45°  Az 180° (S)",
    ]
