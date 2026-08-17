
from tests._planet_marker_support import *


def test_satellite_overlay_draws_below_horizon_marker_when_in_fov(monkeypatch) -> None:
    cross_calls: list[tuple[float, float]] = []

    monkeypatch.setattr(
        render_satellites,
        "draw_gauge_cross",
        lambda _painter, _color, _center, *, scale=1.0, pen_width=1.0: (
            cross_calls.append((scale, pen_width))
        ),
    )
    monkeypatch.setattr(
        render_satellites,
        "project_satellite_records",
        lambda *_args, **_kwargs: [
            SatelliteOverlayPoint(
                group_key="iss",
                satellite_name="ISS (ZARYA)",
                alt_deg=-40.0,
                az_deg=151.0,
                marker_scale=0.42,
            )
        ],
    )

    image = QImage(40, 40, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(20, 20), radius=20),
            satellite_records_by_group={"iss": []},
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(0.0, 151.0),
                observer_height_m=1.7,
            ),
            time_obj=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
            opacity=1.0,
        )
    finally:
        painter.end()

    assert cross_calls == [(0.42, 1.0)]


def test_atlas_satellite_cross_uses_aircraft_weighted_line(monkeypatch) -> None:
    cross_calls: list[float] = []

    monkeypatch.setattr(
        render_satellites,
        "draw_gauge_cross",
        lambda _painter, _color, _center, *, scale=1.0, pen_width=1.0: (
            cross_calls.append(pen_width)
        ),
    )
    monkeypatch.setattr(
        render_satellites,
        "project_satellite_records",
        lambda *_args, **_kwargs: [
            SatelliteOverlayPoint(
                group_key="iss",
                satellite_name="ISS (ZARYA)",
                alt_deg=-40.0,
                az_deg=151.0,
                marker_scale=0.42,
            )
        ],
    )

    image = QImage(40, 40, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(20, 20), radius=20),
            satellite_records_by_group={"iss": []},
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(0.0, 151.0),
                observer_height_m=1.7,
            ),
            time_obj=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
            opacity=1.0,
            theme=THEME_STYLES_BY_PRESET[ATLAS_THEME_PRESET],
        )
    finally:
        painter.end()

    assert cross_calls == [4.0, 2.0]


def test_satellite_overlay_scales_marker_with_window_scale(monkeypatch) -> None:
    cross_calls: list[tuple[float, float]] = []

    monkeypatch.setattr(
        render_satellites,
        "draw_gauge_cross",
        lambda _painter, _color, _center, *, scale=1.0, pen_width=1.0: (
            cross_calls.append((scale, pen_width))
        ),
    )
    monkeypatch.setattr(
        render_satellites,
        "project_satellite_records",
        lambda *_args, **_kwargs: [
            SatelliteOverlayPoint(
                group_key="iss",
                satellite_name="ISS (ZARYA)",
                alt_deg=-40.0,
                az_deg=151.0,
                marker_scale=0.42,
            )
        ],
    )

    image = QImage(40, 40, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(20, 20), radius=20),
            satellite_records_by_group={"iss": []},
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(0.0, 151.0),
                observer_height_m=1.7,
            ),
            time_obj=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
            opacity=1.0,
            marker_scale=2.5,
        )
    finally:
        painter.end()

    assert cross_calls == [(1.05, 1.0)]


def test_satellite_overlay_keeps_overscan_position_beyond_90_deg(monkeypatch) -> None:
    positions: list[tuple[float, float]] = []

    def _record_cross(_painter, _color, center, *, scale=1.0, pen_width=1.0) -> None:
        positions.append((float(center.x()), float(center.y())))

    monkeypatch.setattr(render_satellites, "draw_gauge_cross", _record_cross)
    monkeypatch.setattr(
        render_satellites,
        "project_satellite_records",
        lambda *_args, **_kwargs: [
            SatelliteOverlayPoint(
                group_key="iss",
                satellite_name="ISS (ZARYA)",
                alt_deg=-50.0,
                az_deg=180.0,
                marker_scale=0.42,
            )
        ],
    )

    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(120, 120), radius=80),
            satellite_records_by_group={"iss": []},
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(45.0, 180.0),
                edge_fov_deg=110.0,
                content_fov_deg=110.0,
                observer_height_m=1.7,
            ),
            time_obj=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
            opacity=1.0,
        )
    finally:
        painter.end()

    assert len(positions) == 1
    x, y = positions[0]
    assert math.isclose(x, 120.0, abs_tol=0.2)
    assert math.isclose(y, 189.1, abs_tol=0.3)


def test_satellite_overlay_does_not_hide_old_element_epoch(monkeypatch) -> None:
    cross_calls: list[tuple[float, float]] = []

    monkeypatch.setattr(
        render_satellites,
        "draw_gauge_cross",
        lambda _painter, _color, _center, *, scale=1.0, pen_width=1.0: (
            cross_calls.append((scale, pen_width))
        ),
    )

    image = QImage(40, 40, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(20, 20), radius=20),
            satellite_records_by_group={
                SATELLITE_HORIZONS_CACHE_KEY: [
                    {
                        "OBJECT_NAME": "JWST",
                        "ALT_DEG": 10.0,
                        "AZ_DEG": 151.0,
                    }
                ]
            },
            viewer_data=ViewerData(
                location=(0.0, 0.0),
                timezone_name="UTC",
                city_name="Test",
                view_center=(10.0, 151.0),
                edge_fov_deg=110.0,
                content_fov_deg=110.0,
                observer_height_m=0.0,
            ),
            time_obj=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
            opacity=1.0,
        )
    finally:
        painter.end()

    assert cross_calls == [(0.3, 1.0)]


app = QApplication.instance() or QApplication([])


def test_satellite_overlay_does_not_add_labels(monkeypatch) -> None:
    monkeypatch.setattr(
        render_satellites, "draw_gauge_cross", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        render_satellites,
        "project_satellite_records",
        lambda *_args, **_kwargs: [
            SatelliteOverlayPoint(
                group_key="iss",
                satellite_name="ISS (ZARYA)",
                alt_deg=10.0,
                az_deg=151.0,
                marker_scale=0.42,
            )
        ],
    )

    image = QImage(40, 40, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(20, 20), radius=20),
            satellite_records_by_group={"iss": []},
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(0.0, 151.0),
                observer_height_m=1.7,
            ),
            time_obj=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
            opacity=1.0,
        )
    finally:
        painter.end()


def test_satellite_overlay_draws_simplified_labels(monkeypatch) -> None:
    label_calls: list[tuple[str, QColor, float]] = []

    monkeypatch.setattr(
        render_satellites,
        "draw_gauge_cross",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        render_satellites,
        "project_satellite_records",
        lambda *_args, **_kwargs: [
            SatelliteOverlayPoint(
                group_key="iss",
                satellite_name="ISS (ZARYA)",
                alt_deg=10.0,
                az_deg=151.0,
                marker_scale=0.42,
            )
        ],
    )

    def fake_draw_outlined_text(_painter, text, *_args, **kwargs) -> None:
        style = kwargs["style"]
        label_calls.append((str(text), style.text_color, float(style.outline_width)))

    monkeypatch.setattr(render_text, "draw_outlined_text", fake_draw_outlined_text)

    image = QImage(40, 40, QImage.Format.Format_ARGB32_Premultiplied)
    painter = QPainter(image)
    try:
        render_satellites.draw_satellite_overlay(
            painter=painter,
            geometry=ScreenGeometry(center=(20, 20), radius=20),
            satellite_records_by_group={"iss": []},
            viewer_data=ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(0.0, 151.0),
                observer_height_m=1.7,
            ),
            time_obj=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
            opacity=1.0,
            draw_simplified_labels=True,
            text_font=QFont(),
        )
    finally:
        painter.end()

    assert len(label_calls) == 1
    text, color, outline_width = label_calls[0]
    assert text == "ISS (ZARYA)"
    assert (
        color.red(),
        color.green(),
        color.blue(),
    ) == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
    assert color.alpha() == round(255 * 0.7)
    assert (
        outline_width
        == render_text.resolve_text_style(
            THEME_STYLES_BY_PRESET["night"],
            QFont(),
        ).outline_width
    )


def test_satellite_overlay_info_shows_hover_name(monkeypatch) -> None:
    class DummyPainter:
        def setPen(self, *_args, **_kwargs) -> None:
            pass

        def setBrush(self, *_args, **_kwargs) -> None:
            pass

        def drawEllipse(self, *_args, **_kwargs) -> None:
            raise AssertionError("satellite hover should not draw a circle")

    painter = DummyPainter()
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=12.0,
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    hovered_satellite = (
        SatelliteOverlayPoint(
            group_key="horizons",
            satellite_name="JWST",
            alt_deg=12.0,
            az_deg=220.0,
            marker_scale=0.3,
        ),
        QPointF(120.0, 90.0),
    )
    label_calls: list[tuple[str, QColor]] = []

    def fake_draw_outlined_text(_painter, text, *_args, **kwargs) -> None:
        style = kwargs["style"]
        label_calls.append((str(text), style.text_color))

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        highlighted_dso=None,
        highlighted_object=None,
        highlighted_satellite=hovered_satellite,
        text_font=QFont(),
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert any(
        text == "JWST"
        and (color.red(), color.green(), color.blue())
        == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
        for text, color in label_calls
    )


def test_overlay_info_skips_hover_satellite_label_when_simplified_labels_are_drawn(
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
        city_name="Tokyo",
        view_center=(45.0, 180.0),
        observer_height_m=12.0,
    )
    geometry = ScreenGeometry(center=(120, 90), radius=70)
    hovered_satellite = (
        SatelliteOverlayPoint(
            group_key="horizons",
            satellite_name="JWST",
            alt_deg=12.0,
            az_deg=220.0,
            marker_scale=0.3,
        ),
        QPointF(120.0, 90.0),
    )
    label_calls: list[tuple[str, QColor]] = []

    def fake_draw_outlined_text(_painter, text, *_args, **kwargs) -> None:
        style = kwargs["style"]
        label_calls.append((str(text), style.text_color))

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        highlighted_dso=None,
        highlighted_object=None,
        highlighted_satellite=hovered_satellite,
        text_font=QFont(),
        draw_simplified_satellite_labels=True,
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert all(text != "JWST" for text, _color in label_calls)


def test_overlay_info_shows_star_and_satellite_labels_independently() -> None:
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
    highlighted_object = ({"name": "Vega"}, QPointF(100.0, 80.0))
    highlighted_satellite = (
        SatelliteOverlayPoint(
            group_key="iss",
            satellite_name="ISS",
            alt_deg=12.0,
            az_deg=220.0,
            marker_scale=0.42,
        ),
        QPointF(130.0, 95.0),
    )
    label_calls: list[tuple[str, QColor]] = []

    def fake_draw_outlined_text(_painter, text, *_args, **kwargs) -> None:
        style = kwargs["style"]
        label_calls.append((str(text), style.text_color))

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        highlighted_dso=None,
        highlighted_object=highlighted_object,
        highlighted_satellite=highlighted_satellite,
        text_font=QFont(),
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    label_names = {text for text, _color in label_calls}
    assert "Vega" in label_names
    assert "ISS" in label_names
    assert any(
        text == "ISS"
        and (color.red(), color.green(), color.blue())
        == PALETTE_AIRCRAFT_AND_SATELLITE_RGB
        for text, color in label_calls
    )


def test_overlay_info_colors_dso_labels_like_the_dso_marker() -> None:
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
    highlighted_dso = (
        {
            "name": "Andromeda Galaxy",
            "alt": 45.0,
            "az": 180.0,
            "major_arcmin": 30.0,
            "minor_arcmin": 20.0,
            "pa_deg": 0.0,
        },
        QPointF(120.0, 90.0),
    )
    label_calls: list[tuple[str, QColor]] = []

    def fake_draw_outlined_text(_painter, text, *_args, **kwargs) -> None:
        style = kwargs["style"]
        label_calls.append((str(text), style.text_color))

    render_overlay_info.draw_overlay_info(
        painter,
        geometry,
        _empty_celestial_data([]),
        viewer,
        vmag_limit=6.0,
        highlighted_dso=highlighted_dso,
        highlighted_object=None,
        highlighted_satellite=None,
        text_font=QFont(),
        draw_outlined_text_func=fake_draw_outlined_text,
        text_bounds_at_baseline_func=render_text._text_bounds_at_baseline,
        theme=THEME_STYLES_BY_PRESET["night"],
    )

    assert any(
        text == "Andromeda Galaxy"
        and (color.red(), color.green(), color.blue()) == DSO_LABEL_TEXT_RGB
        for text, color in label_calls
    )


def test_find_highlighted_satellite_prefers_nearest_marker(monkeypatch) -> None:
    monkeypatch.setattr(render_satellites, "is_in_fov", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(
        render_satellites,
        "altaz_to_normalized_xy",
        lambda alt, az, _view_center, **_kwargs: (float(alt), float(az)),
    )
    monkeypatch.setattr(
        render_satellites,
        "normalized_to_screen_xy",
        lambda nx, ny, _geometry: (float(nx), float(ny)),
    )

    geometry = ScreenGeometry(center=(120, 90), radius=70)
    view_center = (45.0, 180.0)
    satellite_points = [
        SatelliteOverlayPoint(
            group_key="iss",
            satellite_name="ISS",
            alt_deg=120.0,
            az_deg=90.0,
            marker_scale=0.42,
        ),
        SatelliteOverlayPoint(
            group_key="horizons",
            satellite_name="JWST",
            alt_deg=123.0,
            az_deg=90.0,
            marker_scale=0.30,
        ),
    ]
    monkeypatch.setattr(
        render_satellites,
        "project_satellite_records",
        lambda *_args, **_kwargs: satellite_points,
    )

    highlighted = render_satellites.find_highlighted_satellite(
        satellite_records_by_group={"iss": []},
        mouse_pos=QPointF(121.0, 90.0),
        geometry=geometry,
        viewer_data=ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=view_center,
            observer_height_m=1.7,
        ),
        time_obj=astropy.time.Time("2026-02-27T00:00:00", scale="utc"),
    )

    assert highlighted is not None
    satellite, _ = highlighted
    assert satellite.satellite_name == "ISS"
