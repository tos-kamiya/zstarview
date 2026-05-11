from __future__ import annotations

import pytest
from PySide6.QtCore import QPointF, QRect, QRectF
from PySide6.QtGui import QColor, QFont, QImage, QPainter
from PySide6.QtWidgets import QApplication

from zstarview.render import search_overlay as render_search_overlay
from zstarview.render import pipeline as render_pipeline
from zstarview.render import text as render_text
from zstarview.render import asterisms as render_asterisms
from zstarview.render import deep_sky_objects as render_deep_sky_objects
from zstarview.render import overlay_info as render_overlay_info
from zstarview.render import guides as render_guides
from zstarview.render import solar_system as render_solar_system
from zstarview.paths import THEME_STYLES_BY_PRESET
from zstarview.search.models import SearchJumpTarget
from zstarview.types import ScreenGeometry

_app = QApplication.instance() or QApplication([])


def test_draw_label_candidates_uses_priority_order_and_offsets(monkeypatch) -> None:
    captured: list[tuple[str, float, float]] = []

    def fake_draw_outlined_text(
        painter,
        text,
        pos,
        font=None,
        text_color=QColor(255, 255, 255),
        outline_color=QColor(0, 0, 0, 0),
        outline_width=3.0,
        *,
        style=None,
    ) -> None:
        captured.append((text, float(pos.x()), float(pos.y())))

    monkeypatch.setattr(render_text, "draw_outlined_text", fake_draw_outlined_text)

    img = QImage(240, 120, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    font = QFont()
    font.setPointSize(14)
    style = render_text.ResolvedTextStyle(
        font=font,
        text_color=QColor(255, 255, 255),
        outline_color=QColor(0, 0, 0),
        outline_width=2.0,
    )

    render_text._draw_label_candidates(
        painter,
        [
            {"text": "Primary", "pos": QPointF(40.0, 60.0), "style": style, "priority": 10},
            {"text": "Secondary", "pos": QPointF(40.0, 60.0), "style": style, "priority": 20},
        ],
        font,
    )
    painter.end()

    assert [item[0] for item in captured] == ["Primary", "Secondary"]
    assert captured[0][1:] != captured[1][1:]
    assert captured[0][1] == 40.0
    assert captured[1][1:] != (40.0, 60.0)


def test_draw_label_candidates_prefers_upper_item_when_pair_clusters(monkeypatch) -> None:
    captured: list[str] = []

    def fake_draw_outlined_text(
        painter,
        text,
        pos,
        font=None,
        text_color=QColor(255, 255, 255),
        outline_color=QColor(0, 0, 0, 0),
        outline_width=3.0,
        *,
        style=None,
    ) -> None:
        del painter, pos, font, text_color, outline_color, outline_width, style
        captured.append(text)

    def fake_text_bounds(text: str, font: QFont, baseline_pos: QPointF) -> QRectF:
        del text, font
        return QRectF(float(baseline_pos.x()), float(baseline_pos.y()) - 4.0, 40.0, 20.0)

    monkeypatch.setattr(render_text, "draw_outlined_text", fake_draw_outlined_text)
    monkeypatch.setattr(render_text, "_text_bounds_at_baseline", fake_text_bounds)

    img = QImage(240, 200, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    font = QFont()
    font.setPointSize(12)
    style = render_text.ResolvedTextStyle(
        font=font,
        text_color=QColor(255, 255, 255),
        outline_color=QColor(0, 0, 0),
        outline_width=2.0,
    )

    render_text._draw_label_candidates(
        painter,
        [
            {"text": "Lower", "pos": QPointF(80.0, 120.0), "style": style, "priority": 10},
            {"text": "Upper", "pos": QPointF(80.0, 100.0), "style": style, "priority": 10},
        ],
        font,
    )
    painter.end()

    assert captured[:2] == ["Upper", "Lower"]


def test_draw_label_candidates_explores_more_offsets_for_dense_clusters(monkeypatch) -> None:
    captured: list[tuple[str, float, float]] = []

    def fake_draw_outlined_text(
        painter,
        text,
        pos,
        font=None,
        text_color=QColor(255, 255, 255),
        outline_color=QColor(0, 0, 0, 0),
        outline_width=3.0,
        *,
        style=None,
    ) -> None:
        captured.append((text, float(pos.x()), float(pos.y())))

    def fake_text_bounds(text: str, font: QFont, baseline_pos: QPointF) -> QRectF:
        del text, font
        return QRectF(float(baseline_pos.x()), float(baseline_pos.y()), 2.0, 2.0)

    monkeypatch.setattr(render_text, "draw_outlined_text", fake_draw_outlined_text)
    monkeypatch.setattr(render_text, "_text_bounds_at_baseline", fake_text_bounds)

    img = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    font = QFont()
    font.setPointSize(12)
    style = render_text.ResolvedTextStyle(
        font=font,
        text_color=QColor(255, 255, 255),
        outline_color=QColor(0, 0, 0),
        outline_width=2.0,
    )

    render_text._draw_label_candidates(
        painter,
        [
            {
                "text": f"Label {idx}",
                "pos": QPointF(120.0, 120.0),
                "style": style,
                "priority": 10,
                "hide_on_overlap": True,
            }
            for idx in range(16)
        ],
        font,
    )
    painter.end()

    assert len(captured) == 16
    assert captured[-1][0] == "Label 15"
    assert max(abs(captured[-1][1] - 120.0), abs(captured[-1][2] - 120.0)) >= 42.0


def test_draw_search_target_overlay_appends_label_candidate() -> None:
    img = QImage(240, 120, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    label_candidates: list[dict[str, object]] = []
    target = SearchJumpTarget(
        label="Mars",
        kind="planet",
        sort_key=(0.0, "mars"),
        alt_deg=45.0,
        az_deg=120.0,
    )

    render_search_overlay.draw_search_target_overlay(
        painter,
        ScreenGeometry(center=(120, 60), radius=80),
        target,
        view_center=(45.0, 180.0),
        edge_fov_deg=95.0,
        content_fov_deg=110.0,
        text_font=QFont(),
        draw_marker=False,
        label_candidates=label_candidates,
        theme=THEME_STYLES_BY_PRESET["night"],
    )
    painter.end()

    assert len(label_candidates) == 1
    assert label_candidates[0]["text"] == "Mars"
    assert label_candidates[0]["priority"] == 15


def test_draw_search_target_overlay_keeps_label_at_edge_of_content_fov() -> None:
    img = QImage(240, 120, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    draw_calls: list[str] = []
    target = SearchJumpTarget(
        label="Mars",
        kind="planet",
        sort_key=(0.0, "mars"),
        alt_deg=45.0,
        az_deg=120.0,
    )

    def fake_is_in_fov(alt, az, view_center, *, fov_deg):
        del alt, az, view_center
        return float(fov_deg) == 95.0

    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(render_search_overlay, "is_in_fov", fake_is_in_fov)
    monkeypatch.setattr(
        render_search_overlay.render_guides,
        "draw_gauge_cross",
        lambda *_args, **_kwargs: draw_calls.append("marker"),
    )
    monkeypatch.setattr(
        render_search_overlay.render_text,
        "draw_outlined_text",
        lambda *_args, **_kwargs: draw_calls.append("label"),
    )
    try:
        render_search_overlay.draw_search_target_overlay(
            painter,
            ScreenGeometry(center=(120, 60), radius=80),
            target,
            view_center=(45.0, 180.0),
            edge_fov_deg=95.0,
            content_fov_deg=110.0,
            text_font=QFont(),
            draw_marker=True,
            draw_label=True,
            theme=THEME_STYLES_BY_PRESET["night"],
        )
    finally:
        monkeypatch.undo()
        painter.end()

    assert draw_calls == ["marker", "label"]


def test_hover_overlay_passes_label_candidates_to_asterisms(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def fake_draw_asterisms(
        painter,
        geometry,
        celestial_data,
        viewer_data,
        highlighted_object,
        text_font,
        label_reservations=None,
        label_candidates=None,
        **kwargs,
    ) -> None:
        captured["label_candidates"] = label_candidates

    monkeypatch.setattr(render_asterisms, "draw_asterisms", fake_draw_asterisms)
    monkeypatch.setattr(render_solar_system, "draw_hovered_moon_overlay", lambda *args, **kwargs: None)
    monkeypatch.setattr(render_deep_sky_objects, "draw_dso_hover_info", lambda *args, **kwargs: None)
    monkeypatch.setattr(render_guides, "draw_direction_hover_guide", lambda *args, **kwargs: None)
    monkeypatch.setattr(render_overlay_info, "draw_overlay_info", lambda *args, **kwargs: None)
    monkeypatch.setattr(render_pipeline, "_draw_dso_hover_layer", lambda *args, **kwargs: None)

    scene = type(
        "Scene",
        (),
        {
            "celestial_data": object(),
            "viewer": type(
                "Viewer",
                (),
                {"view_center": (45.0, 180.0), "edge_fov_deg": 95.0, "content_fov_deg": 110.0},
            )(),
        },
    )()
    style = type(
        "Style",
        (),
        {
            "show_asterisms": True,
            "text_font": QFont(),
            "visual_preset": "night",
            "theme": THEME_STYLES_BY_PRESET["night"],
            "bright_bodies_mode": "outline",
            "star_render_expected_width": 128,
            "show_guidelines": False,
            "vmag_limit": 6.0,
        },
    )()

    img = QImage(64, 64, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    painter = QPainter(img)
    try:
        render_pipeline._draw_hover_overlay_layer(
            painter=painter,
            geometry=ScreenGeometry(center=(32, 32), radius=32),
            viewport_rect=QRect(0, 0, 64, 64),
            scene=scene,
            style=style,
            mouse_pos=None,
            highlighted_object=({"source_id": "HIP1", "name": "Star A"}, QPointF(32.0, 32.0)),
            highlighted_dso=None,
            label_candidates=[],
        )
    finally:
        painter.end()

    assert captured["label_candidates"] == []
