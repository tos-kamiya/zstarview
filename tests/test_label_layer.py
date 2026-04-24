from __future__ import annotations

from PySide6.QtCore import QPointF
from PySide6.QtGui import QColor, QFont, QImage, QPainter
from PySide6.QtWidgets import QApplication

from zstarview.render import search_overlay as render_search_overlay
from zstarview.render import text as render_text
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
        visual_preset="night",
        text_font=QFont(),
        draw_marker=False,
        label_candidates=label_candidates,
    )
    painter.end()

    assert len(label_candidates) == 1
    assert label_candidates[0]["text"] == "Mars"
    assert label_candidates[0]["priority"] == 15
