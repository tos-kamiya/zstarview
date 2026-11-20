from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import create_autospec

import pytest

import zstarview.render.atlas_pipeline as atlas_pipeline_module
import zstarview.render.pipeline as pipeline_module


def test_instrument_context_layers_pass_viewer_to_urban_outline_layer(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    urban_outline_layer = create_autospec(pipeline_module._draw_urban_outline_layer)
    monkeypatch.setattr(
        pipeline_module, "_draw_urban_outline_layer", urban_outline_layer
    )
    monkeypatch.setattr(
        pipeline_module, "compute_star_render_upscale_factor", lambda *_args: 1.0
    )
    monkeypatch.setattr(
        pipeline_module, "_draw_main_terrain_profile_layer", lambda *_a, **_k: None
    )
    monkeypatch.setattr(
        pipeline_module, "_should_draw_water_overlay", lambda *_args: False
    )
    painter = object()
    geometry = SimpleNamespace(radius=10)
    scene = object()
    viewer = object()
    style = SimpleNamespace(
        star_render_expected_width=None,
        show_dso=False,
        show_asterisms=False,
        earth_guide_opacity=0.0,
    )

    atlas_pipeline_module._draw_instrument_context_layers(
        painter,
        geometry=geometry,
        scene=scene,
        viewer=viewer,
        style=style,
        label_candidates=[],
    )

    urban_outline_layer.assert_called_once_with(
        painter, geometry=geometry, scene=scene, viewer=viewer, style=style
    )
