from __future__ import annotations

import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


def _load_module():
    module_path = Path(__file__).resolve().parents[1] / "dev-samples" / "render_coastline_inner_rings_svg.py"
    spec = spec_from_file_location("render_coastline_inner_rings_svg", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("failed to load render_coastline_inner_rings_svg.py")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_build_svg_draws_inner_rings_filled_only() -> None:
    module = _load_module()
    footprints = (
        module.Footprint(  # noqa: SLF001
            index=0,
            water_id="coastline/0",
            kind="coastline",
            outer_rings_lonlat=(
                ((0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)),
            ),
            inner_rings_lonlat=(
                ((0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5), (0.5, 0.5)),
            ),
            source="coastline",
            tags={},
        ),
    )

    svg = module.build_svg(  # noqa: SLF001
        footprints,
        width=400,
        height=300,
        padding_ratio=0.0,
        background="#ffffff",
        show_labels=False,
        show_legend=False,
        ele_only=False,
    )

    assert 'fill="none"' not in svg
    assert 'fill="#4c8cff"' in svg
    assert '0.0,0.0' not in svg
    assert svg.count("<path ") == 1


def test_load_footprints_keeps_coastline_only(tmp_path) -> None:
    module = _load_module()
    payload_path = tmp_path / "cache.json"
    payload_path.write_text(
        """{
  "water_polygons": [
    {
      "osm_id": "coastline/0",
      "kind": "coastline",
      "source": "coastline",
      "tags": {},
      "outer_rings": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]],
      "inner_rings": [[[0.2, 0.2], [0.8, 0.2], [0.8, 0.8], [0.2, 0.8], [0.2, 0.2]]]
    },
    {
      "osm_id": "way/2",
      "kind": "natural_water",
      "source": "way",
      "tags": {},
      "outer_rings": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]],
      "inner_rings": []
    }
  ]
}
""",
        encoding="utf-8",
    )

    footprints = module.load_footprints(payload_path)  # noqa: SLF001

    assert len(footprints) == 1
    assert footprints[0].kind == "coastline"
    assert len(footprints[0].inner_rings_lonlat) == 1
