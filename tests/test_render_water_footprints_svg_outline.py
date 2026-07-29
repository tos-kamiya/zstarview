from __future__ import annotations

import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


def _load_outline_module():
    module_path = Path(__file__).resolve().parents[1] / "dev-samples" / "render_water_footprints_svg_outline.py"
    spec = spec_from_file_location("render_water_footprints_svg_outline", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("failed to load render_water_footprints_svg_outline.py")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_build_svg_draws_outer_rings_as_strokes_only() -> None:
    outline = _load_outline_module()
    footprints = (
        outline.Footprint(
            index=0,
            water_id="way/1",
            kind="natural_water",
            outer_rings_lonlat=(
                ((0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)),
            ),
            inner_rings_lonlat=(
                ((0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5), (0.5, 0.5)),
            ),
            source="way",
            tags={},
        ),
    )

    svg = outline.build_svg(
        footprints,
        width=400,
        height=300,
        padding_ratio=0.0,
        background="#ffffff",
        show_labels=False,
        show_legend=False,
        ele_only=False,
    )

    assert 'fill="none"' in svg
    assert 'fill-rule="evenodd"' not in svg
    assert svg.count("<path ") == 1
    assert "0.5,0.5" not in svg


def test_load_footprints_accepts_water_polygon_payload(tmp_path) -> None:
    outline = _load_outline_module()
    payload_path = tmp_path / "footprints.json"
    payload_path.write_text(
        """{
  "water_polygons": [
    {
      "osm_id": "way/1",
      "kind": "coastline",
      "source": "coastline",
      "tags": {"water_level": "12"},
      "outer_rings": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]],
      "inner_rings": [[[0.2, 0.2], [0.8, 0.2], [0.8, 0.8], [0.2, 0.8], [0.2, 0.2]]]
    },
    {
      "osm_id": "way/2",
      "kind": "natural_water",
      "source": "way",
      "tags": {"water_level": "12"},
      "outer_rings": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]],
      "inner_rings": [[[0.2, 0.2], [0.8, 0.2], [0.8, 0.8], [0.2, 0.8], [0.2, 0.2]]]
    }
  ]
}
""",
        encoding="utf-8",
    )

    footprints = outline.load_footprints(payload_path)

    assert len(footprints) == 1
    assert footprints[0].kind == "coastline"
    assert footprints[0].has_explicit_height
    assert len(footprints[0].outer_rings_lonlat) == 1
    assert len(footprints[0].inner_rings_lonlat) == 1


def test_main_accepts_input_cache(monkeypatch, tmp_path) -> None:
    outline = _load_outline_module()
    cache_path = tmp_path / "cache.json"
    output_path = tmp_path / "out.svg"
    cache_path.write_text(
        """{
  "cache_format_version": 2,
  "query": {
    "bbox": {
      "west": 0.0,
      "south": 0.0,
      "east": 1.0,
      "north": 1.0
    }
  },
  "footprints": [
    {
      "water_id": "way/1",
      "kind": "coastline",
      "source": "coastline",
      "tags": {},
      "outer_rings_lonlat": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]],
      "inner_rings_lonlat": []
    },
    {
      "water_id": "way/2",
      "kind": "natural_water",
      "source": "way",
      "tags": {},
      "outer_rings_lonlat": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]],
      "inner_rings_lonlat": []
    }
  ]
}
""",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "render_water_footprints_svg_outline.py",
            "--input-cache",
            str(cache_path),
            "--output",
            str(output_path),
        ],
    )

    assert outline.main() == 0
    assert output_path.exists()
    assert 'stroke="#8c5a00"' in output_path.read_text(encoding="utf-8")
