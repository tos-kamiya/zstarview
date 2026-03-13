from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Iterable


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "src" / "zstarview" / "data" / "build_plateau_derived_tiles.py"
    spec = importlib.util.spec_from_file_location("build_plateau_derived_tiles", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_convert_tile_emits_expected_schema(tmp_path: Path) -> None:
    mod = _load_module()
    city_dir = tmp_path / "13100_tokyo23" / "udx" / "bldg"
    city_dir.mkdir(parents=True)
    path = city_dir / "53393599_bldg_6697_2_op.gml"
    path.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<core:CityModel xmlns:core="http://www.opengis.net/citygml/2.0"
                xmlns:bldg="http://www.opengis.net/citygml/building/2.0"
                xmlns:gml="http://www.opengis.net/gml">
  <gml:boundedBy>
    <gml:Envelope srsName="http://www.opengis.net/def/crs/EPSG/0/6697" srsDimension="3">
      <gml:lowerCorner>35.0 139.0 0</gml:lowerCorner>
      <gml:upperCorner>35.01 139.01 30</gml:upperCorner>
    </gml:Envelope>
  </gml:boundedBy>
  <core:cityObjectMember>
    <bldg:Building gml:id="bldg-high">
      <bldg:measuredHeight>80</bldg:measuredHeight>
      <bldg:lod0RoofEdge>
        <gml:MultiSurface>
          <gml:surfaceMember>
            <gml:Polygon><gml:exterior><gml:LinearRing>
              <gml:posList>35.0 139.0 0 35.0 139.001 0 35.001 139.001 0 35.001 139.0 0 35.0 139.0 0</gml:posList>
            </gml:LinearRing></gml:exterior></gml:Polygon>
          </gml:surfaceMember>
        </gml:MultiSurface>
      </bldg:lod0RoofEdge>
    </bldg:Building>
  </core:cityObjectMember>
</core:CityModel>
""",
        encoding="utf-8",
    )

    payload = mod.convert_tile(
        path,
        min_building_height_m=40.0,
        city_code="13100",
        source_root=city_dir,
    )

    assert payload is not None
    assert payload["schema_version"] == 1
    assert payload["source"]["city_code"] == "13100"
    assert payload["tile"]["id"] == "53393599_bldg_6697_2_op"
    buildings = payload["buildings"]
    assert len(buildings) == 1
    assert buildings[0]["height_source"] == "measuredHeight"
    assert buildings[0]["height_m"] == 80.0


def test_convert_tile_falls_back_to_lod0_footprint_when_roof_edge_is_missing(tmp_path: Path) -> None:
    mod = _load_module()
    city_dir = tmp_path / "23100_nagoya" / "udx" / "bldg"
    city_dir.mkdir(parents=True)
    path = city_dir / "52376021_bldg_6697_op.gml"
    path.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<core:CityModel xmlns:core="http://www.opengis.net/citygml/2.0"
                xmlns:bldg="http://www.opengis.net/citygml/building/2.0"
                xmlns:gml="http://www.opengis.net/gml">
  <gml:boundedBy>
    <gml:Envelope srsName="http://www.opengis.net/def/crs/EPSG/0/6697" srsDimension="3">
      <gml:lowerCorner>35.0 136.0 0</gml:lowerCorner>
      <gml:upperCorner>35.01 136.01 80</gml:upperCorner>
    </gml:Envelope>
  </gml:boundedBy>
  <core:cityObjectMember>
    <bldg:Building gml:id="bldg-high">
      <bldg:measuredHeight>66.1</bldg:measuredHeight>
      <bldg:lod0FootPrint>
        <gml:MultiSurface>
          <gml:surfaceMember>
            <gml:Polygon><gml:exterior><gml:LinearRing>
              <gml:posList>35.0 136.0 0 35.0 136.001 0 35.001 136.001 0 35.001 136.0 0 35.0 136.0 0</gml:posList>
            </gml:LinearRing></gml:exterior></gml:Polygon>
          </gml:surfaceMember>
        </gml:MultiSurface>
      </bldg:lod0FootPrint>
    </bldg:Building>
  </core:cityObjectMember>
</core:CityModel>
""",
        encoding="utf-8",
    )

    payload = mod.convert_tile(
        path,
        min_building_height_m=40.0,
        city_code="23100",
        source_root=city_dir,
    )

    assert payload is not None
    assert payload["source"]["city_code"] == "23100"
    assert payload["buildings"][0]["height_m"] == 66.1


def test_main_writes_json_file(tmp_path: Path) -> None:
    mod = _load_module()
    city_dir = tmp_path / "27100_osaka" / "udx" / "bldg"
    city_dir.mkdir(parents=True)
    path = city_dir / "tile.gml"
    path.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<core:CityModel xmlns:core="http://www.opengis.net/citygml/2.0"
                xmlns:bldg="http://www.opengis.net/citygml/building/2.0"
                xmlns:gml="http://www.opengis.net/gml">
  <gml:boundedBy>
    <gml:Envelope srsName="http://www.opengis.net/def/crs/EPSG/0/6697" srsDimension="3">
      <gml:lowerCorner>35.0 139.0 0</gml:lowerCorner>
      <gml:upperCorner>35.01 139.01 30</gml:upperCorner>
    </gml:Envelope>
  </gml:boundedBy>
  <core:cityObjectMember>
    <bldg:Building gml:id="estimated">
      <bldg:measuredHeight>-9999</bldg:measuredHeight>
      <bldg:storeysAboveGround>20</bldg:storeysAboveGround>
      <bldg:lod0RoofEdge>
        <gml:MultiSurface>
          <gml:surfaceMember>
            <gml:Polygon><gml:exterior><gml:LinearRing>
              <gml:posList>35.0 139.0 0 35.0 139.001 0 35.001 139.001 0 35.001 139.0 0 35.0 139.0 0</gml:posList>
            </gml:LinearRing></gml:exterior></gml:Polygon>
          </gml:surfaceMember>
        </gml:MultiSurface>
      </bldg:lod0RoofEdge>
    </bldg:Building>
  </core:cityObjectMember>
</core:CityModel>
""",
        encoding="utf-8",
    )
    output_dir = tmp_path / "derived"

    rc = mod.main(
        [
            "--citygml-dir",
            str(city_dir),
            "--output-dir",
            str(output_dir),
        ]
    )

    assert rc == 0
    written = output_dir / "tile.json"
    payload = json.loads(written.read_text(encoding="utf-8"))
    assert payload["buildings"][0]["height_source"] == "storeysAboveGround*3.5"
    assert payload["buildings"][0]["height_m"] == 70.0


def test_main_uses_process_pool_when_workers_gt_one(tmp_path: Path, monkeypatch) -> None:
    mod = _load_module()
    city_dir = tmp_path / "13100_tokyo23" / "udx" / "bldg"
    city_dir.mkdir(parents=True)
    (city_dir / "tile-a.gml").write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<core:CityModel xmlns:core="http://www.opengis.net/citygml/2.0"
                xmlns:bldg="http://www.opengis.net/citygml/building/2.0"
                xmlns:gml="http://www.opengis.net/gml">
  <gml:boundedBy>
    <gml:Envelope srsName="http://www.opengis.net/def/crs/EPSG/0/6697" srsDimension="3">
      <gml:lowerCorner>35.0 139.0 0</gml:lowerCorner>
      <gml:upperCorner>35.01 139.01 30</gml:upperCorner>
    </gml:Envelope>
  </gml:boundedBy>
  <core:cityObjectMember>
    <bldg:Building gml:id="bldg-high">
      <bldg:measuredHeight>80</bldg:measuredHeight>
      <bldg:lod0RoofEdge>
        <gml:MultiSurface>
          <gml:surfaceMember>
            <gml:Polygon><gml:exterior><gml:LinearRing>
              <gml:posList>35.0 139.0 0 35.0 139.001 0 35.001 139.001 0 35.001 139.0 0 35.0 139.0 0</gml:posList>
            </gml:LinearRing></gml:exterior></gml:Polygon>
          </gml:surfaceMember>
        </gml:MultiSurface>
      </bldg:lod0RoofEdge>
    </bldg:Building>
  </core:cityObjectMember>
</core:CityModel>
""",
        encoding="utf-8",
    )
    output_dir = tmp_path / "derived"
    used_workers: list[int] = []

    class FakeExecutor:
        def __init__(self, *, max_workers: int) -> None:
            used_workers.append(max_workers)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

        def map(
            self,
            fn,
            path_iter: Iterable[str],
            min_height_iter: Iterable[float],
            city_code_iter: Iterable[str | None],
            source_root_iter: Iterable[str],
        ):
            return [
                fn(path_str, min_height, city_code, source_root)
                for path_str, min_height, city_code, source_root in zip(
                    path_iter,
                    min_height_iter,
                    city_code_iter,
                    source_root_iter,
                    strict=True,
                )
            ]

    monkeypatch.setattr(mod, "ProcessPoolExecutor", FakeExecutor)

    rc = mod.main(
        [
            "--citygml-dir",
            str(city_dir),
            "--output-dir",
            str(output_dir),
            "--workers",
            "2",
        ]
    )

    assert rc == 0
    assert used_workers == [2]
    assert (output_dir / "tile-a.json").exists()
