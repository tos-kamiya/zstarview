from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
import zipfile


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "src" / "zstarview" / "data" / "import_plateau_citygml_zip.py"
    spec = importlib.util.spec_from_file_location("import_plateau_citygml_zip", mod_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_test_zip(path: Path, *, with_top_level_dir: bool) -> None:
    gml = """<?xml version="1.0" encoding="UTF-8"?>
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
"""
    prefix = "12100_chiba-shi_city_2024_citygml_1_op/" if with_top_level_dir else ""
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(prefix + "udx/bldg/tile-a.gml", gml)


def test_detect_output_city_name_normalizes_zip_filename() -> None:
    mod = _load_module()

    assert mod.detect_output_city_name(Path("12100_chiba-shi_city_2024_citygml_1_op.zip")) == "12100_chiba"
    assert mod.detect_output_city_name(Path("13100_tokyo23-ku_2022_citygml_1_2_op.zip")) == "13100_tokyo23"


def test_main_imports_zip_with_top_level_directory(tmp_path: Path) -> None:
    mod = _load_module()
    zip_path = tmp_path / "12100_chiba-shi_city_2024_citygml_1_op.zip"
    _write_test_zip(zip_path, with_top_level_dir=True)
    derived_root = tmp_path / "derived-root"

    rc = mod.main(
        [
            str(zip_path),
            "--derived-root-dir",
            str(derived_root),
            "--workers",
            "1",
        ]
    )

    assert rc == 0
    derived_dir = derived_root / "12100_chiba" / "bldg"
    payload = json.loads((derived_dir / "tile-a.json").read_text(encoding="utf-8"))
    index_payload = json.loads((derived_dir / "tile_index.json").read_text(encoding="utf-8"))
    assert payload["source"]["city_code"] == "12100"
    assert index_payload["city_name"] == "chiba"
    assert index_payload["tile_count"] == 1


def test_main_imports_zip_without_top_level_directory_and_can_write_outline_json(
    tmp_path: Path,
    monkeypatch,
) -> None:
    mod = _load_module()
    zip_path = tmp_path / "14100_yokohama-shi_city_2024_citygml_2_op.zip"
    _write_test_zip(zip_path, with_top_level_dir=False)
    derived_root = tmp_path / "derived-root"
    outline_root = tmp_path / "outline-root"
    called_with: list[list[str]] = []

    def fake_outline_main(argv: list[str]) -> int:
        called_with.append(list(argv))
        output_path = Path(argv[argv.index("--output-json") + 1])
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("[]", encoding="utf-8")
        return 0

    monkeypatch.setattr(mod.urban_debug_layer_from_citygml, "main", fake_outline_main)

    rc = mod.main(
        [
            str(zip_path),
            "--derived-root-dir",
            str(derived_root),
            "--urban-outline-root-dir",
            str(outline_root),
            "--write-urban-outline-json",
            "--workers",
            "1",
        ]
    )

    assert rc == 0
    derived_dir = derived_root / "14100_yokohama" / "bldg"
    assert (derived_dir / "tile-a.json").exists()
    assert called_with
    assert "--all-covered-towers" in called_with[0]
    assert (outline_root / "yokohama_urban_outline_layer.json").exists()
