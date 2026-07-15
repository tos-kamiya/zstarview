from __future__ import annotations

import json
from pathlib import Path
from zipfile import ZipFile

from zstarview.data.plateau_buildings import (
    build_download_url,
    catalog_archive_url,
    catalog_file_entries,
    catalog_year,
    format_binary_size,
    find_building_files,
    main,
    parse_citygml_buildings,
)


CITYGML = """<?xml version="1.0" encoding="UTF-8"?>
<core:CityModel xmlns:core="http://www.opengis.net/citygml/2.0"
    xmlns:gml="http://www.opengis.net/gml"
    xmlns:bldg="http://www.opengis.net/citygml/building/2.0">
  <gml:boundedBy>
    <gml:Envelope srsDimension="3">
      <gml:lowerCorner>35.4500 133.0500 0</gml:lowerCorner>
      <gml:upperCorner>35.4510 133.0510 100</gml:upperCorner>
    </gml:Envelope>
  </gml:boundedBy>
  <core:cityObjectMember>
    <bldg:Building gml:id="BLD_TEST_1">
      <bldg:measuredHeight>12.5</bldg:measuredHeight>
      <bldg:lod0RoofEdge>
        <gml:LinearRing>
          <gml:posList>35.4500 133.0500 0 35.4500 133.0501 0 35.4501 133.0501 0 35.4501 133.0500 0 35.4500 133.0500 0</gml:posList>
        </gml:LinearRing>
      </bldg:lod0RoofEdge>
    </bldg:Building>
  </core:cityObjectMember>
</core:CityModel>
"""


def _write_citygml_zip(path: Path) -> None:
    with ZipFile(path, "w") as archive:
        archive.writestr(
            "32201_2024_citygml_1/udx/bldg/53394500_bldg_6697.gml", CITYGML
        )


def test_build_download_url() -> None:
    assert build_download_url("32201", "latest").endswith(
        "/datacatalog/citygml/32201-latest/citygml.zip"
    )


def test_format_binary_size_uses_binary_units() -> None:
    assert format_binary_size(1024) == "1.00 KiB"
    assert format_binary_size(1024 * 1024) == "1.00 MiB"
    assert format_binary_size(3 * 1024**3) == "3.00 GiB"


def test_catalog_file_entries_reads_city_wrapped_bldg_files() -> None:
    payload = {
        "cities": [
            {
                "year": 2024,
                "files": {
                    "bldg": [{"code": "53394500", "url": "https://example.test/a.gml"}]
                },
            }
        ]
    }

    entries = catalog_file_entries(payload)

    assert len(entries) == 1
    assert entries[0]["code"] == "53394500"
    assert catalog_year(payload, "latest") == "2024"
    assert catalog_archive_url(payload) is None


def test_catalog_archive_url_reads_city_zip_url() -> None:
    payload = {"cities": [{"url": "https://example.test/citygml.zip"}]}

    assert catalog_archive_url(payload) == "https://example.test/citygml.zip"


def test_parse_citygml_buildings_reads_height_and_ring(tmp_path: Path) -> None:
    path = tmp_path / "building.gml"
    path.write_text(CITYGML, encoding="utf-8")

    buildings = parse_citygml_buildings(path)

    assert len(buildings) == 1
    assert buildings[0]["id"] == "BLD_TEST_1"
    assert buildings[0]["height_m"] == 12.5
    assert buildings[0]["height_source"] == "measuredHeight"
    assert len(buildings[0]["rings"]) == 1


def test_find_building_files_ignores_non_building_citygml(tmp_path: Path) -> None:
    building = tmp_path / "dataset" / "udx" / "bldg" / "building.gml"
    terrain = tmp_path / "dataset" / "udx" / "dem" / "terrain.gml"
    building.parent.mkdir(parents=True)
    terrain.parent.mkdir(parents=True)
    building.write_text(CITYGML, encoding="utf-8")
    terrain.write_text(CITYGML, encoding="utf-8")

    assert find_building_files(tmp_path) == (building,)


def test_main_converts_local_zip_to_overture_shape(tmp_path: Path) -> None:
    zip_path = tmp_path / "matsue.zip"
    output_root = tmp_path / "cache"
    _write_citygml_zip(zip_path)

    assert (
        main(
            [
                "--city-code",
                "32201",
                "--year",
                "2024",
                "--input-zip",
                str(zip_path),
                "--output-root",
                str(output_root),
            ]
        )
        == 0
    )

    derived_dir = output_root / "32201_2024" / "bldg"
    tile_path = derived_dir / "53394500_bldg_6697.json"
    payload = json.loads(tile_path.read_text(encoding="utf-8"))
    assert payload["source"]["format"] == "PLATEAU-CityGML"
    assert payload["buildings"][0]["height_m"] == 12.5
    assert (derived_dir / "tile_index.json").exists()
    metadata = json.loads((output_root / "32201_2024" / "cache_meta.json").read_text())
    assert metadata["status"] == "complete"
    assert metadata["building_count"] == 1
    assert metadata["metadata_schema_version"] == 1
    assert metadata["derived_tile_schema_version"] == 1
    assert metadata["converter"] == "zstarview-plateau-buildings"
