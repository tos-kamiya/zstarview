from __future__ import annotations

import json
from pathlib import Path
from urllib.error import HTTPError
from zipfile import ZipFile

import zstarview.data.plateau_buildings as plateau_module
from zstarview.data.plateau_buildings import (
    build_download_url,
    catalog_archive_url,
    catalog_file_entries,
    catalog_file_size_bytes,
    catalog_registration_year,
    catalog_year,
    format_binary_size,
    find_building_files,
    main,
    parse_city_codes,
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

CITYGML_WITH_ROOF_SURFACE = CITYGML.replace(
    "    </bldg:Building>",
    """      <bldg:boundedBy>
        <bldg:RoofSurface>
          <bldg:lod2MultiSurface>
            <gml:MultiSurface>
              <gml:surfaceMember>
                <gml:Polygon>
                  <gml:exterior>
                    <gml:LinearRing>
                      <gml:posList>35.4500 133.0500 10 35.4500 133.0501 10 35.4501 133.0501 12 35.4501 133.0500 12 35.4500 133.0500 10</gml:posList>
                    </gml:LinearRing>
                  </gml:exterior>
                </gml:Polygon>
              </gml:surfaceMember>
            </gml:MultiSurface>
          </bldg:lod2MultiSurface>
        </bldg:RoofSurface>
      </bldg:boundedBy>
    </bldg:Building>""",
)


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


def test_parse_city_codes_expands_ranges_and_comma_separated_values() -> None:
    assert parse_city_codes("13100-13102, 32201,13102") == (
        "13100",
        "13101",
        "13102",
        "32201",
    )


def test_parse_city_codes_rejects_invalid_ranges() -> None:
    import pytest

    with pytest.raises(ValueError, match="reversed"):
        parse_city_codes("13102-13100")
    with pytest.raises(ValueError, match="five-digit"):
        parse_city_codes("13100,1234")


def test_main_skips_missing_city_code_in_multi_code_request(
    monkeypatch, tmp_path: Path
) -> None:
    prepared: list[str] = []

    def fake_prepare(_args, city_code: str) -> int:
        if city_code == "13100":
            raise HTTPError("https://example.test", 404, "not found", {}, None)
        prepared.append(city_code)
        return 0

    monkeypatch.setattr(plateau_module, "_prepare_city_code", fake_prepare)
    monkeypatch.setattr(
        plateau_module,
        "_preflight_city_codes",
        lambda _args, city_codes: city_codes,
    )

    assert (
        plateau_module.main(
            [
                "--city-code",
                "13100,13101",
                "--output-root",
                str(tmp_path),
                "--yes",
            ]
        )
        == 0
    )
    assert prepared == ["13101"]


def test_main_reports_missing_city_code_without_traceback(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    monkeypatch.setattr(
        plateau_module,
        "_prepare_city_code",
        lambda _args, _city_code: (_ for _ in ()).throw(
            HTTPError("https://example.test", 404, "not found", {}, None)
        ),
    )

    assert (
        plateau_module.main(
            ["--city-code", "32203", "--output-root", str(tmp_path), "--yes"]
        )
        == 1
    )
    assert capsys.readouterr().out == (
        "PLATEAU catalog request failed for city code 32203: HTTP 404. "
        "PLATEAU building data may not be available for this municipality.\n"
    )


def test_main_reports_no_available_city_codes_without_traceback(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    monkeypatch.setattr(
        plateau_module,
        "_preflight_city_codes",
        lambda _args, _city_codes: (_ for _ in ()).throw(
            ValueError(
                "No PLATEAU building catalogs found for the requested city codes"
            )
        ),
    )

    assert (
        plateau_module.main(
            [
                "--city-code",
                "27100,27101",
                "--output-root",
                str(tmp_path),
                "--yes",
            ]
        )
        == 1
    )
    captured = capsys.readouterr()
    assert captured.err == (
        "Error: No PLATEAU building catalogs found for the requested city codes\n"
    )
    assert "Traceback" not in captured.err


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
    assert catalog_registration_year(payload, "unknown") == "unknown"
    assert catalog_file_size_bytes(entries) is None
    assert catalog_archive_url(payload) is None


def test_cache_metadata_matches_catalog_fields() -> None:
    payload = {
        "cities": [
            {
                "year": 2025,
                "registrationYear": 2026,
                "spec": "3.0",
                "files": {
                    "bldg": [
                        {
                            "code": "53394500",
                            "url": "https://example.test/a.gml",
                            "fileSize": 120,
                        },
                        {
                            "code": "53394501",
                            "url": "https://example.test/b.gml",
                            "fileSize": 80,
                        },
                    ]
                },
            }
        ]
    }
    entries = catalog_file_entries(payload)
    metadata = {
        "metadata_schema_version": 1,
        "derived_tile_schema_version": 3,
        "preparation_year": "2025",
        "registration_year": "2026",
        "source_spec": "3.0",
        "source_file_size_bytes": 200,
        "source_file_count": 2,
    }

    assert plateau_module._cache_matches_catalog(metadata, payload, entries, "latest")
    metadata["source_file_size_bytes"] = 201
    assert not plateau_module._cache_matches_catalog(
        metadata, payload, entries, "latest"
    )


def test_parse_citygml_buildings_extracts_lod2_roof_surfaces(tmp_path: Path) -> None:
    path = tmp_path / "roof.gml"
    path.write_text(CITYGML_WITH_ROOF_SURFACE, encoding="utf-8")

    buildings = parse_citygml_buildings(path)

    assert len(buildings) == 1
    assert buildings[0]["geometry_lod"] == 2
    assert buildings[0]["roof_surfaces"] == [
        [
            [133.05, 35.45, 10.0],
            [133.0501, 35.45, 10.0],
            [133.0501, 35.4501, 12.0],
            [133.05, 35.4501, 12.0],
            [133.05, 35.45, 10.0],
        ]
    ]


def test_main_replaces_outdated_remote_cache(monkeypatch, tmp_path: Path) -> None:
    zip_path = tmp_path / "matsue.zip"
    _write_citygml_zip(zip_path)
    output_root = tmp_path / "cache"
    dataset_dir = output_root / "32201_2024"
    dataset_dir.mkdir(parents=True)
    (dataset_dir / "cache_meta.json").write_text(
        json.dumps({"city_code": "32201", "status": "complete"}),
        encoding="utf-8",
    )
    catalog = {
        "cities": [
            {
                "year": 2024,
                "registrationYear": 2025,
                "spec": "3.0",
                "url": "https://example.test/citygml.zip",
                "files": {
                    "bldg": [
                        {
                            "code": "53394500",
                            "url": "https://example.test/a.gml",
                            "fileSize": 10,
                        }
                    ]
                },
            }
        ]
    }
    monkeypatch.setattr(
        plateau_module, "fetch_catalog", lambda *_args, **_kwargs: catalog
    )
    monkeypatch.setattr(plateau_module, "_content_length", lambda _url: None)
    monkeypatch.setattr(
        plateau_module,
        "download_file",
        lambda _url, destination, **_kwargs: (
            zip_path.read_bytes()
            and destination.write_bytes(zip_path.read_bytes())
            or zip_path.stat().st_size
        ),
    )

    assert (
        main(
            [
                "--city-code",
                "32201",
                "--yes",
                "--output-root",
                str(output_root),
            ]
        )
        == 0
    )
    assert (output_root / "32201_2024" / "cache_meta.json").exists()
    assert tuple(output_root.glob("32201_2024.outdated-*"))


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


def test_parse_citygml_buildings_prefers_lod1_over_lod0(tmp_path: Path) -> None:
    path = tmp_path / "lod1-building.gml"
    path.write_text(
        CITYGML.replace("lod0RoofEdge", "lod1Solid"),
        encoding="utf-8",
    )

    buildings = parse_citygml_buildings(path)

    assert buildings[0]["geometry_lod"] == 1


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
    assert metadata["derived_tile_schema_version"] == 3
    assert metadata["converter"] == "zstarview-plateau-buildings"
    assert metadata["geometry_mode"] == "lod0-footprint"
    assert metadata["max_geometry_lod"] == 0
    assert metadata["lod0_building_count"] == 1
