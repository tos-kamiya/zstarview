"""Tests for the PLATEAU cache layout and the subcommand-based preparation CLI."""

from __future__ import annotations

import errno
import json
from pathlib import Path
from urllib.error import HTTPError
from zipfile import ZipFile

import pytest

import zstarview.data.plateau_buildings as plateau_module
from zstarview.data.plateau_building_cache import (
    classify_plateau_caches,
    find_plateau_building_dirs,
    is_valid_plateau_cache,
)
from zstarview.data.plateau_buildings import hash_source_zip, main

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
      <bldg:measuredHeight>HEIGHT</bldg:measuredHeight>
      <bldg:lod0RoofEdge>
        <gml:LinearRing>
          <gml:posList>35.4500 133.0500 0 35.4500 133.0501 0 35.4501 133.0501 0 35.4501 133.0500 0 35.4500 133.0500 0</gml:posList>
        </gml:LinearRing>
      </bldg:lod0RoofEdge>
    </bldg:Building>
  </core:cityObjectMember>
</core:CityModel>
"""

CATALOG = {
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
MATCHING_METADATA = {
    "preparation_year": "2024",
    "registration_year": "2025",
    "source_spec": "3.0",
    "source_file_size_bytes": 10,
    "source_file_count": 1,
}


def _write_zip(path: Path, *, height: str = "12.5", with_building: bool = True) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(path, "w") as archive:
        if with_building:
            archive.writestr(
                "32201_2024_citygml_1/udx/bldg/53394500_bldg_6697.gml",
                CITYGML.replace("HEIGHT", height),
            )
        else:
            archive.writestr("readme.txt", "no buildings")
    return path


def _write_cache(
    dataset_dir: Path,
    *,
    city_code: str = "32201",
    lat: float = 35.45,
    lon: float = 133.05,
    **metadata: object,
) -> Path:
    derived_dir = dataset_dir / "bldg"
    derived_dir.mkdir(parents=True)
    bbox = {
        "min_lat": lat,
        "min_lon": lon,
        "max_lat": lat + 0.01,
        "max_lon": lon + 0.01,
    }
    (derived_dir / "tile.json").write_text(
        json.dumps({"schema_version": 3, "tile": {"id": "t", "bbox": bbox}}),
        encoding="utf-8",
    )
    (derived_dir / "tile_index.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "bbox": bbox,
                "tiles": [{"path": "tile.json", "bbox": bbox}],
            }
        ),
        encoding="utf-8",
    )
    payload: dict[str, object] = {
        "metadata_schema_version": 1,
        "derived_tile_schema_version": 3,
        "source": "PLATEAU-CityGML",
        "city_code": city_code,
        "status": "complete",
    }
    payload.update(metadata)
    (dataset_dir / "cache_meta.json").write_text(json.dumps(payload), encoding="utf-8")
    return dataset_dir


def _names(root: Path) -> set[str]:
    return {path.name for path in root.iterdir()}


@pytest.fixture
def fake_remote(monkeypatch, tmp_path: Path) -> list[Path]:
    """Serve CATALOG and a CityGML ZIP without network access."""

    source_zip = _write_zip(tmp_path / "remote" / "citygml.zip")
    downloads: list[Path] = []
    monkeypatch.setattr(plateau_module, "fetch_catalog", lambda *_a, **_k: CATALOG)
    monkeypatch.setattr(plateau_module, "_content_length", lambda _url: 1234)

    def fake_download(_url, destination: Path, **_kwargs):
        downloads.append(destination)
        destination.write_bytes(source_zip.read_bytes())
        return hash_source_zip(destination)

    monkeypatch.setattr(plateau_module, "download_file", fake_download)
    return downloads


def _no_prompt(monkeypatch) -> None:
    def fail(_prompt: str) -> str:
        raise AssertionError("unexpected confirmation prompt")

    monkeypatch.setattr("builtins.input", fail)


# --- cache layout ---------------------------------------------------------


def test_one_cache_is_read_per_municipality(tmp_path: Path) -> None:
    root = tmp_path / "cache"
    current = _write_cache(root / "32201", **MATCHING_METADATA)
    _write_cache(root / "32201_2025")
    _write_cache(root / "32201_2024.outdated-20260715142530")
    _write_cache(root / "32201_2024.replaced-20260715142530")
    _write_cache(root / "32201_latest")
    _write_cache(root / ".32201-abc")
    legacy_old = _write_cache(root / "13103_2024")
    legacy_new = _write_cache(root / "13103_2025")
    _write_cache(root / "13103_2026.incomplete-20260101000000")
    (root / "unrelated").mkdir()

    entries = {entry.city_code: entry for entry in classify_plateau_caches(root)}

    assert entries["32201"].selected == current
    assert not entries["32201"].legacy
    assert len(entries["32201"].removable) == 5
    assert entries["13103"].selected == legacy_new
    assert entries["13103"].legacy
    assert legacy_old in entries["13103"].removable
    assert find_plateau_building_dirs(
        observer_lat_deg=35.455,
        observer_lon_deg=133.055,
        radius_km=2.0,
        root_dir=root,
    ) == (legacy_new / "bldg", current / "bldg")


def test_invalid_current_cache_falls_back_to_legacy(tmp_path: Path) -> None:
    root = tmp_path / "cache"
    broken = _write_cache(root / "32201")
    (broken / "cache_meta.json").write_text("{}", encoding="utf-8")
    legacy = _write_cache(root / "32201_2024")

    (entry,) = classify_plateau_caches(root)

    assert entry.selected == legacy
    assert broken in entry.removable


def test_dot_prefixed_or_retired_caches_are_never_read(tmp_path: Path) -> None:
    root = tmp_path / "cache"
    _write_cache(root / ".32201-retired-abc")
    _write_cache(root / "32201_2024.outdated-20260715142530")

    (entry,) = classify_plateau_caches(root)

    assert entry.selected is None
    assert len(entry.removable) == 2


# --- prepare via catalog --------------------------------------------------


def test_prepare_downloads_into_municipality_directory(
    fake_remote: list[Path], tmp_path: Path
) -> None:
    root = tmp_path / "cache"
    temp_dir = tmp_path / "large-temp"

    assert (
        main(
            [
                "prepare",
                "32201",
                "--yes",
                "--output-root",
                str(root),
                "--temp-dir",
                str(temp_dir),
            ]
        )
        == 0
    )

    assert fake_remote[0].parent.parent == temp_dir
    assert _names(root) == {"32201"}
    assert is_valid_plateau_cache(root / "32201")
    metadata = json.loads((root / "32201" / "cache_meta.json").read_text())
    assert metadata["year"] == "2024"
    assert metadata["preparation_year"] == "2024"
    assert metadata["source_zip_sha256"] == hash_source_zip(
        tmp_path / "remote" / "citygml.zip"
    ).sha256
    assert metadata["source_zip_size_bytes"] > 0
    assert metadata["archive_url"] == "https://example.test/citygml.zip"
    assert not (root / "32201" / "source-citygml.zip").exists()
    assert not tuple(temp_dir.iterdir())


def test_prepare_replaces_legacy_layouts_and_leftovers_only_for_that_code(
    fake_remote: list[Path], tmp_path: Path, capsys
) -> None:
    root = tmp_path / "cache"
    _write_cache(root / "32201_2023")
    _write_cache(root / "32201_2023.outdated-20260715142530")
    _write_cache(root / "32201_latest")
    (root / ".32201-leftover").mkdir()
    other = _write_cache(root / "13103_2025.outdated-20260715142530", city_code="13103")

    assert main(["prepare", "32201", "--yes", "--output-root", str(root)]) == 0

    out = capsys.readouterr().out
    assert "remove after success:" in out
    assert f"Removed: {root / '32201_2023'}" in out
    assert _names(root) == {"32201", other.name}
    assert len(fake_remote) == 1


def test_prepare_skips_download_when_catalog_matches(
    fake_remote: list[Path], monkeypatch, tmp_path: Path, capsys
) -> None:
    root = tmp_path / "cache"
    _write_cache(root / "32201", **MATCHING_METADATA)
    _no_prompt(monkeypatch)

    assert main(["prepare", "32201", "--output-root", str(root)]) == 0

    assert fake_remote == []
    assert "cache is up to date" in capsys.readouterr().out


def test_catalog_check_ignores_retired_metadata(
    fake_remote: list[Path], tmp_path: Path
) -> None:
    # Regression: an .outdated-* copy with stale metadata used to be compared
    # against the catalog and forced a re-download of an up-to-date cache.
    root = tmp_path / "cache"
    _write_cache(root / "32201", **MATCHING_METADATA)
    _write_cache(
        root / "32201_2023.outdated-20260715142530",
        **{**MATCHING_METADATA, "preparation_year": "2023"},
        derived_tile_schema_version=1,
    )

    assert main(["prepare", "32201", "--yes", "--output-root", str(root)]) == 0

    assert fake_remote == []
    assert _names(root) == {"32201"}


def test_catalog_check_does_not_trust_retired_matching_metadata(
    fake_remote: list[Path], tmp_path: Path
) -> None:
    root = tmp_path / "cache"
    _write_cache(root / "32201_2024.outdated-20260715142530", **MATCHING_METADATA)

    assert main(["prepare", "32201", "--yes", "--output-root", str(root)]) == 0

    assert len(fake_remote) == 1
    assert _names(root) == {"32201"}


def test_prepare_migrates_matching_legacy_cache_without_download(
    fake_remote: list[Path], tmp_path: Path
) -> None:
    root = tmp_path / "cache"
    _write_cache(root / "32201_2024", **MATCHING_METADATA, marker="legacy")
    _write_cache(root / "32201_2023")

    assert main(["prepare", "32201", "--yes", "--output-root", str(root)]) == 0

    assert fake_remote == []
    assert _names(root) == {"32201"}
    metadata = json.loads((root / "32201" / "cache_meta.json").read_text())
    assert metadata["marker"] == "legacy"


def test_prepare_cancel_keeps_everything(
    fake_remote: list[Path], monkeypatch, tmp_path: Path, capsys
) -> None:
    root = tmp_path / "cache"
    _write_cache(root / "32201_2023")
    monkeypatch.setattr("builtins.input", lambda _prompt: "n")

    assert main(["prepare", "32201", "--output-root", str(root)]) == 0

    assert fake_remote == []
    assert _names(root) == {"32201_2023"}
    assert "Cancelled." in capsys.readouterr().out


def test_prepare_failure_keeps_existing_cache(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    root = tmp_path / "cache"
    existing = _write_cache(root / "32201", marker="old")
    _write_cache(root / "32201_2023")
    empty_zip = _write_zip(tmp_path / "empty.zip", with_building=False)

    assert (
        main(
            [
                "prepare",
                "32201",
                "--input-zip",
                str(empty_zip),
                "--yes",
                "--output-root",
                str(root),
            ]
        )
        == 1
    )

    assert "existing cache was kept" in capsys.readouterr().err
    assert _names(root) == {"32201", "32201_2023"}
    assert json.loads((existing / "cache_meta.json").read_text())["marker"] == "old"


def test_install_restores_old_cache_when_final_rename_fails(
    monkeypatch, tmp_path: Path
) -> None:
    root = tmp_path / "cache"
    _write_cache(root / "32201", marker="old")
    staging = _write_cache(root / ".32201-staging", marker="new")
    original_rename = Path.rename

    def flaky_rename(self: Path, target):  # type: ignore[no-untyped-def]
        if self == staging:
            raise OSError("simulated failure")
        return original_rename(self, target)

    monkeypatch.setattr(Path, "rename", flaky_rename)

    with pytest.raises(OSError, match="simulated"):
        plateau_module._install_cache(root, "32201", staging)

    metadata = json.loads((root / "32201" / "cache_meta.json").read_text())
    assert metadata["marker"] == "old"


# --- prepare via local ZIP ------------------------------------------------


def test_prepare_zip_confirms_then_reconverts_same_source_without_prompt(
    monkeypatch, tmp_path: Path
) -> None:
    root = tmp_path / "cache"
    zip_path = _write_zip(tmp_path / "matsue.zip")
    answers = iter(["y"])
    monkeypatch.setattr("builtins.input", lambda _prompt: next(answers))

    assert (
        main(
            [
                "prepare",
                "32201",
                "--input-zip",
                str(zip_path),
                "--keep-zip",
                "--output-root",
                str(root),
            ]
        )
        == 0
    )
    kept_zip = root / "32201" / "source-citygml.zip"
    assert kept_zip.is_file()
    first = json.loads((root / "32201" / "cache_meta.json").read_text())
    assert first["source_zip_sha256"] == hash_source_zip(zip_path).sha256
    assert first["year"] is None

    # Reconvert from the ZIP kept inside the cache: same checksum, no prompt.
    _no_prompt(monkeypatch)
    assert (
        main(
            [
                "prepare",
                "32201",
                "--input-zip",
                str(kept_zip),
                "--keep-zip",
                "--output-root",
                str(root),
            ]
        )
        == 0
    )
    assert _names(root) == {"32201"}
    assert (root / "32201" / "source-citygml.zip").is_file()


def test_prepare_zip_keeps_catalog_fields_for_same_source(tmp_path: Path) -> None:
    root = tmp_path / "cache"
    zip_path = _write_zip(tmp_path / "matsue.zip")
    _write_cache(
        root / "32201",
        **MATCHING_METADATA,
        source_zip_sha256=hash_source_zip(zip_path).sha256,
        archive_url="https://example.test/citygml.zip",
    )

    assert (
        main(["prepare", "32201", "--input-zip", str(zip_path), "--output-root", str(root)])
        == 0
    )

    metadata = json.loads((root / "32201" / "cache_meta.json").read_text())
    assert metadata["preparation_year"] == "2024"
    assert metadata["archive_url"] == "https://example.test/citygml.zip"
    assert metadata["building_count"] == 1


def test_prepare_zip_with_different_source_asks_before_replacing(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    root = tmp_path / "cache"
    _write_cache(root / "32201", source_zip_sha256="0" * 64)
    zip_path = _write_zip(tmp_path / "other.zip", height="30")
    prompts: list[str] = []

    def answer(prompt: str) -> str:
        prompts.append(prompt)
        return "no"

    monkeypatch.setattr("builtins.input", answer)

    assert (
        main(["prepare", "32201", "--input-zip", str(zip_path), "--output-root", str(root)])
        == 0
    )

    assert prompts
    assert "replace current cache" in capsys.readouterr().out
    assert json.loads((root / "32201" / "cache_meta.json").read_text())[
        "source_zip_sha256"
    ] == "0" * 64


def test_prepare_zip_rejects_multiple_codes(tmp_path: Path, capsys) -> None:
    zip_path = _write_zip(tmp_path / "matsue.zip")

    assert main(["prepare", "32201,32202", "--input-zip", str(zip_path)]) == 2
    assert "one city code" in capsys.readouterr().err


# --- errors -------------------------------------------------------------------


def test_old_flag_interface_reports_migration(capsys) -> None:
    assert main(["--city-code", "32201", "--yes"]) == 2
    assert "subcommands" in capsys.readouterr().err


def test_missing_single_catalog_reports_without_traceback(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    def not_found(*_args, **_kwargs):
        raise HTTPError("https://example.test", 404, "not found", {}, None)

    monkeypatch.setattr(plateau_module, "fetch_catalog", not_found)

    assert main(["prepare", "32203", "--yes", "--output-root", str(tmp_path)]) == 1
    captured = capsys.readouterr()
    assert "HTTP 404" in captured.err
    assert "Traceback" not in captured.err


def test_missing_catalog_in_multi_code_request_is_skipped(
    fake_remote: list[Path], monkeypatch, tmp_path: Path
) -> None:
    def catalog(city_code: str, **_kwargs):
        if city_code == "13100":
            raise HTTPError("https://example.test", 404, "not found", {}, None)
        return CATALOG

    monkeypatch.setattr(plateau_module, "fetch_catalog", catalog)
    root = tmp_path / "cache"

    assert main(["prepare", "13100,32201", "--yes", "--output-root", str(root)]) == 0
    assert _names(root) == {"32201"}


def test_prepare_suggests_temp_dir_when_temporary_filesystem_is_full(
    fake_remote: list[Path], monkeypatch, tmp_path: Path, capsys
) -> None:
    def full(*_args, **_kwargs):
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.setattr(plateau_module, "download_file", full)

    assert main(["prepare", "32201", "--yes", "--output-root", str(tmp_path)]) == 1
    err = capsys.readouterr().err
    assert "existing cache was kept" in err
    assert (
        '  zstarview-download-plateau-buildings prepare 32201 '
        '--temp-dir "$HOME/zstarview-tmp"'
    ) in err


# --- list and remove ------------------------------------------------------


def test_list_shows_current_legacy_and_removable(tmp_path: Path, capsys) -> None:
    root = tmp_path / "cache"
    _write_cache(root / "32201", preparation_year="2024")
    _write_cache(root / "32201_2023")
    _write_cache(root / "13103_2025", city_code="13103", preparation_year="2025")

    assert main(["list", "--output-root", str(root)]) == 0

    assert capsys.readouterr().out == (
        f"13103 2025 legacy {root / '13103_2025'}\n"
        f"32201 2024 current {root / '32201'}\n"
        f"32201 - removable {root / '32201_2023'}\n"
    )

    assert main(["list", "--jsonl", "--output-root", str(root)]) == 0
    rows = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert rows[1]["status"] == "current"
    assert rows[1]["preparation_year"] == "2024"
    assert rows[2] == {
        "city_code": "32201",
        "path": str(root / "32201_2023"),
        "status": "removable",
    }


def test_remove_unused_keeps_selected_caches(tmp_path: Path) -> None:
    root = tmp_path / "cache"
    _write_cache(root / "32201")
    _write_cache(root / "32201_2023")
    (root / ".13103-staging").mkdir()
    _write_cache(root / "13103_2025", city_code="13103")

    assert main(["remove", "--unused", "--yes", "--output-root", str(root)]) == 0

    assert _names(root) == {"32201", "13103_2025"}


def test_remove_city_code_removes_all_of_its_directories(
    monkeypatch, tmp_path: Path
) -> None:
    root = tmp_path / "cache"
    _write_cache(root / "32201")
    _write_cache(root / "32201_2023")
    _write_cache(root / "13103_2025", city_code="13103")
    monkeypatch.setattr("builtins.input", lambda _prompt: "y")

    assert main(["remove", "32201", "--output-root", str(root)]) == 0

    assert _names(root) == {"13103_2025"}


def test_remove_requires_target(capsys) -> None:
    assert main(["remove"]) == 2
    assert "--unused" in capsys.readouterr().err
