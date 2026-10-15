from __future__ import annotations

import hashlib
import io
import json
import socket
from dataclasses import asdict
from datetime import timedelta
from pathlib import Path

import pytest
from PIL import Image

from zstarview.aircraft.opensky import build_observer_bbox
from zstarview.cli.export_replay import query_cache_ranges
from zstarview.cli.export_replay_args import parse_replay_args
from zstarview.clouddisc.types import CloudSourceData, SourceKey
from zstarview.replay.inventory import Inventory
from zstarview.replay.models import MINUTE, HistoryRecord, parse_utc
from zstarview.replay.report import build_report, format_report
from zstarview.replay.sequence import _chunk_frame_indices
from zstarview.replay.selection import HistoryIndex, coverage

T = parse_utc("2026-10-14T00:00:00Z")


def record(name="cloud-native", start=T, source=T, end=T + 30 * MINUTE, **kwargs):
    return HistoryRecord(name, Path(str(start)), source, start, end, **kwargs)


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))


def aircraft(root, fetched=T, contacts=(), suffix="one"):
    payload = {
        "bbox": asdict(build_observer_bbox(35, 133)),
        "fetched_at_utc": fetched.isoformat(),
        "snapshots": [
            dict(
                icao24=str(i),
                latitude=35,
                longitude=133,
                last_contact_unix=int(c.timestamp()),
            )
            for i, c in enumerate(contacts)
        ],
    }
    path = (
        root
        / "aircraft/opensky/history"
        / fetched.strftime("%Y/%m/%d/%H/%M")
        / f"bbox-{suffix}.json"
    )
    write_json(path, payload)
    return path


def test_future_fetch_cannot_be_used_before_acquisition():
    index = HistoryIndex([record(start=T + 2 * MINUTE)])
    assert index.select(T) is None
    assert coverage(index, T, T + 4 * MINUTE) == [(T + 2 * MINUTE, T + 4 * MINUTE)]


def test_inclusive_freshness_and_exclusive_end():
    index = HistoryIndex([record()])
    assert index.select(T + 30 * MINUTE) is not None
    assert index.select(T + 30 * MINUTE + timedelta(microseconds=1)) is None
    assert coverage(index, T, T + 2 * MINUTE) == [(T, T + 2 * MINUTE)]


def test_gaps_are_preserved_and_adjacent_ranges_merge():
    index = HistoryIndex(
        [
            record("satellite-iss", end=T),
            record(
                "satellite-iss",
                start=T + 2 * MINUTE,
                source=T + 2 * MINUTE,
                end=T + 3 * MINUTE,
            ),
        ]
    )
    assert coverage(index, T, T + 5 * MINUTE) == [
        (T, T + 4 * MINUTE),
    ]


def test_latest_aircraft_must_not_fall_back_to_old_snapshot():
    old = record("aircraft", end=T + 10 * MINUTE, empty=True)
    new = record("aircraft", start=T + MINUTE, end=T, contact_times=(T - 10 * MINUTE,))
    assert HistoryIndex([old, new]).select(T + 2 * MINUTE) is None


def test_aircraft_empty_snapshot_is_valid_and_deduplicated(tmp_path):
    path = aircraft(tmp_path)
    write_json(tmp_path / "aircraft/opensky/latest.json", json.loads(path.read_text()))
    inventory = Inventory(tmp_path, 35, 133, 1.7)
    inventory.scan()
    assert len(inventory.records["aircraft"]) == 1
    index = HistoryIndex(inventory.records["aircraft"])
    assert index.select(T).empty
    assert index.select(T + 11 * MINUTE) is None


def test_aircraft_contacts_and_scope(tmp_path):
    aircraft(tmp_path, fetched=T + MINUTE, contacts=(T,))
    inventory = Inventory(tmp_path, 35, 133, 1.7)
    inventory.scan()
    index = HistoryIndex(inventory.records["aircraft"])
    assert index.select(T) is None
    assert index.select(T + MINUTE) is not None
    assert index.select(T + 11 * MINUTE) is None
    wrong = Inventory(tmp_path, 36, 133, 1.7)
    wrong.scan()
    assert not wrong.records["aircraft"]
    assert wrong.diagnostics["aircraft"]["scope-mismatch"] == 1


def geo(root, fetched=T + 2 * MINUTE):
    path = root / "geosatellite/raw/history/2026/10/14/00/00/infrared.json"
    buffer = io.BytesIO()
    Image.new("RGB", (4, 4)).save(buffer, format="PNG")
    data = buffer.getvalue()
    image = path.with_suffix(".png")
    image.parent.mkdir(parents=True, exist_ok=True)
    image.write_bytes(data)
    write_json(
        path,
        dict(
            schema_version=1,
            captured_at_utc=T.isoformat(),
            fetched_at_utc=fetched.isoformat(),
            kind="infrared",
            area="europe",
            size="normal",
            sha256=hashlib.sha256(data).hexdigest(),
            image_path=str(image.relative_to(root / "geosatellite")),
        ),
    )
    return image


def test_geo_metadata_does_not_verify_payload_digest(tmp_path):
    image = geo(tmp_path)
    inventory = Inventory(tmp_path, 51, 0, 1.7)
    inventory.scan()
    index = HistoryIndex(inventory.records["cloud-geo"])
    assert index.select(T) is None
    assert index.select(T + 2 * MINUTE) is not None
    assert index.select(T + 60 * MINUTE) is not None
    image.write_bytes(b"broken")
    inventory = Inventory(tmp_path, 51, 0, 1.7)
    inventory.scan()
    assert len(inventory.records["cloud-geo"]) == 1
    image.unlink()
    inventory = Inventory(tmp_path, 51, 0, 1.7)
    inventory.scan()
    assert not inventory.records["cloud-geo"]
    assert inventory.diagnostics["cloud-geo"]["corrupt-or-missing"] == 1


def test_path_escape_is_rejected(tmp_path):
    image = geo(tmp_path)
    path = image.with_suffix(".json")
    payload = json.loads(path.read_text())
    payload["image_path"] = "../../outside.png"
    write_json(path, payload)
    inventory = Inventory(tmp_path, 51, 0, 1.7)
    inventory.scan()
    assert not inventory.records["cloud-geo"]


def test_common_intersection_and_disabled_layers(tmp_path):
    inventory = Inventory(tmp_path, 35, 133, 1.7)
    inventory.records["cloud-native"] = [record(end=T + MINUTE)]
    inventory.records["aircraft"] = [record("aircraft", start=T + MINUTE, empty=True)]
    result = build_report(
        inventory, {"cloud-native", "aircraft"}, start=T, end=T + 3 * MINUTE
    )
    assert result["common_ranges"] == [
        {"start": "2026-10-14T00:01:00Z", "end": "2026-10-14T00:02:00Z"}
    ]
    assert build_report(inventory, set())["common_status"] == "not-applicable"
    assert (
        build_report(inventory, {"satellite-iss"})["common_status"] == "missing-history"
    )


@pytest.mark.parametrize(
    "options",
    [
        [],
        ["--start", "2026-01-01T00:00:00"],
        ["--end", "2026-01-01T00:00:01Z"],
        ["-a", "nan"],
        ["-c", "2"],
        ["--tropical-cyclone-opacity", "0.2"],
        ["--jobs", "0"],
        ["--jobs", "9"],
        ["--start", "2026-01-02T00:00Z", "--end", "2026-01-01T00:00Z"],
    ],
)
def test_cli_errors(options):
    with pytest.raises(SystemExit) as error:
        parse_replay_args(options if not options else ["--list-cache-ranges", *options])
    assert error.value.code == 1


def test_cli_shared_options():
    args = parse_replay_args(
        [
            "matsue",
            "--list-cache-ranges",
            "-A",
            "30",
            "-Z",
            "S",
            "-a",
            "0.5",
            "-B",
            "--height-add-m",
            "12",
            "--image-size",
            "1920,1080",
            "--jobs",
            "3",
        ]
    )
    assert (args.view_center_alt, args.view_center_az, args.aircraft_opacity) == (
        30,
        180,
        0.5,
    )
    assert args.use_building_top and args.observer_height_m == 12
    assert args.image_size == (1920, 1080)
    assert args.jobs == 3


def test_cli_replay_parallel_jobs_default_to_two():
    args = parse_replay_args(["matsue", "--list-cache-ranges"])
    assert args.jobs == 2


def test_render_chunks_split_at_cloud_record_changes_and_ten_frames():
    cloud_a = CloudSourceData(SourceKey("sat", T), None, "sat", "product", T, [])
    cloud_b = CloudSourceData(
        SourceKey("sat", T + 12 * MINUTE), None, "sat", "product", T, []
    )
    records = [
        HistoryRecord(
            "cloud-native",
            Path("record-a" if index < 12 else "record-b"),
            T,
            T,
            T + 24 * MINUTE,
        )
        for index in range(22)
    ]
    plans = [{"cloud": record} for record in records]
    decoded = {
        ("cloud-native", "record-a"): cloud_a,
        ("cloud-native", "record-b"): cloud_b,
    }

    chunks = _chunk_frame_indices(list(range(22)), plans, decoded, voxel_enabled=True)

    assert chunks == [list(range(10)), [10, 11], list(range(12, 22))]


def test_query_empty_cache_has_no_side_effects(tmp_path, monkeypatch):
    import zstarview.location_resolver.resolve as resolver

    def forbidden(*args, **kwargs):
        raise AssertionError("Unexpected I/O")

    monkeypatch.setattr(socket, "create_connection", forbidden)
    monkeypatch.setattr(resolver, "fetch_copernicus_dem", forbidden)
    monkeypatch.setattr(resolver, "save_last_city", forbidden)
    root = tmp_path / "absent"
    args = parse_replay_args(
        ["35,133", "--list-cache-ranges", "--cache-dir", str(root)]
    )
    result = query_cache_ranges(args)
    assert not root.exists()
    assert result["common_status"] == "missing-history"
    assert not result["complete_output_ready"]
    json.loads(json.dumps(result))
    assert "Coordinates: 35.0, 133.0" in format_report(result)


@pytest.mark.parametrize("offline", [False, True])
def test_place_resolution_matches_gui_online_first_and_cache_fallback(
    tmp_path, monkeypatch, offline
):
    import urllib.error

    from zstarview.location_resolver import place_search
    from zstarview.location_resolver.place_cache import (
        build_place_cache_key,
        load_place_cache,
        save_place_cache,
    )

    key = build_place_cache_key("station", None, "en")
    save_place_cache(
        key,
        (
            {
                "lat": 35,
                "lon": 133,
                "name": "Old Station",
                "display_name": "Old Station",
                "importance": 1,
            },
        ),
        original_query="station",
        fetched_at_utc=T,
        cache_root=tmp_path / "nominatim/schema-1",
    )
    write_json(
        tmp_path / "replay/v1/launches/launch.json",
        {
            "schema_version": 1,
            "started_at_utc": T.isoformat(),
            "location": {
                "input": {"query": "station", "language": "en", "countrycode": None},
                "resolved": {"latitude_deg": 30, "longitude_deg": 130},
            },
        },
    )
    calls = []

    def fetch(*args, **kwargs):
        calls.append(1)
        if offline:
            raise urllib.error.URLError("offline")
        return [
            {
                "lat": "34",
                "lon": "132",
                "display_name": "Lower Station",
                "importance": 0.5,
            },
            {
                "lat": "36",
                "lon": "134",
                "display_name": "Fresh Station",
                "importance": 1,
            },
        ]

    monkeypatch.setattr(place_search.nominatim, "_fetch", fetch)
    args = parse_replay_args(
        ["-p", "station", "--list-cache-ranges", "--cache-dir", str(tmp_path)]
    )
    result = query_cache_ranges(args)
    assert calls == [1]
    assert result["location"]["lat"] == (35 if offline else 36)
    assert result["location"]["lon"] == (133 if offline else 134)
    cached = load_place_cache(key, cache_root=tmp_path / "nominatim/schema-1")
    assert cached is not None
    assert cached.results[0]["lat"] == (35 if offline else 36)
    assert cached.fetched_at_utc == T if offline else cached.fetched_at_utc > T


def test_native_cloud_artifact_and_incompatible_schema(tmp_path):
    import numpy as np

    from zstarview.clouddisc.altaz_grid import CloudAltAzGrid
    from zstarview.clouddisc.history_cache import save_cloud_source_history
    from zstarview.clouddisc.types import CloudSourceData, SourceKey

    key = SourceKey("HIMAWARI", T)
    grid = CloudAltAzGrid(
        np.zeros((2, 2)),
        np.zeros((2, 2), dtype=np.uint8),
        0,
        90,
        0,
        360,
        35,
        133,
        "HIMAWARI",
        "test",
        T,
        (1,),
        key,
        1.0,
    )
    source = CloudSourceData(
        key, np.zeros((2, 2)), "HIMAWARI", "test", T, [], altaz_grid=grid
    )
    path = save_cloud_source_history(
        source, fetched_at_utc=T + MINUTE, cache_root=tmp_path
    )
    inventory = Inventory(tmp_path, 35, 133, 1.7)
    inventory.scan()
    index = HistoryIndex(inventory.records["cloud-native"])
    assert index.select(T) is None
    assert index.select(T + MINUTE) is not None
    payload = json.loads(path.read_text())
    payload["schema_version"] = 999
    write_json(path, payload)
    inventory = Inventory(tmp_path, 35, 133, 1.7)
    inventory.scan()
    assert inventory.diagnostics["cloud-native"]["incompatible"] == 1


def test_horizons_vector_scope_and_partial_group_diagnostic(tmp_path):
    from zstarview.satellite_constants import SATELLITE_CACHE_FORMAT_VERSION
    from zstarview.satellites.cache import satellite_cache_scope_key
    from zstarview.satellites.fetch import HORIZONS_TARGETS_BY_KEY

    scope = satellite_cache_scope_key(
        observer_lat=35, observer_lon=133, observer_height_m=1.7
    )
    records = [
        dict(
            OBJECT_NAME=target.label,
            EPOCH=T.isoformat(),
            HORIZONS_X_KM=1,
            HORIZONS_Y_KM=2,
            HORIZONS_Z_KM=3,
            HORIZONS_VX_KM_S=1,
            HORIZONS_VY_KM_S=2,
            HORIZONS_VZ_KM_S=3,
        )
        for target in HORIZONS_TARGETS_BY_KEY["horizons"]
    ]
    path = tmp_path / "satellites/celestrak" / scope / "horizons.json"
    write_json(
        path,
        dict(
            cache_format_version=SATELLITE_CACHE_FORMAT_VERSION,
            group_key="horizons",
            fetched_at_utc=(T + MINUTE).isoformat(),
            element_epoch_utc=T.isoformat(),
            records=records,
        ),
    )
    inventory = Inventory(tmp_path, 35, 133, 1.7)
    inventory.scan()
    index = HistoryIndex(inventory.records["satellite-horizons"])
    assert index.select(T) is None
    assert index.select(T + MINUTE) is not None
    wrong = Inventory(tmp_path, 35, 133, 2.0)
    wrong.scan()
    assert wrong.diagnostics["satellite-horizons"]["scope-mismatch"] == 1
    payload = json.loads(path.read_text())
    payload["records"] = payload["records"][:-1]
    write_json(path, payload)
    inventory = Inventory(tmp_path, 35, 133, 1.7)
    inventory.scan()
    assert inventory.diagnostics["satellite-horizons"]["partial-group"] == 1


def test_json_stdout_and_error_exit(tmp_path, monkeypatch, capsys):
    from zstarview.cli.export_replay import main

    monkeypatch.setattr(
        "sys.argv",
        [
            "zstarview-export-replay",
            "35,133",
            "--list-cache-ranges",
            "--json",
            "--cache-dir",
            str(tmp_path),
        ],
    )
    main()
    assert json.loads(capsys.readouterr().out)["schema_version"] == 1
    monkeypatch.setattr(
        "sys.argv", ["zstarview-export-replay", "999,133", "--list-cache-ranges"]
    )
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 1


def test_inventory_permission_error_is_not_empty(tmp_path, monkeypatch):
    import zstarview.replay.inventory as module

    def forbidden(*args, **kwargs):
        raise PermissionError("denied")

    monkeypatch.setattr(module.os, "walk", forbidden)
    with pytest.raises(PermissionError):
        module.json_paths(tmp_path)


def test_readonly_policy_is_reset_after_failure(monkeypatch, tmp_path):
    import zstarview.location_resolver.resolve as module

    def fail(*args, **kwargs):
        raise ValueError("test failure")

    policy = module.LocationDataPolicy(fail, fail, tmp_path)
    with pytest.raises(ValueError, match="test failure"):
        module.resolve_launch_location("35,133", data_policy=policy)
    assert module._location_data_policy.get() is None


def test_query_in_fresh_process_does_not_create_cache(tmp_path):
    import os
    import subprocess
    import sys

    cache_home = tmp_path / "cache"
    cache_home.mkdir()
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "zstarview.cli.export_replay",
            "35,133",
            "--list-cache-ranges",
            "--json",
        ],
        env=dict(os.environ, XDG_CACHE_HOME=str(cache_home)),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)["schema_version"] == 1
    assert not list(cache_home.rglob("*"))


def test_history_directory_timestamp_mismatch_is_rejected(tmp_path):
    path = aircraft(tmp_path)
    payload = json.loads(path.read_text())
    payload["fetched_at_utc"] = (T + MINUTE).isoformat()
    write_json(path, payload)
    inventory = Inventory(tmp_path, 35, 133, 1.7)
    inventory.scan()
    assert not inventory.records["aircraft"]
    assert inventory.diagnostics["aircraft"]["corrupt-or-missing"] == 1


def test_iss_omm_utc_epoch_without_suffix(tmp_path):
    from zstarview.satellite_constants import SATELLITE_CACHE_FORMAT_VERSION

    record = {
        "OBJECT_NAME": "ISS (ZARYA)",
        "OBJECT_ID": "1998-067A",
        "EPOCH": "2026-10-14T00:00:00.000000",
        "MEAN_MOTION": "15.5",
        "ECCENTRICITY": "0.0006703",
        "INCLINATION": "51.64",
        "RA_OF_ASC_NODE": "257.7",
        "ARG_OF_PERICENTER": "130.5",
        "MEAN_ANOMALY": "325.0",
        "EPHEMERIS_TYPE": "0",
        "CLASSIFICATION_TYPE": "U",
        "NORAD_CAT_ID": "25544",
        "ELEMENT_SET_NO": "999",
        "REV_AT_EPOCH": "12345",
        "BSTAR": "0.00010234",
        "MEAN_MOTION_DOT": "0.00002182",
        "MEAN_MOTION_DDOT": "0.00000000",
    }
    write_json(
        tmp_path / "satellites/celestrak/iss.json",
        {
            "cache_format_version": SATELLITE_CACHE_FORMAT_VERSION,
            "group_key": "iss",
            "fetched_at_utc": (T + MINUTE).isoformat(),
            "element_epoch_utc": T.isoformat(),
            "records": [record],
        },
    )
    inventory = Inventory(tmp_path, 35, 133, 1.7)
    inventory.scan()
    index = HistoryIndex(inventory.records["satellite-iss"])
    assert index.select(T) is None
    assert index.select(T + MINUTE) is not None
    assert index.select(T + timedelta(days=1, minutes=2)) is None


def test_query_reads_only_json_metadata_even_with_all_layers_enabled(
    tmp_path, monkeypatch
):
    import builtins
    import pickle

    import rasterio

    from zstarview.location_resolver import resolve as resolver

    artifact = tmp_path / "source_history/2026/10/14/00/00/source.pkl"
    artifact.parent.mkdir(parents=True)
    artifact.write_bytes(b"payload contents must never be inspected")
    write_json(
        artifact.with_suffix(".json"),
        {
            "schema_version": 1,
            "source_time_utc": T.isoformat(),
            "fetched_at_utc": T.isoformat(),
            "observer_latitude_deg": 35,
            "observer_longitude_deg": 133,
            "artifact_path": str(artifact.relative_to(tmp_path)),
            "artifact_sha256": "a" * 64,
        },
    )
    (tmp_path / "de442s.bsp").write_bytes(
        b"not a valid SPK; metadata listing must not open it"
    )
    original_path_open = Path.open
    original_open = builtins.open
    forbidden_suffixes = {".pkl", ".png", ".tif", ".npz", ".bsp"}

    def path_open(path, *args, **kwargs):
        assert path.suffix not in forbidden_suffixes, f"Unexpected payload read: {path}"
        return original_path_open(path, *args, **kwargs)

    def file_open(path, *args, **kwargs):
        if isinstance(path, (str, Path)):
            assert Path(path).suffix not in forbidden_suffixes
        return original_open(path, *args, **kwargs)

    def forbidden(*args, **kwargs):
        raise AssertionError("Payload validation is forbidden in metadata-only queries")

    monkeypatch.setattr(Path, "open", path_open)
    monkeypatch.setattr(builtins, "open", file_open)
    monkeypatch.setattr(pickle, "loads", forbidden)
    monkeypatch.setattr(rasterio, "open", forbidden)
    monkeypatch.setattr(Image, "open", forbidden)
    monkeypatch.setattr(resolver, "fetch_copernicus_dem", forbidden)
    args = parse_replay_args(
        ["35,133", "-B", "--list-cache-ranges", "--cache-dir", str(tmp_path)]
    )
    result = query_cache_ranges(args)
    assert result["layers"]["cloud-native"]["record_count"] == 1
    assert result["verification_level"] == "metadata-only"
    assert not result["payloads_verified"]
    assert not result["complete_output_ready"]
    assert result["ephemeris_ranges"] is None
    assert result["celestial_common_ranges"] is None
    assert result["location"]["ground_elevation_m"] is None
    assert result["location"]["observer_height_m"] is None
    assert "Metadata only" in format_report(result)


@pytest.mark.parametrize("layer", ["cloud-native", "cloud-geo"])
def test_cloud_hold_uses_acquisition_order_and_replaces_with_next_success(layer):
    first = record(layer, source=T - 26 * MINUTE, end=T + 30 * MINUTE)
    second = record(
        layer, start=T + 42 * MINUTE, source=T - 30 * MINUTE, end=T + 72 * MINUTE
    )
    index = HistoryIndex([second, first])
    assert index.select(T - MINUTE) is None
    assert index.select(T + 41 * MINUTE) is first
    # Even an older observation replaces its predecessor when acquired later,
    # matching the GUI's accepted-source replacement semantics.
    assert index.select(T + 42 * MINUTE) is second
    assert index.select(T + 72 * MINUTE) is second
    assert index.select(T + 73 * MINUTE) is None
    assert coverage(index, T, T + 80 * MINUTE) == [(T, T + 73 * MINUTE)]


def test_native_cloud_tail_is_measured_from_fetch_not_observation(tmp_path):
    path = tmp_path / "source_history/2026/10/14/00/26/source.json"
    artifact = path.with_suffix(".pkl")
    artifact.parent.mkdir(parents=True)
    artifact.write_bytes(b"metadata only")
    write_json(
        path,
        dict(
            schema_version=1,
            source_time_utc=T.isoformat(),
            fetched_at_utc=(T + 26 * MINUTE).isoformat(),
            observer_latitude_deg=35,
            observer_longitude_deg=133,
            artifact_path=str(artifact.relative_to(tmp_path)),
            artifact_sha256="a" * 64,
        ),
    )
    inventory = Inventory(tmp_path, 35, 133, 1.7)
    inventory.scan()
    index = HistoryIndex(inventory.records["cloud-native"])
    assert index.select(T + 25 * MINUTE) is None
    assert index.select(T + 40 * MINUTE) is not None
    assert index.select(T + 56 * MINUTE) is not None
    assert index.select(T + 57 * MINUTE) is None


def test_geo_tail_is_measured_from_fetch(tmp_path):
    geo(tmp_path, fetched=T + 22 * MINUTE)
    inventory = Inventory(tmp_path, 51, 0, 1.7)
    inventory.scan()
    index = HistoryIndex(inventory.records["cloud-geo"])
    assert index.select(T + 21 * MINUTE) is None
    assert index.select(T + 82 * MINUTE) is not None
    assert index.select(T + 83 * MINUTE) is None
