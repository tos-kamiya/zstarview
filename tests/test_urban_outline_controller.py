from __future__ import annotations

from pathlib import Path

from zstarview.gui.urban_outline_controller import UrbanOutlineController
from zstarview.types import UrbanOutlinePolyline, ViewerData


def test_run_update_keeps_base_outlines_when_skyscraper_phase_fails(monkeypatch, tmp_path: Path) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        observer_height_m=1.7,
    )
    controller = UrbanOutlineController(
        derived_root_dir=tmp_path / "normal",
        skyscraper_derived_root_dir=tmp_path / "skyscraper",
    )
    (tmp_path / "normal" / "bldg").mkdir(parents=True)
    ready_payloads = []
    controller.urban_ready.connect(lambda payload: ready_payloads.append(payload))

    monkeypatch.setattr(
        controller,
        "_required_derived_dirs",
        lambda _viewer: (("building", tmp_path / "normal" / "bldg"),),
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.resolve_urban_outline_layer_for_viewer",
        lambda *_args, **kwargs: (
            ["base-outline"]
            if kwargs.get("derived_root_dir") == controller._derived_root_dir
            else None
        ),
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.resolve_overture_release_for_cache_root",
        lambda **_kwargs: None,
    )
    monkeypatch.setattr(
        controller,
        "_selected_skyscraper_tiles",
        lambda _viewer: (
            type(
                "SeedTile",
                (),
                {
                    "zoom": 14,
                    "x": 1,
                    "y": 2,
                    "cache_key": "z14_1_2",
                    "envelope": type(
                        "Envelope",
                        (),
                        {
                            "min_lon_deg": 139.80,
                            "min_lat_deg": 35.65,
                            "max_lon_deg": 139.82,
                            "max_lat_deg": 35.67,
                        },
                    )(),
                },
            )(),
        ),
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.skyscraper_tile_derived_dir",
        lambda *_args, **_kwargs: tmp_path / "skyscraper" / "z14_1_2" / "bldg",
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.import_overture_buildings_for_bbox",
        lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
    )

    controller._run_update(viewer_data=viewer, dataset_name="tokyo-test", reason="manual")

    assert ready_payloads == [
        {
            "outlines": ["base-outline"],
            "source": "Urban: cache",
            "base_outline_count": 1,
            "skyscraper_outline_count": None,
        }
    ]


def test_run_update_skips_base_outlines_in_skyscraper_only_mode(
    monkeypatch,
    tmp_path: Path,
) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        observer_height_m=1.7,
    )
    controller = UrbanOutlineController(
        derived_root_dir=tmp_path / "normal",
        skyscraper_derived_root_dir=tmp_path / "skyscraper",
        skyscraper_only=True,
    )
    ready_payloads = []
    controller.urban_ready.connect(lambda payload: ready_payloads.append(payload))

    monkeypatch.setattr(
        controller,
        "_selected_skyscraper_tiles",
        lambda _viewer: (
            type(
                "SeedTile",
                (),
                {
                    "zoom": 14,
                    "x": 1,
                    "y": 2,
                    "cache_key": "z14_1_2",
                    "envelope": type(
                        "Envelope",
                        (),
                        {
                            "min_lon_deg": 139.80,
                            "min_lat_deg": 35.65,
                            "max_lon_deg": 139.82,
                            "max_lat_deg": 35.67,
                        },
                    )(),
                },
            )(),
        ),
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.skyscraper_tile_derived_dir",
        lambda *_args, **_kwargs: tmp_path / "skyscraper" / "z14_1_2" / "bldg",
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.import_overture_buildings_for_bbox",
        lambda **_kwargs: None,
    )

    def _resolve(*_args, **kwargs):
        if kwargs.get("derived_root_dir") == controller._derived_root_dir:
            raise AssertionError("base outlines should not be resolved in skyscraper-only mode")
        return [
            UrbanOutlinePolyline(
                points=[(180.0, 12.0), (181.0, 12.5)],
                height_m=220.0,
                distance_km=0.75,
            )
        ]

    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.resolve_urban_outline_layer_for_viewer",
        _resolve,
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.resolve_overture_release_for_cache_root",
        lambda **_kwargs: None,
    )

    controller._run_update(viewer_data=viewer, dataset_name="tokyo-test", reason="manual")

    assert ready_payloads == [
        {
            "outlines": [
                UrbanOutlinePolyline(
                    points=[(180.0, 12.0), (181.0, 12.5)],
                    height_m=220.0,
                    distance_km=0.75,
                    source="skyscraper",
                )
            ],
            "source": "Urban: cache",
            "base_outline_count": None,
            "skyscraper_outline_count": 1,
        }
    ]


def test_run_update_refreshes_stale_base_cache(monkeypatch, tmp_path: Path) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        observer_height_m=1.7,
    )
    controller = UrbanOutlineController(
        derived_root_dir=tmp_path / "normal",
        skyscraper_derived_root_dir=tmp_path / "skyscraper",
    )
    derived_dir = tmp_path / "normal" / "bldg"
    derived_dir.mkdir(parents=True)
    ready_payloads = []
    refresh_calls: list[str] = []
    controller.urban_ready.connect(lambda payload: ready_payloads.append(payload))

    monkeypatch.setattr(
        controller,
        "_required_derived_dirs",
        lambda _viewer: (("building", derived_dir),),
    )
    monkeypatch.setattr(
        controller,
        "_selected_skyscraper_tiles",
        lambda _viewer: (),
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.is_derived_dataset_stale",
        lambda path, **_kwargs: path == derived_dir,
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.import_overture_buildings",
        lambda **_kwargs: refresh_calls.append("refresh") or derived_dir,
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.resolve_urban_outline_layer_for_viewer",
        lambda *_args, **kwargs: ["outline"] if kwargs.get("derived_root_dir") == controller._derived_root_dir else None,
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.resolve_overture_release_for_cache_root",
        lambda **_kwargs: None,
    )

    controller._run_update(viewer_data=viewer, dataset_name="tokyo-test", reason="manual")

    assert refresh_calls == ["refresh"]
    assert ready_payloads == [
        {
            "outlines": ["outline"],
            "source": "Urban: cache",
            "base_outline_count": 1,
            "skyscraper_outline_count": None,
        }
    ]


def test_run_update_refreshes_when_release_changes(monkeypatch, tmp_path: Path) -> None:
    viewer = ViewerData(
        location=(35.67, 139.76),
        timezone_name="Asia/Tokyo",
        city_name="Tokyo",
        observer_height_m=1.7,
    )
    controller = UrbanOutlineController(
        derived_root_dir=tmp_path / "normal",
        skyscraper_derived_root_dir=tmp_path / "skyscraper",
    )
    derived_dir = tmp_path / "normal" / "bldg"
    derived_dir.mkdir(parents=True)
    (tmp_path / "normal" / "cache_meta.json").write_text(
        "{\"dataset_name\": \"normal\", \"fetched_at_utc\": \"2026-04-13T00:00:00+00:00\", \"overture_release\": \"2026-03-18.0\"}",
        encoding="utf-8",
    )
    ready_payloads = []
    refresh_calls: list[str] = []
    controller.urban_ready.connect(lambda payload: ready_payloads.append(payload))

    monkeypatch.setattr(
        controller,
        "_required_derived_dirs",
        lambda _viewer: (("building", derived_dir),),
    )
    monkeypatch.setattr(
        controller,
        "_selected_skyscraper_tiles",
        lambda _viewer: (),
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.resolve_overture_release_for_cache_root",
        lambda **_kwargs: "2026-04-01.0",
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.import_overture_buildings",
        lambda **_kwargs: refresh_calls.append("refresh") or derived_dir,
    )
    monkeypatch.setattr(
        "zstarview.gui.urban_outline_controller.resolve_urban_outline_layer_for_viewer",
        lambda *_args, **kwargs: ["outline"] if kwargs.get("derived_root_dir") == controller._derived_root_dir else None,
    )

    controller._run_update(viewer_data=viewer, dataset_name="tokyo-test", reason="manual")

    assert refresh_calls == ["refresh"]
    assert ready_payloads == [
        {
            "outlines": ["outline"],
            "source": "Urban: cache",
            "base_outline_count": 1,
            "skyscraper_outline_count": None,
        }
    ]
