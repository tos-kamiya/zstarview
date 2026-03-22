from __future__ import annotations

import astropy.time

import zstarview.satellites.project as project_module


class _FakeAngle:
    def __init__(self, degrees: float) -> None:
        self.degrees = degrees


class _FakeTopocentric:
    def __init__(self, alt_deg: float, az_deg: float) -> None:
        self._alt_deg = alt_deg
        self._az_deg = az_deg

    def altaz(self):
        return _FakeAngle(self._alt_deg), _FakeAngle(self._az_deg), object()


class _FakeDifference:
    def __init__(self, alt_deg: float, az_deg: float) -> None:
        self._alt_deg = alt_deg
        self._az_deg = az_deg

    def at(self, _time_obj):
        return _FakeTopocentric(self._alt_deg, self._az_deg)


class _FakeSatellite:
    def __init__(self, name: str, alt_deg: float, az_deg: float) -> None:
        self.name = name
        self._alt_deg = alt_deg
        self._az_deg = az_deg

    def __sub__(self, _observer):
        return _FakeDifference(self._alt_deg, self._az_deg)


def test_project_satellite_records_marks_iss_with_label(monkeypatch) -> None:
    monkeypatch.setattr(
        project_module,
        "build_earth_satellites",
        lambda records, *, ts=None: [_FakeSatellite("ISS (ZARYA)", 45.0, 120.0)],
    )

    points = project_module.project_satellite_records(
        {"iss": [{"OBJECT_NAME": "ISS (ZARYA)"}]},
        observer_lat=35.47,
        observer_lon=133.05,
        observer_height_m=0.0,
        time_obj=astropy.time.Time("2026-03-22T12:00:00Z"),
    )

    assert len(points) == 1
    assert points[0].group_key == "iss"
    assert points[0].show_label is True
    assert points[0].marker_scale == 0.3


def test_project_satellite_records_marks_css_tianhe_like_iss(monkeypatch) -> None:
    monkeypatch.setattr(
        project_module,
        "build_earth_satellites",
        lambda records, *, ts=None: [_FakeSatellite("CSS (TIANHE)", 42.0, 130.0)],
    )

    points = project_module.project_satellite_records(
        {"iss": [{"OBJECT_NAME": "CSS (TIANHE)"}]},
        observer_lat=35.47,
        observer_lon=133.05,
        observer_height_m=0.0,
        time_obj=astropy.time.Time("2026-03-22T12:00:00Z"),
    )

    assert len(points) == 1
    assert points[0].group_key == "iss"
    assert points[0].satellite_name == "CSS (TIANHE)"
    assert points[0].show_label is True
    assert points[0].marker_scale == 0.3


def test_project_satellite_records_sorts_by_group_then_altitude(monkeypatch) -> None:
    def fake_builder(records, *, ts=None):
        names = [record["OBJECT_NAME"] for record in records]
        if names == ["ISS (ZARYA)"]:
            return [_FakeSatellite("ISS (ZARYA)", 20.0, 100.0)]
        return [
            _FakeSatellite("STARLINK A", 10.0, 150.0),
            _FakeSatellite("STARLINK B", 30.0, 160.0),
        ]

    monkeypatch.setattr(project_module, "build_earth_satellites", fake_builder)

    points = project_module.project_satellite_records(
        {
            "starlink": [{"OBJECT_NAME": "STARLINK GROUP"}],
            "iss": [{"OBJECT_NAME": "ISS (ZARYA)"}],
        },
        observer_lat=35.47,
        observer_lon=133.05,
        observer_height_m=0.0,
        time_obj=astropy.time.Time("2026-03-22T12:00:00Z"),
    )

    assert [point.group_key for point in points] == ["iss", "starlink", "starlink"]
    assert [point.satellite_name for point in points[1:]] == ["STARLINK B", "STARLINK A"]
    assert all(point.show_label is False for point in points[1:])


def test_project_satellite_records_keeps_starlink_smaller_and_unlabeled(monkeypatch) -> None:
    def fake_builder(records, *, ts=None):
        names = [record["OBJECT_NAME"] for record in records]
        if names == ["ISS (ZARYA)"]:
            return [_FakeSatellite("ISS (ZARYA)", 50.0, 120.0)]
        return [_FakeSatellite("STARLINK-1", 55.0, 140.0)]

    monkeypatch.setattr(project_module, "build_earth_satellites", fake_builder)

    points = project_module.project_satellite_records(
        {
            "iss": [{"OBJECT_NAME": "ISS (ZARYA)"}],
            "starlink": [{"OBJECT_NAME": "STARLINK-1"}],
        },
        observer_lat=35.47,
        observer_lon=133.05,
        observer_height_m=0.0,
        time_obj=astropy.time.Time("2026-03-22T12:00:00Z"),
    )

    assert [point.group_key for point in points] == ["iss", "starlink"]
    assert points[0].show_label is True
    assert points[1].show_label is False
    assert points[1].marker_scale == 0.156
    assert points[1].marker_scale < points[0].marker_scale
