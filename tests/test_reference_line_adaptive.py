from astropy.time import Time

from zstarview.render.guides import draw_sky_reference_lines
from zstarview.types import CelestialData, ScreenGeometry


class _PolylineCountingPainter:
    def __init__(self) -> None:
        self.polylines: list[list[tuple[float, float]]] = []

    def save(self) -> None:
        pass

    def restore(self) -> None:
        pass

    def setPen(self, *_args, **_kwargs) -> None:
        pass

    def drawPolyline(self, polyline) -> None:
        self.polylines.append([(point.x(), point.y()) for point in polyline])


def _make_celestial_data() -> CelestialData:
    return CelestialData(
        time=Time("2026-04-27T00:00:00", scale="utc"),
        planets=[],
        stars={
            "star_index": [],
            "alt": [],
            "az": [],
            "vmag": [],
            "bv": [],
            "size_factor": [],
            "color_factor_base": [],
        },
        deep_sky_objects={
            "id": [],
            "name": [],
            "type": [],
            "alt": [],
            "az": [],
            "vmag": [],
            "major_arcmin": [],
            "minor_arcmin": [],
            "pa_deg": [],
        },
        celestial_equator_points=[(0.0, 0.0), (90.0, 90.0)],
        ecliptic_points=[],
        horizon_points=[],
    )


def test_reference_line_fixed_sampling_keeps_polyline_length_stable() -> None:
    def fake_projection(alt_deg: float, az_deg: float, view_center, *, edge_fov_deg: float = 95.0):
        _ = view_center, edge_fov_deg
        return (0.0001 * az_deg, 0.012 * ((alt_deg / 90.0) ** 2))

    small_painter = _PolylineCountingPainter()
    large_painter = _PolylineCountingPainter()
    draw_sky_reference_lines(
        small_painter,  # type: ignore[arg-type]
        geometry=ScreenGeometry(center=(100, 100), radius=50),
        celestial_data=_make_celestial_data(),
        viewer_data=type(
            "Viewer",
            (),
            {
                "view_center": (0.0, 0.0),
                "content_fov_deg": 180.0,
                "edge_fov_deg": 95.0,
            },
        )(),
        altaz_to_normalized_xy_func=fake_projection,
    )
    draw_sky_reference_lines(
        large_painter,  # type: ignore[arg-type]
        geometry=ScreenGeometry(center=(400, 400), radius=400),
        celestial_data=_make_celestial_data(),
        viewer_data=type(
            "Viewer",
            (),
            {
                "view_center": (0.0, 0.0),
                "content_fov_deg": 180.0,
                "edge_fov_deg": 95.0,
            },
        )(),
        altaz_to_normalized_xy_func=fake_projection,
    )

    assert all(len(polyline) == 2 for polyline in small_painter.polylines)
    assert all(len(polyline) == 2 for polyline in large_painter.polylines)
    assert len(small_painter.polylines) == len(large_painter.polylines)
