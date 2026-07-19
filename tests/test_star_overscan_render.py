from __future__ import annotations

import astropy.time
import numpy as np
from PySide6.QtGui import QImage, QPainter

from zstarview.render import stars as render_stars
from zstarview.render.qt_image import qimage_to_np_rgba
from zstarview.types import CelestialData, ScreenGeometry, StarCatalogMeta, ViewerData


def _single_star_celestial_data(
    *,
    alt: float,
    az: float,
    vmag: float = 5.5,
    bv: float = 0.0,
    size_factor: float = 0.18,
    color_factor_base: float = 0.7,
) -> CelestialData:
    return CelestialData(
        time=astropy.time.Time("2026-03-20T00:00:00", scale="utc"),
        planets=[],
        stars={
            "star_index": np.array([0], dtype=np.int32),
            "alt": np.array([alt], dtype=float),
            "az": np.array([az], dtype=float),
            "vmag": np.array([vmag], dtype=float),
            "bv": np.array([bv], dtype=float),
            "size_factor": np.array([size_factor], dtype=float),
            "color_factor_base": np.array([color_factor_base], dtype=float),
        },
        deep_sky_objects={
            "id": np.array([], dtype=object),
            "name": np.array([], dtype=object),
            "type": np.array([], dtype=object),
            "alt": np.array([], dtype=float),
            "az": np.array([], dtype=float),
            "vmag": np.array([], dtype=float),
            "major_arcmin": np.array([], dtype=float),
            "minor_arcmin": np.array([], dtype=float),
            "pa_deg": np.array([], dtype=float),
        },
        celestial_equator_points=[],
        ecliptic_points=[],
        horizon_points=[],
        star_catalog_meta=StarCatalogMeta(
            name_indices=np.array([0], dtype=np.int32),
            names=np.array(["Faint"], dtype=object),
            source_id_indices=np.array([0], dtype=np.int32),
            source_ids=np.array(["HIP_FAINT"], dtype=object),
        ),
    )


def _alpha_bbox_width(image: QImage) -> int:
    arr = qimage_to_np_rgba(image)
    xs = np.nonzero(arr[:, :, 3])[1]
    if xs.size == 0:
        return 0
    return int(xs.max() - xs.min() + 1)


def test_draw_stars_keeps_faint_overscan_star_outside_90_deg_background() -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(120, 120), radius=80)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            edge_fov_deg=90.0,
            content_fov_deg=110.0,
        )
        celestial_data = _single_star_celestial_data(
            # Angular distance is 95 deg from the view center, so this is outside 90 deg
            # but still within the requested 110 deg content FOV.
            alt=-50.0,
            az=180.0,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=4.0,
            viewport_size=(240, 240),
            content_fov_deg=110.0,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    expected_y = round(geometry.center[1] + geometry.radius * 95.0 / viewer.edge_fov_deg)
    assert int(arr[expected_y, geometry.center[0], 3]) > 0


def test_light_background_star_render_skips_subpixel_stars() -> None:
    image = QImage(120, 120, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0xFFFFFFFF)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(60, 60), radius=50)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            content_fov_deg=110.0,
        )
        celestial_data = _single_star_celestial_data(
            alt=45.0,
            az=180.0,
            size_factor=0.1,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=4.0,
            viewport_size=(120, 120),
            light_background_outline=True,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    assert np.all(arr[:, :, :3] == 255)


def test_scenic_bright_star_underlay_includes_fourth_magnitude_only() -> None:
    def render_underlay(vmag: float) -> np.ndarray:
        image = QImage(120, 120, QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(0xFF808080)
        painter = QPainter(image)
        try:
            geometry = ScreenGeometry(center=(60, 60), radius=50)
            viewer = ViewerData(
                location=(35.0, 139.0),
                timezone_name="UTC",
                city_name="Tokyo",
                view_center=(90.0, 180.0),
                content_fov_deg=90.0,
            )
            render_stars.draw_bright_star_underlay(
                painter,
                geometry,
                _single_star_celestial_data(
                    alt=90.0,
                    az=180.0,
                    vmag=vmag,
                    size_factor=0.5,
                ),
                viewer,
                star_base_radius=12.0,
                outline_bright_bodies=False,
                viewport_size=(120, 120),
                content_fov_deg=90.0,
            )
        finally:
            painter.end()
        return qimage_to_np_rgba(image)

    assert np.any(render_underlay(4.0)[55:66, 55:66, :3] < 128)
    assert np.all(render_underlay(4.0)[60, 60, :3] == 128)
    assert np.all(render_underlay(4.01)[55:66, 55:66, :3] == 128)


def test_light_background_star_render_draws_outline_before_body() -> None:
    image = QImage(120, 120, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0xFFFFFFFF)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(60, 60), radius=50)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(45.0, 180.0),
            content_fov_deg=110.0,
        )
        celestial_data = _single_star_celestial_data(
            alt=45.0,
            az=180.0,
            bv=1.0,
            size_factor=0.5,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=4.0,
            viewport_size=(120, 120),
            light_background_outline=True,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    center = arr[60, 60, :3]
    neighborhood = arr[58:63, 58:63, :3].reshape(-1, 3)

    assert np.any(np.max(neighborhood, axis=1) < 200)
    assert not np.max(center) < 100
    assert not np.all(center == (255, 255, 255))


def test_draw_stars_uses_peak_channel_as_alpha_for_faint_pixels() -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(120, 120), radius=80)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(90.0, 180.0),
            content_fov_deg=90.0,
        )
        celestial_data = _single_star_celestial_data(
            alt=90.0,
            az=180.0,
            vmag=5.8,
            bv=0.35,
            size_factor=0.18,
            color_factor_base=0.18,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=4.0,
            viewport_size=(240, 240),
            content_fov_deg=90.0,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    center = arr[120, 120, :]

    assert 0 < int(center[3]) < 255
    assert int(np.max(center[:3])) >= 250
    assert int(np.min(center[:3])) > 0


def test_draw_stars_renders_bright_outline_only_when_requested() -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(120, 120), radius=80)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(90.0, 180.0),
            content_fov_deg=90.0,
        )
        celestial_data = _single_star_celestial_data(
            alt=90.0,
            az=180.0,
            vmag=1.0,
            bv=0.0,
            size_factor=3.0,
            color_factor_base=1.0,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=20.0,
            outline_bright_bodies=True,
            viewport_size=(240, 240),
            content_fov_deg=90.0,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    center = arr[120, 120, :]

    assert int(center[3]) == 0
    assert int(np.count_nonzero(arr[:, :, 3])) > 0


def test_draw_stars_renders_7px_and_larger_stars_as_outline_rectangles() -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(120, 120), radius=80)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(90.0, 180.0),
            content_fov_deg=90.0,
        )
        celestial_data = _single_star_celestial_data(
            alt=90.0,
            az=180.0,
            vmag=5.5,
            bv=0.0,
            size_factor=0.27,
            color_factor_base=1.0,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=20.0,
            viewport_size=(240, 240),
            content_fov_deg=90.0,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    center = arr[120, 120, :]
    top_edge = arr[117, 120, :]

    assert int(center[3]) == 0
    assert int(top_edge[3]) > 0


def test_draw_stars_fast_mode_renders_7px_and_larger_stars_as_filled_rectangles() -> None:
    image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(120, 120), radius=80)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(90.0, 180.0),
            content_fov_deg=90.0,
        )
        celestial_data = _single_star_celestial_data(
            alt=90.0,
            az=180.0,
            vmag=5.5,
            bv=0.0,
            size_factor=0.27,
            color_factor_base=1.0,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=20.0,
            viewport_size=(240, 240),
            content_fov_deg=90.0,
            fast_mode=True,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    center = arr[120, 120, :]
    top_edge = arr[117, 120, :]

    assert int(center[3]) > 0
    assert int(top_edge[3]) > 0
    assert int(np.count_nonzero(arr[:, :, 3])) > 0


def test_draw_stars_keeps_bright_diamonds_no_smaller_than_outline_rectangles_at_high_scale() -> None:
    geometry = ScreenGeometry(center=(120, 120), radius=80)
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        view_center=(90.0, 180.0),
        content_fov_deg=90.0,
    )

    bright_image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    bright_image.fill(0)
    bright_painter = QPainter(bright_image)
    try:
        bright_data = _single_star_celestial_data(
            alt=90.0,
            az=180.0,
            vmag=1.9,
            bv=0.0,
            size_factor=0.27,
            color_factor_base=1.0,
        )
        render_stars.draw_stars(
            bright_painter,
            geometry,
            bright_data,
            viewer,
            star_base_radius=20.0,
            outline_bright_bodies=True,
            outline_render_scale=2.5,
            viewport_size=(240, 240),
            content_fov_deg=90.0,
        )
    finally:
        bright_painter.end()

    faint_image = QImage(240, 240, QImage.Format.Format_ARGB32_Premultiplied)
    faint_image.fill(0)
    faint_painter = QPainter(faint_image)
    try:
        faint_data = _single_star_celestial_data(
            alt=90.0,
            az=180.0,
            vmag=5.5,
            bv=0.0,
            size_factor=0.27,
            color_factor_base=1.0,
        )
        render_stars.draw_stars(
            faint_painter,
            geometry,
            faint_data,
            viewer,
            star_base_radius=20.0,
            outline_bright_bodies=True,
            outline_render_scale=2.5,
            viewport_size=(240, 240),
            content_fov_deg=90.0,
        )
    finally:
        faint_painter.end()

    assert _alpha_bbox_width(bright_image) >= _alpha_bbox_width(faint_image)


def test_light_background_bright_outline_uses_diamond_underlay_and_marker() -> None:
    image = QImage(120, 120, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0xFFFFFFFF)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(60, 60), radius=50)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(90.0, 180.0),
            content_fov_deg=90.0,
        )
        celestial_data = _single_star_celestial_data(
            alt=90.0,
            az=180.0,
            vmag=1.0,
            bv=0.0,
            size_factor=0.5,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=12.0,
            outline_bright_bodies=True,
            viewport_size=(120, 120),
            content_fov_deg=90.0,
            light_background_outline=True,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    neighborhood = arr[52:69, 52:69, :3].reshape(-1, 3)
    assert np.all(arr[56, 56, :3] == 255)  # diagonal corner is not a square outline
    assert np.any(np.max(neighborhood, axis=1) < 200)  # dark outer diamond underlay
    assert np.any(
        (np.min(neighborhood, axis=1) < 250) & (np.max(neighborhood, axis=1) > 100)
    )
    # colored inner diamond marker
    assert np.all(arr[60, 60, :3] == 255)  # outline mode leaves the center open


def test_light_background_bright_fill_uses_filled_diamond() -> None:
    image = QImage(120, 120, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0xFFFFFFFF)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(60, 60), radius=50)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(90.0, 180.0),
            content_fov_deg=90.0,
        )
        celestial_data = _single_star_celestial_data(
            alt=90.0,
            az=180.0,
            vmag=1.0,
            bv=0.0,
            size_factor=0.5,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=12.0,
            outline_bright_bodies=False,
            viewport_size=(120, 120),
            content_fov_deg=90.0,
            light_background_outline=True,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    neighborhood = arr[52:69, 52:69, :3].reshape(-1, 3)
    assert not np.all(arr[60, 60, :3] == 255)
    assert np.all(arr[56, 56, :3] == 255)  # diamond corners remain outside the fill
    assert np.any(np.max(neighborhood, axis=1) < 200)  # dark underlay remains visible


def test_light_background_magnitude_boundary_keeps_square_rendering_at_two() -> None:
    image = QImage(120, 120, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0xFFFFFFFF)
    painter = QPainter(image)
    try:
        geometry = ScreenGeometry(center=(60, 60), radius=50)
        viewer = ViewerData(
            location=(35.0, 139.0),
            timezone_name="UTC",
            city_name="Tokyo",
            view_center=(90.0, 180.0),
            content_fov_deg=90.0,
        )
        celestial_data = _single_star_celestial_data(
            alt=90.0,
            az=180.0,
            vmag=2.0,
            bv=0.0,
            size_factor=0.5,
        )

        render_stars.draw_stars(
            painter,
            geometry,
            celestial_data,
            viewer,
            star_base_radius=12.0,
            outline_bright_bodies=True,
            viewport_size=(120, 120),
            content_fov_deg=90.0,
            light_background_outline=True,
        )
    finally:
        painter.end()

    arr = qimage_to_np_rgba(image)
    assert np.all(arr[56, 56, :3] == 255)  # vmag == 2.0 uses the diamond path
