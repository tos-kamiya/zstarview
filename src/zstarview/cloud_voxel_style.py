"""Application-level appearance settings for voxel-rendered clouds."""

from __future__ import annotations

import itertools
import math
from contextlib import AbstractContextManager, nullcontext
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

# Ambient night cloud color at full cloud amount, in normalized RGB.
CLOUD_VOXEL_NIGHT_COLOR_RGB = (0.36, 0.385, 0.44)

# Voxel cloud images use a 513x513 source surface (256px radius).
CLOUD_VOXEL_CUTOUT_REFERENCE_SIZE = 513
CLOUD_VOXEL_CUTOUT_WIDTH_PX = 1
CLOUD_VOXEL_CUTOUT_COUNT = 16

# Clear holes are drawn at the vertices of a truncated icosahedron fixed in
# alt/az.  The radius is in pixels of the 513x513 reference surface, so the
# hole size on screen does not change with the field of view.
CLOUD_VOXEL_HOLE_RADIUS_PX = 4

CLOUD_VOXEL_MARK_HOLES = "holes"
CLOUD_VOXEL_MARK_CUTOUT = "cutout"
CLOUD_VOXEL_MARK_MODES = (CLOUD_VOXEL_MARK_HOLES, CLOUD_VOXEL_MARK_CUTOUT)


@dataclass(frozen=True)
class CloudVoxelMarkStyle:
    """Transparent marks that make thin cloud layers easier to recognize.

    ``value`` is the hole radius in reference pixels for ``holes`` and the
    number of diagonal lines for ``cutout``.  ``0`` disables the marks.
    """

    mode: str = CLOUD_VOXEL_MARK_HOLES
    value: int = CLOUD_VOXEL_HOLE_RADIUS_PX

    @classmethod
    def from_cli(cls, cloud_stripe: tuple[str, int, float]) -> CloudVoxelMarkStyle:
        """Build a style from the parsed ``--cloud-stripe`` value."""
        mode, value, _ = cloud_stripe
        return cls(mode=str(mode), value=max(0, int(value)))


DEFAULT_CLOUD_VOXEL_MARK_STYLE = CloudVoxelMarkStyle()


def render_cloud_voxel_image(
    source,
    *,
    lat: float,
    lon: float,
    view_center: tuple[float, float],
    edge_fov_deg: float,
    content_fov_deg: float,
    radius_px: int,
    sun_alt_deg: float,
    sun_az_deg: float,
    request_id: int,
    timeout_s: float = 120.0,
    mark_style: CloudVoxelMarkStyle = DEFAULT_CLOUD_VOXEL_MARK_STYLE,
    moon_altaz_phase: tuple[float, float, float] | None = None,
    native_work_lock: AbstractContextManager | None = None,
) -> np.ndarray:
    """Render and style the shared voxel cloud overlay image."""
    with native_work_lock if native_work_lock is not None else nullcontext():
        from .clouddisc.workers.cloud_voxel_worker import (
            render_cloud_voxels_in_subprocess,
        )
        from .night_lights import night_light_strength_factor

        sunlight_mix = 1.0 - night_light_strength_factor(float(sun_alt_deg))

    cloud_rgba = render_cloud_voxels_in_subprocess(
        source,
        lat=lat,
        lon=lon,
        view_center=view_center,
        edge_fov_deg=edge_fov_deg,
        content_fov_deg=content_fov_deg,
        radius_px=radius_px,
        sun_alt_deg=sun_alt_deg,
        sun_az_deg=sun_az_deg,
        sunlight_mix=sunlight_mix,
        night_color_rgb=CLOUD_VOXEL_NIGHT_COLOR_RGB,
        request_id=request_id,
        timeout_s=timeout_s,
        moon_altaz_phase=moon_altaz_phase,
        native_work_lock=native_work_lock,
    )
    with native_work_lock if native_work_lock is not None else nullcontext():
        return apply_cloud_voxel_marks(
            cloud_rgba,
            mark_style,
            view_center=view_center,
            content_fov_deg=content_fov_deg,
        )


def apply_cloud_voxel_marks(
    cloud_rgba: np.ndarray,
    style: CloudVoxelMarkStyle,
    *,
    view_center: tuple[float, float],
    content_fov_deg: float,
) -> np.ndarray:
    """Apply the configured transparent marks to a voxel cloud raster."""
    if style.mode == CLOUD_VOXEL_MARK_CUTOUT:
        return apply_cloud_voxel_cutout(cloud_rgba, count=style.value)
    if style.mode == CLOUD_VOXEL_MARK_HOLES:
        return apply_cloud_voxel_holes(
            cloud_rgba,
            radius_px=style.value,
            view_center=view_center,
            content_fov_deg=content_fov_deg,
        )
    raise ValueError(f"Unknown cloud mark mode: {style.mode!r}")


def apply_cloud_voxel_cutout(
    cloud_rgba: np.ndarray,
    *,
    count: int = CLOUD_VOXEL_CUTOUT_COUNT,
) -> np.ndarray:
    """Return a voxel cloud raster with sparse, 45-degree cutout lines."""
    image = np.asarray(cloud_rgba)
    if image.ndim != 3 or image.shape[2] != 4 or int(count) <= 0:
        return image

    height, width = image.shape[:2]
    if height < 3 or width < 1:
        return image

    line_count = min(int(count), height + width - 2)
    line_width = max(
        1,
        int(
            round(
                CLOUD_VOXEL_CUTOUT_WIDTH_PX
                * min(height, width)
                / CLOUD_VOXEL_CUTOUT_REFERENCE_SIZE
            )
        ),
    )
    center_x = (width - 1) * 0.5
    center_y = (height - 1) * 0.5
    x = np.arange(width, dtype=np.float32)[None, :] - center_x
    y = np.arange(height, dtype=np.float32)[:, None] - center_y
    diagonal_coord = x - y
    max_offset = 0.5 * (min(height, width) - 1) * np.sqrt(2.0)
    line_offsets = np.rint(
        np.linspace(-max_offset, max_offset, line_count + 2)[1:-1]
    )
    half_width_in_diagonal = line_width / np.sqrt(2.0)
    cutout_mask = np.zeros((height, width), dtype=bool)
    for offset in line_offsets:
        cutout_mask |= np.abs(diagonal_coord - offset) <= half_width_in_diagonal

    output = image.copy()
    output[cutout_mask, 3] = 0
    return output


@lru_cache(maxsize=1)
def cloud_hole_altaz_points() -> tuple[np.ndarray, np.ndarray]:
    """Return (alt_deg, az_deg) of the 60 truncated-icosahedron vertices.

    A pentagonal face is centred on the zenith, so the 30 vertices above the
    horizon lie on rings at about 69.9, 46.6, 31.0 and 9.9 degrees altitude
    and none sits on the horizon itself.  Neighbouring vertices are about
    23.3 degrees apart.  One vertex of the zenith pentagon faces north.
    """
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    bases = ((0.0, 1.0, 3.0 * phi), (1.0, 2.0 + phi, 2.0 * phi), (phi, 2.0, 2.0 * phi + 1.0))
    vertices: set[tuple[float, float, float]] = set()
    for base in bases:
        for signs in itertools.product((1.0, -1.0), repeat=3):
            signed = [sign * component for sign, component in zip(signs, base)]
            for shift in range(3):
                x, y, z = (round(signed[(i + shift) % 3], 12) for i in range(3))
                vertices.add((x, y, z))
    points = np.array(sorted(vertices), dtype=np.float64)
    points /= np.linalg.norm(points, axis=1)[:, None]

    # Rotate a pentagon centre (an icosahedron vertex direction) to +z.
    zenith = np.array([0.0, 1.0, phi]) / math.hypot(1.0, phi)
    north_ref = np.array([1.0, 0.0, 0.0])
    north_ref = north_ref - np.dot(north_ref, zenith) * zenith
    north_ref /= np.linalg.norm(north_ref)
    east_ref = np.cross(zenith, north_ref)
    local = points @ np.stack([north_ref, east_ref, zenith], axis=1)

    alt = np.degrees(np.arcsin(np.clip(local[:, 2], -1.0, 1.0)))
    az = np.degrees(np.arctan2(local[:, 1], local[:, 0]))
    top = np.argmax(alt)
    az = (az - az[top]) % 360.0
    alt.setflags(write=False)
    az.setflags(write=False)
    return alt, az


def project_altaz_to_cloud_image(
    alt_deg: np.ndarray,
    az_deg: np.ndarray,
    *,
    width: int,
    height: int,
    view_center: tuple[float, float],
    content_fov_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project alt/az directions to voxel cloud image pixels.

    This is the forward counterpart of ``render.ground_mask.inverse_project_disc``
    as used by the voxel worker: an azimuthal-equidistant projection whose
    image radius spans ``content_fov_deg``.  Returns ``(x, y, distance_deg)``.
    """
    eps = 1.0e-3
    phi1 = math.radians(float(np.clip(view_center[0], -90.0 + eps, 90.0 - eps)))
    lam1 = math.radians(float(view_center[1]))
    phi2 = np.radians(np.asarray(alt_deg, dtype=np.float64))
    dlam = np.radians(np.asarray(az_deg, dtype=np.float64)) - lam1

    cos_theta = np.clip(
        math.sin(phi1) * np.sin(phi2) + math.cos(phi1) * np.cos(phi2) * np.cos(dlam),
        -1.0,
        1.0,
    )
    theta_deg = np.degrees(np.arccos(cos_theta))
    psi = np.arctan2(
        np.sin(dlam) * np.cos(phi2),
        math.cos(phi1) * np.sin(phi2) - math.sin(phi1) * np.cos(phi2) * np.cos(dlam),
    )
    radius = max(1.0, (min(int(width), int(height)) - 1) * 0.5)
    r_px = theta_deg / max(1.0e-6, float(content_fov_deg)) * radius
    x = (int(width) - 1) * 0.5 + r_px * np.sin(psi)
    y = (int(height) - 1) * 0.5 - r_px * np.cos(psi)
    return x, y, theta_deg


def apply_cloud_voxel_holes(
    cloud_rgba: np.ndarray,
    *,
    radius_px: int = CLOUD_VOXEL_HOLE_RADIUS_PX,
    view_center: tuple[float, float],
    content_fov_deg: float,
) -> np.ndarray:
    """Return a voxel cloud raster with small clear holes fixed in alt/az."""
    image = np.asarray(cloud_rgba)
    if image.ndim != 3 or image.shape[2] != 4 or int(radius_px) <= 0:
        return image
    height, width = image.shape[:2]
    if height < 3 or width < 3:
        return image

    radius = max(
        1.0,
        float(radius_px) * min(height, width) / CLOUD_VOXEL_CUTOUT_REFERENCE_SIZE,
    )
    alt, az = cloud_hole_altaz_points()
    xs, ys, theta = project_altaz_to_cloud_image(
        alt,
        az,
        width=width,
        height=height,
        view_center=view_center,
        content_fov_deg=content_fov_deg,
    )
    reach = radius + 1.0
    visible = (
        (alt > 0.0)
        & (theta <= float(content_fov_deg) + 1.0)
        & (xs > -reach)
        & (xs < width - 1 + reach)
        & (ys > -reach)
        & (ys < height - 1 + reach)
    )
    if not np.any(visible):
        return image

    output = image.copy()
    alpha = output[..., 3]
    for cx, cy in zip(xs[visible], ys[visible]):
        x0 = max(0, int(math.floor(cx - reach)))
        x1 = min(width, int(math.ceil(cx + reach)) + 1)
        y0 = max(0, int(math.floor(cy - reach)))
        y1 = min(height, int(math.ceil(cy + reach)) + 1)
        if x0 >= x1 or y0 >= y1:
            continue
        gx = np.arange(x0, x1, dtype=np.float32)[None, :] - np.float32(cx)
        gy = np.arange(y0, y1, dtype=np.float32)[:, None] - np.float32(cy)
        distance = np.sqrt(gx * gx + gy * gy)
        # One-pixel anti-aliased edge: 0 inside the hole, 1 outside.
        keep = np.clip(distance - (np.float32(radius) - 0.5), 0.0, 1.0)
        window = alpha[y0:y1, x0:x1]
        alpha[y0:y1, x0:x1] = np.rint(window.astype(np.float32) * keep).astype(
            alpha.dtype
        )
    return output
