"""Cache contract tests for the low-level star renderer.

Contract: a `StarRenderCache` entry is reused only when every input that can
change the rendered pixels is unchanged.  A cached result must therefore be
pixel-identical to a fresh render of the same inputs.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace

import astropy.time
import numpy as np
import pytest
from PySide6.QtGui import QImage, QPainter

from zstarview.render import stars as render_stars
from zstarview.render.qt_image import qimage_to_np_rgba
from zstarview.types import CelestialData, ScreenGeometry, ViewerData

_SIZE = (180, 140)


def _stars(seed: int = 7, count: int = 160) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    vmag = np.concatenate([[-1.0, 0.5, 1.9, 2.0, 4.0], rng.uniform(2.1, 7.0, count - 5)])
    return {
        "star_index": np.arange(count, dtype=np.int32) * 3 + 11,
        "alt": rng.uniform(5.0, 85.0, count),
        "az": rng.uniform(120.0, 240.0, count),
        "vmag": vmag,
        "bv": rng.uniform(-0.3, 1.8, count),
        # Spans 1 px dots, 2 px squares, 3..6 px fills and >= 7 px outlines.
        "size_factor": np.clip(10.0 ** (-0.15 * (vmag - 1.0)), 0.05, 4.0),
        "color_factor_base": rng.uniform(0.1, 1.0, count),
    }


def _data(stars: dict[str, np.ndarray]) -> CelestialData:
    return CelestialData(
        time=astropy.time.Time("2026-03-20T00:00:00", scale="utc"),
        planets=[],
        stars=stars,  # type: ignore[arg-type]
        deep_sky_objects={},  # type: ignore[arg-type]
        celestial_equator_points=[],
        ecliptic_points=[],
        horizon_points=[],
    )


def _viewer(**changes: object) -> ViewerData:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Test",
        view_center=(45.0, 180.0),
        edge_fov_deg=95.0,
        content_fov_deg=115.0,
    )
    return replace(viewer, **changes) if changes else viewer  # type: ignore[arg-type]


def _geometry(size: tuple[int, int]) -> ScreenGeometry:
    return ScreenGeometry(center=(size[0] // 2, size[1] // 2), radius=min(size) // 2 - 4)


def _render(
    stars: dict[str, np.ndarray],
    *,
    cache: render_stars.StarRenderCache | None,
    viewer: ViewerData | None = None,
    size: tuple[int, int] = _SIZE,
    geometry: ScreenGeometry | None = None,
    **kwargs: object,
) -> np.ndarray:
    image = QImage(size[0], size[1], QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(0)
    painter = QPainter(image)
    try:
        options: dict[str, object] = {
            "star_base_radius": 4.0,
            "draw_vmag_limit": 6.5,
            "viewport_size": size,
            "render_cache": cache,
        }
        options.update(kwargs)
        draw = render_stars.draw_stars_fast if options.pop("fast_mode", False) else render_stars.draw_stars_normal
        draw(
            painter,
            geometry or _geometry(size),
            _data(stars),
            viewer or _viewer(),
            **options,  # type: ignore[arg-type]
        )
    finally:
        painter.end()
    return qimage_to_np_rgba(image)


def _with(name: str, update: Callable[[np.ndarray], np.ndarray]) -> Callable[[dict], dict]:
    def apply(stars: dict) -> dict:
        changed = dict(stars)
        changed[name] = update(stars[name].copy())
        return changed

    return apply


_COLOR_FACTOR_CHANGE = _with("color_factor_base", lambda a: a * 0.5)

# Each case: (id, star-table change, viewer change, call option change).
_INVALIDATING_CASES: list[tuple[str, Callable[[dict], dict] | None, dict, dict]] = [
    ("view-center", None, {"view_center": (45.0, 183.0)}, {}),
    ("edge-fov", None, {"edge_fov_deg": 70.0}, {}),
    ("viewport-size", None, {}, {"size": (_SIZE[0] + 20, _SIZE[1]), "geometry": _geometry(_SIZE)}),
    ("geometry", None, {}, {"geometry": ScreenGeometry(center=(80, 70), radius=60)}),
    ("vmag-limit", None, {}, {"draw_vmag_limit": 4.5}),
    ("vmag-min-exclusive", None, {}, {"draw_vmag_min_exclusive": 4.0}),
    ("star-base-radius", None, {}, {"star_base_radius": 5.0}),
    ("visibility-boost", None, {}, {"visibility_boost": 1.6}),
    ("outline", None, {}, {"outline_bright_bodies": True, "outline_render_scale": 1.5}),
    ("fast-mode", None, {}, {"fast_mode": True}),
    ("twinkle", None, {}, {"twinkle_targets": ((11, 0.8),)}),
    ("alt", _with("alt", lambda a: a + 1.0), {}, {}),
    ("az", _with("az", lambda a: a + 1.0), {}, {}),
    ("bv", _with("bv", lambda a: a + 0.8), {}, {}),
    ("vmag", _with("vmag", lambda a: a - 0.4), {}, {}),
    ("size-factor", _with("size_factor", lambda a: a * 1.4), {}, {}),
    ("color-factor", _COLOR_FACTOR_CHANGE, {}, {}),
    ("star-index", _with("star_index", lambda a: a[::-1].copy()), {}, {"twinkle_targets": ((11, 0.8),)}),
]


@pytest.mark.parametrize("light", [False, True], ids=["dark", "light"])
@pytest.mark.parametrize(
    ("change_stars", "viewer_change", "option_change"),
    [case[1:] for case in _INVALIDATING_CASES],
    ids=[case[0] for case in _INVALIDATING_CASES],
)
def test_star_cache_never_returns_stale_image_after_input_change(
    light: bool,
    change_stars: Callable[[dict], dict] | None,
    viewer_change: dict,
    option_change: dict,
) -> None:
    stars = _stars()
    base_options: dict[str, object] = {"light_background_outline": light}
    if "twinkle_targets" in option_change:
        base_options["twinkle_targets"] = ()
    changed_stars = change_stars(stars) if change_stars is not None else stars
    changed_viewer = _viewer(**viewer_change)
    changed_options = {**base_options, **option_change}

    cache = render_stars.StarRenderCache()
    base = _render(stars, cache=cache, **base_options)  # type: ignore[arg-type]
    assert cache.entry is not None
    cached = _render(changed_stars, cache=cache, viewer=changed_viewer, **changed_options)  # type: ignore[arg-type]
    fresh = _render(changed_stars, cache=None, viewer=changed_viewer, **changed_options)  # type: ignore[arg-type]

    np.testing.assert_array_equal(cached, fresh)
    light_ignored = {"twinkle_targets", "fast_mode", "visibility_boost"}
    if not (light and (light_ignored & option_change.keys() or change_stars is _COLOR_FACTOR_CHANGE)):
        # The light-background renderer ignores twinkle, fast mode, and the
        # dark-sky brightness terms (visibility boost and color factor).
        assert not np.array_equal(base, fresh), "test case must change the image"


@pytest.mark.parametrize("light", [False, True], ids=["dark", "light"])
def test_star_cache_reuses_image_for_identical_inputs(light: bool) -> None:
    stars = _stars()
    cache = render_stars.StarRenderCache()
    first = _render(stars, cache=cache, light_background_outline=light)
    entry = cache.entry
    assert entry is not None
    # Equal content in new array objects must still hit the cache.
    copied = {name: values.copy() for name, values in stars.items()}
    second = _render(copied, cache=cache, light_background_outline=light)
    assert cache.entry is entry
    np.testing.assert_array_equal(first, second)


def test_star_cache_detects_in_place_array_mutation() -> None:
    stars = _stars()
    cache = render_stars.StarRenderCache()
    _render(stars, cache=cache)
    stars["az"] += 2.0
    cached = _render(stars, cache=cache)
    fresh = _render(stars, cache=None)
    np.testing.assert_array_equal(cached, fresh)


def test_star_cache_owners_are_independent() -> None:
    stars = _stars()
    first_cache = render_stars.StarRenderCache()
    second_cache = render_stars.StarRenderCache()
    _render(stars, cache=first_cache)
    moved = _viewer(view_center=(45.0, 186.0))
    second = _render(stars, cache=second_cache, viewer=moved)
    again = _render(stars, cache=first_cache)
    np.testing.assert_array_equal(second, _render(stars, cache=None, viewer=moved))
    np.testing.assert_array_equal(again, _render(stars, cache=None))


def test_star_cache_is_not_reused_or_poisoned_by_screen_positions() -> None:
    stars = _stars()
    viewer = _viewer()
    mask = stars["vmag"] <= 6.5
    plain = _render(stars, cache=None)
    cache = render_stars.StarRenderCache()
    _render(stars, cache=cache)
    positions = np.column_stack(
        [np.linspace(20.0, 160.0, int(mask.sum())), np.linspace(15.0, 125.0, int(mask.sum()))]
    )
    with_positions = _render(stars, cache=cache, viewer=viewer, screen_positions=positions)
    np.testing.assert_array_equal(
        with_positions, _render(stars, cache=None, viewer=viewer, screen_positions=positions)
    )
    assert not np.array_equal(with_positions, plain)
    # A later call without explicit positions must not reuse the custom image.
    np.testing.assert_array_equal(_render(stars, cache=cache), plain)


def test_star_cache_vmag_boundaries_are_inclusive_limit_and_exclusive_minimum() -> None:
    stars = _stars()
    # Index 4 is exactly vmag 4.0, index 3 exactly 2.0.
    only_four = {name: values[4:5].copy() for name, values in stars.items()}
    only_four["alt"][:] = 45.0
    only_four["az"][:] = 180.0
    cache = render_stars.StarRenderCache()
    included = _render(only_four, cache=cache, draw_vmag_limit=4.0)
    assert included[:, :, 3].any()
    excluded = _render(only_four, cache=cache, draw_vmag_limit=6.5, draw_vmag_min_exclusive=4.0)
    assert not excluded[:, :, 3].any()


def test_star_cache_after_resize_has_no_stale_pixels() -> None:
    stars = _stars()
    cache = render_stars.StarRenderCache()
    for size in [(180, 140), (240, 160), (180, 140), (120, 200)]:
        cached = _render(stars, cache=cache, size=size)
        np.testing.assert_array_equal(cached, _render(stars, cache=None, size=size))


def _frozen(stars: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    frozen = {name: values.copy() for name, values in stars.items()}
    for values in frozen.values():
        values.flags.writeable = False
    return frozen


def test_star_cache_reuses_frozen_snapshot_by_identity() -> None:
    stars = _frozen(_stars())
    cache = render_stars.StarRenderCache()
    first = _render(stars, cache=cache)
    entry = cache.entry
    second = _render(dict(stars), cache=cache)
    assert cache.entry is entry
    np.testing.assert_array_equal(first, second)


def test_star_cache_treats_new_frozen_snapshot_as_new_data() -> None:
    stars = _frozen(_stars())
    cache = render_stars.StarRenderCache()
    _render(stars, cache=cache)
    moved = dict(stars)
    shifted = stars["az"] + 2.0
    shifted.flags.writeable = False
    moved["az"] = shifted
    np.testing.assert_array_equal(_render(moved, cache=cache), _render(moved, cache=None))


def test_star_cache_hashes_read_only_views() -> None:
    base = _stars()
    view_stars = dict(base)
    az_view = base["az"][:]
    az_view.flags.writeable = False
    view_stars["az"] = az_view
    cache = render_stars.StarRenderCache()
    _render(view_stars, cache=cache)
    base["az"] += 2.0  # mutates the view's base in place
    np.testing.assert_array_equal(_render(view_stars, cache=cache), _render(view_stars, cache=None))


def test_calculate_visible_stars_returns_frozen_snapshot() -> None:
    from zstarview.astro import calculate_visible_stars

    catalog = {
        "catalog_index": np.arange(3, dtype=np.int32),
        "unit_vectors": np.eye(3),
        "vmag": np.array([1.0, 3.0, 5.0]),
        "bv": np.zeros(3),
        "size_scale": np.ones(3),
        "color_base": np.ones(3),
    }
    time = astropy.time.Time("2026-03-20T00:00:00", scale="utc")
    for kwargs in ({}, {"max_vmag": 0.0}, {"subset_indices": np.array([], dtype=np.int32)}):
        stars, _location = calculate_visible_stars(catalog, 35.0, 139.0, 0.0, time, **kwargs)  # type: ignore[arg-type]
        for name, values in stars.items():
            assert not values.flags.writeable, name
