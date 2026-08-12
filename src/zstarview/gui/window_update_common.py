from __future__ import annotations

from typing import Protocol

from ..night_lights import is_night_light_enabled


class _CloudProjectionUpdateOwner(Protocol):
    def _geo_satellite_mode_active(self) -> bool: ...

    def reproject_geo_satellite_overlay(self, reason: str = "manual") -> None: ...

    def reproject_cloud_overlay(self, reason: str = "manual") -> None: ...


class _ViewportInteractionState(Protocol):
    viewport_interaction_release_pending: bool
    viewport_interaction_completion_reason: str | None
    viewport_interaction_mode: bool
    viewport_interaction_stars: object | None


class _ViewportInteractionWaitOwner(Protocol):
    state: _ViewportInteractionState

    def _sync_viewport_interaction_chrome_visibility(self) -> None: ...

    def request_client_update(self) -> None: ...


def _request_cloud_projection_update(
    obj: _CloudProjectionUpdateOwner, *, reason: str
) -> None:
    if obj._geo_satellite_mode_active():
        obj.reproject_geo_satellite_overlay(reason=reason)
        return
    obj.reproject_cloud_overlay(reason=reason)


def _initial_data_load_active(obj: object) -> bool:
    return bool(obj._startup_initial_load_started) and not bool(
        obj._startup_initial_data_loaded
    )


def _extract_sun_altitude_deg(celestial_data: object) -> float | None:
    planets = getattr(celestial_data, "planets", None)
    if not isinstance(planets, (list, tuple)):
        return None
    for body in planets:
        if getattr(body, "name", "").strip().casefold() == "sun":
            alt = getattr(body, "alt", None)
            if isinstance(alt, (int, float)):
                return float(alt)
            try:
                return float(alt)
            except (TypeError, ValueError):
                return None
    return None


def _startup_night_light_requires_warmup(obj: object, payload: dict) -> bool:
    if float(getattr(obj, "terrain_horizon_opacity", 0.0)) <= 0.0:
        return False
    if (
        float(getattr(obj, "night_light_opacity", 0.0)) <= 0.0
        and float(getattr(obj, "ridge_glow_opacity", 0.0)) <= 0.0
    ):
        return False
    celestial_data = payload.get("celestial")
    sun_alt_deg = _extract_sun_altitude_deg(celestial_data)
    return sun_alt_deg is not None and is_night_light_enabled(float(sun_alt_deg))


def _clear_viewport_interaction_wait(obj: _ViewportInteractionWaitOwner) -> None:
    state = obj.state
    state.viewport_interaction_release_pending = False
    state.viewport_interaction_completion_reason = None
    state.viewport_interaction_mode = False
    state.viewport_interaction_stars = None
    obj._sync_viewport_interaction_chrome_visibility()
    obj.request_client_update()
