"""Static rendering dependencies are outside metadata-only history listing."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..location_resolver.resolve import ResolvedLocation


def static_source_checks(
    root: Path, args: Any, location: ResolvedLocation
) -> list[dict[str, str]]:
    """Describe enabled dependencies without reading any rendering source."""
    return [
        {
            "name": name,
            "status": "unverified",
            "reason": "Metadata-only listing; rendering source contents are not inspected",
        }
        for name, attribute in (
            ("terrain-source", "terrain_horizon_opacity"),
            ("urban-source", "urban_outline_opacity"),
            ("water-source", "water_surface_opacity"),
            ("night-lights-source", "night_light_opacity"),
            ("ridge-glow-source", "ridge_glow_opacity"),
            ("road-lights-source", "road_light_opacity"),
            ("diffuse-sky-source", "diffuse_sky_opacity"),
        )
        if getattr(args, attribute, 0) > 0
    ]
