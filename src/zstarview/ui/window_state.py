from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from PySide6.QtCore import QPoint
from PySide6.QtGui import QImage

from ..types import CelestialData, StarsTable, UrbanOutlinePolyline, ViewCenterAltAz


@dataclass
class SkyWindowState:
    """Mutable UI/runtime state kept separate from the window shell."""

    render_view_center: ViewCenterAltAz
    rotation_step: float = 5.0
    interaction_idle_ms: int = 300
    interaction_mode: bool = False
    viewport_interaction_idle_ms: int = 700
    viewport_interaction_mode: bool = False
    viewport_interaction_stars: Optional[StarsTable] = None
    mouse_pos: Optional[QPoint] = None
    jump_highlight_name: Optional[str] = None
    jump_highlight_altaz: Optional[ViewCenterAltAz] = None
    jump_highlight_until_ms: float = 0.0
    sky_update_pending: bool = False
    pending_star_vmag_limit: Optional[float] = None
    cloud_repaint_deferred: bool = False
    last_star_render_stats: Optional[tuple[int, int, int, int]] = None
    celestial_data: Optional[CelestialData] = None
    sky_disc_base_size: int = 1024
    sky_disc_image: Optional[QImage] = None
    cloud_base_size: int = 256
    terrain_horizon_profile: Optional[list[tuple[float, float]]] = None
    urban_outlines: Optional[list[UrbanOutlinePolyline]] = None
