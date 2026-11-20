"""Method contract shared by the SkyWindow mixins.

The SkyWindow mixins call methods that another mixin or
``SkyWindowCoreMixin`` defines. This class declares those signatures so each
mixin can be type-checked on its own.

The stubs exist only for the type checker (``TYPE_CHECKING``). They must not
exist at runtime: ``SkyWindowClientWidget`` composes the render mixins and
forwards missing attributes to the owner window through ``__getattr__``, and a
runtime stub would shadow that forwarding.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    import astropy.time

    from ..types import ViewerData


class SkyWindowSharedMethods:
    """Type-only method declarations for the SkyWindow mixins."""

    if TYPE_CHECKING:
        # Defined in window_update_cloud.py.
        def _geo_satellite_mode_active(self) -> bool: ...

        # Defined in window.py.
        def _predicted_cloud_satellite(self) -> str: ...
        def _aircraft_layer_enabled(self) -> bool: ...
        def _current_time_obj(self) -> astropy.time.Time: ...
        def _queue_periodic_debug_snapshot(self, payload: dict) -> None: ...
        @staticmethod
        def _resolve_periodic_debug_snapshot_dir() -> Path | None: ...
        def _satellite_layer_enabled(self) -> bool: ...
        def _start_persistent_search_refresh(
            self, *, reason: str = "timer"
        ) -> bool: ...
        def request_client_update(self) -> None: ...

        # Defined in window_render.py.
        def _viewer_data_for_render(self) -> ViewerData: ...
