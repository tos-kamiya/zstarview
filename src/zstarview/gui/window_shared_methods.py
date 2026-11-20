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


class SkyWindowSharedMethods:
    """Type-only method declarations for the SkyWindow mixins."""

    if TYPE_CHECKING:
        # Defined in window_update_cloud.py.
        def _geo_satellite_mode_active(self) -> bool: ...

        # Defined in window.py.
        def _predicted_cloud_satellite(self) -> str: ...
