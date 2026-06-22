# -*- coding: utf-8 -*-
"""Constants for the observer-centric (altitude, azimuth) cloud grid.

These values govern the resolution of the cached `CloudAltAzGrid` and the
appearance of the debug / MVP circle renderer. They are grouped here so
tuning does not scatter across modules.
"""

# --- Grid geometry ---------------------------------------------------------
# Default angular resolution of the alt/az grid.  The grid covers
# altitude [0, 90] deg and azimuth [0, 360) deg.
ALT_AZ_GRID_AZ_BINS: int = 720
ALT_AZ_GRID_ALT_BINS: int = 90

ALT_AZ_GRID_AZ_MIN_DEG: float = 0.0
ALT_AZ_GRID_AZ_MAX_DEG: float = 360.0
ALT_AZ_GRID_ALT_MIN_DEG: float = 0.0
ALT_AZ_GRID_ALT_MAX_DEG: float = 90.0

# --- Ingestion sampling ----------------------------------------------------
# Step size in degrees for the geographic lat/lon sample grid used when
# transforming satellite BT data into the alt/az grid.  Smaller values give
# finer detail but increase ingestion cost.
ALT_AZ_GEO_SAMPLE_STEP_DEG: float = 0.2

# Half-width of the geographic sample box centred on the observer.  The
# visible cloud shells are within a few degrees of the observer, so a small
# box is sufficient; this value is generous to handle high viewpoints and
# slight projection margin.
ALT_AZ_GEO_SAMPLE_EXTENT_DEG: float = 5.0

# Number of grid cells by which to expand the coverage region when marking
# missing alt/az cells.  A cell is considered missing only if it lies within
# this expanded neighborhood of any valid sample.
ALT_AZ_MISSING_NEIGHBORHOOD_CELLS: int = 2

# --- MVP circle renderer ---------------------------------------------------
# Radius and opacity range for drawing cloud cells as white circles.
ALT_AZ_CIRCLE_BASE_RADIUS_PX: float = 2.0
ALT_AZ_CIRCLE_MAX_RADIUS_PX: float = 12.0
ALT_AZ_CIRCLE_OPACITY_SCALE: float = 0.6

# Cloud-amount threshold below which a grid cell is not drawn at all.
ALT_AZ_CIRCLE_AMOUNT_THRESHOLD: float = 0.03

# --- Debug export ----------------------------------------------------------
# Sub-directory under the cache root where debug images of the alt/az grid
# are written when debugging is enabled.
ALT_AZ_DEBUG_SUBDIR: str = "altaz_debug"
