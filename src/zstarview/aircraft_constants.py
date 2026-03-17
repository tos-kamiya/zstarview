"""Shared constants for the planned OpenSky aircraft overlay."""

AIRCRAFT_REFRESH_INTERVAL_SECONDS = 5 * 60
AIRCRAFT_BBOX_DELTA_DEG = 1.0

# OpenSky credits for `/states/all` are based on bbox area in square degrees.
AIRCRAFT_BBOX_AREA_SQ_DEG = (AIRCRAFT_BBOX_DELTA_DEG * 2.0) ** 2
