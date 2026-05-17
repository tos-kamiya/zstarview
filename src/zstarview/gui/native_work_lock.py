from __future__ import annotations

import threading

# Keep heavy native-backed background work from overlapping across controllers.
#
# The sky worker and cloud pipeline both call into native extension modules
# (NumPy, pyproj, xarray/netCDF, HDF5-backed readers). Serializing those
# sections avoids concurrent native access patterns that can segfault on some
# platforms.
HEAVY_NATIVE_WORK_LOCK = threading.Lock()
