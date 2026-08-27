"""Cloud source fetch task boundary for future out-of-process workers."""

from __future__ import annotations

import datetime as dt
import logging
import threading
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol

from ..diagnostics import DiagnosticSink, emit_diagnostic
from ..providers.goes import GoesProvider
from ..providers.hima import HimaProvider
from ..providers.select import GOES_SATELLITES, visible_satellites
from ..types import (
    CloudBandData,
    CloudSourceData,
    DataNotFoundError,
    DownloadCancelledError,
    SourceKey,
    VisibilityError,
    round_down_utc_to_slot,
)
from .constants import DEFAULT_CLOUD_SHELLS_KM

logger = logging.getLogger(__name__)


class CloudSourceFetchContext(Protocol):
    """Minimal context required to run a cloud source fetch task."""

    goes: GoesProvider
    hima: HimaProvider

    def make_source_key(
        self,
        *,
        lat: float,
        lon: float,
        when_utc: dt.datetime | None = None,
    ) -> SourceKey: ...


@dataclass(frozen=True, slots=True)
class CloudSourceFetchRequest:
    """Normalized request payload for a cloud source fetch worker."""

    lat: float
    lon: float
    when_utc: dt.datetime | None = None
    cloud_shells_km: tuple[float, ...] = DEFAULT_CLOUD_SHELLS_KM

    def __post_init__(self) -> None:
        object.__setattr__(self, "lat", float(self.lat))
        object.__setattr__(self, "lon", float(self.lon))
        if self.when_utc is not None:
            object.__setattr__(self, "when_utc", round_down_utc_to_slot(self.when_utc))
        object.__setattr__(self, "cloud_shells_km", tuple(float(v) for v in self.cloud_shells_km))


def build_cloud_source_fetch_request(
    *,
    lat: float,
    lon: float,
    when_utc: dt.datetime | None = None,
    cloud_shells_km: Sequence[float] = DEFAULT_CLOUD_SHELLS_KM,
) -> CloudSourceFetchRequest:
    """Build a normalized cloud source fetch request."""
    return CloudSourceFetchRequest(
        lat=lat,
        lon=lon,
        when_utc=when_utc,
        cloud_shells_km=tuple(float(v) for v in cloud_shells_km),
    )


def fetch_cloud_source(
    context: CloudSourceFetchContext,
    request: CloudSourceFetchRequest,
    *,
    abort_event: threading.Event | None = None,
    diagnostic_sink: DiagnosticSink | None = None,
) -> CloudSourceData:
    """Fetch cloud source data for the given observer request."""
    logger.info("Cloud source lookup started...")
    source_key = context.make_source_key(lat=request.lat, lon=request.lon, when_utc=request.when_utc)
    sat = source_key.satellite
    when = source_key.timeslot_utc
    sat_used = sat
    shell_max_km = max(float(v) for v in request.cloud_shells_km) if request.cloud_shells_km else (6371.0 + 5.0)
    logger.debug(
        "Cloud source request resolved: sat=%s provider=%s timeslot=%s",
        sat,
        source_key.provider,
        when.isoformat(),
    )
    emit_diagnostic(
        diagnostic_sink,
        "resolve_source",
        "ok",
        "Cloud source request resolved",
        satellite=sat,
        provider=source_key.provider,
        timeslot_utc=when,
        lat=request.lat,
        lon=request.lon,
    )
    if sat in GOES_SATELLITES:
        goes_visible = tuple(visible_satellites(request.lat, request.lon, GOES_SATELLITES))
        logger.debug(
            "GOES request context: visible=%s",
            ",".join(goes_visible) if goes_visible else "(none)",
        )
        res, sat_used = context.goes.fetch_bt_c13_with_failover(
            sat=sat,
            when_utc=when,
            allowed_sats=goes_visible,
            abort_event=abort_event,
            diagnostic_sink=diagnostic_sink,
        )
        da, used_time, src_paths = res
        product = "CMIPF-C13"
    elif sat == "HIMAWARI":
        da, used_time, src_paths = context.hima.fetch_bt_c13(
            when_utc=when,
            observer_lat=request.lat,
            observer_lon=request.lon,
            cloud_shell_km=shell_max_km,
            abort_event=abort_event,
        )
        product = "ISatSS-B13"
    else:
        raise VisibilityError(f"No suitable satellite provider found for '{sat}'")
    auxiliary_bands: dict[str, CloudBandData] = {}
    try:
        if sat_used in GOES_SATELLITES:
            fetch_b16 = getattr(context.goes, "fetch_bt_c16_for_c13", None)
            if fetch_b16 is None:
                raise DataNotFoundError("GOES B16 companion provider is unavailable")
            b16_da, b16_time, b16_paths = fetch_b16(
                sat_used,
                used_time,
                src_paths[0],
                abort_event=abort_event,
                diagnostic_sink=diagnostic_sink,
            )
            b16_product = "CMIPF-C16"
        else:
            fetch_b16 = getattr(context.hima, "fetch_bt_b16_for_b13", None)
            if fetch_b16 is None:
                raise DataNotFoundError("Himawari B16 companion provider is unavailable")
            b16_da, b16_time, b16_paths = fetch_b16(
                used_time,
                src_paths,
                abort_event=abort_event,
            )
            b16_product = "ISatSS-B16"
        if b16_time != used_time:
            raise DataNotFoundError("Cloud companion band time mismatch")
        auxiliary_bands["B16"] = CloudBandData(
            band=16,
            data_array=b16_da,
            product=b16_product,
            time_utc=b16_time,
            src_paths=b16_paths,
        )
        logger.info("Using collocated %s companion data", b16_product)
    except DownloadCancelledError:
        raise
    except Exception as exc:
        logger.info("B16 companion unavailable; using B13-only cloud path: %s", exc)

    logger.info("Using %s (%s) data from time=%s", sat_used, product, used_time.isoformat())
    logger.info("Cloud source lookup ready.")
    emit_diagnostic(
        diagnostic_sink,
        "resolve_source",
        "ok",
        "Cloud source lookup ready",
        satellite=sat_used,
        product=product,
        time_utc=used_time,
        src_paths=src_paths,
    )
    source_expected_count = getattr(da, "attrs", {}).get("source_expected_count")
    source_available_count = getattr(da, "attrs", {}).get("source_available_count")
    source_completeness_ratio = getattr(da, "attrs", {}).get("source_completeness_ratio")
    return CloudSourceData(
        source_key=SourceKey(
            satellite=sat_used,
            provider=("GOES" if sat_used in GOES_SATELLITES else "HIMAWARI"),
            timeslot_utc=source_key.timeslot_utc,
            sat_priority=source_key.sat_priority,
        ),
        data_array=da,
        satellite=sat_used,
        product=product,
        time_utc=used_time,
        src_paths=src_paths,
        source_expected_count=int(source_expected_count) if source_expected_count is not None else None,
        source_available_count=int(source_available_count) if source_available_count is not None else None,
        source_completeness_ratio=(
            float(source_completeness_ratio) if source_completeness_ratio is not None else None
        ),
        auxiliary_bands=auxiliary_bands,
    )
