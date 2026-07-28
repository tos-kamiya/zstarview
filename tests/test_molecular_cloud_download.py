from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from zstarview.molecular_cloud_download import (
    _normalize_display,
    _repair_short_zero_runs,
    prepare_akari_data,
)


def _fits_payload(value: float) -> bytes:
    data = np.full((8, 12), value, dtype=np.float32)
    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = data.shape[1]
    header["NAXIS2"] = data.shape[0]
    header["CTYPE1"] = "GLON-MOL"
    header["CTYPE2"] = "GLAT-MOL"
    header["CRPIX1"] = 6.5
    header["CRPIX2"] = 4.5
    header["CRVAL1"] = 0.0
    header["CRVAL2"] = 0.0
    header["CDELT1"] = -30.0
    header["CDELT2"] = 30.0
    stream = io.BytesIO()
    fits.PrimaryHDU(data=data, header=header).writeto(stream)
    return stream.getvalue()


def test_normalize_display_preserves_empty_pixels() -> None:
    output, parameters = _normalize_display(np.array([[0.0, 1.0, 2.0, 3.0]], dtype=np.float32))
    assert output[0, 0] == 0
    assert output[0, -1] > output[0, 1]
    assert parameters["high"] > parameters["low"]


def test_repair_short_zero_runs_interpolates_and_wraps_longitude() -> None:
    data = np.array([[10.0, 0.0, 0.0, 40.0, 50.0, 0.0, 60.0, 0.0]], dtype=np.float32)

    output = _repair_short_zero_runs(data, max_width=2)

    np.testing.assert_allclose(output, [[10.0, 20.0, 30.0, 40.0, 50.0, 55.0, 60.0, 35.0]])


def test_repair_short_zero_runs_includes_values_below_threshold() -> None:
    data = np.array([[10.0, 1.0, 1.0, 40.0]], dtype=np.float32)

    output = _repair_short_zero_runs(data, max_width=2, value_threshold=5.0)

    np.testing.assert_allclose(output, [[10.0, 20.0, 30.0, 40.0]])


def test_repair_short_zero_runs_leaves_wide_and_unbounded_runs() -> None:
    data = np.array([[10.0, 0.0, 0.0, 0.0, 0.0, 50.0, 0.0, 0.0]], dtype=np.float32)

    output = _repair_short_zero_runs(data, max_width=1)

    np.testing.assert_array_equal(output, data)


def test_prepare_akari_data_downloads_bands_and_writes_manifest(tmp_path: Path) -> None:
    payloads = {"WideS": _fits_payload(2.0), "WideL": _fits_payload(4.0), "160": _fits_payload(8.0)}
    existing_root = tmp_path / "akari-far-infrared-all-sky" / "release-1" / "schema-1"
    existing_root.mkdir(parents=True)
    sentinel = existing_root / "kept-for-inspection.txt"
    sentinel.write_text("keep", encoding="ascii")

    class Response:
        def __init__(self, payload: bytes) -> None:
            self._stream = io.BytesIO(payload)
            self.headers = {"Content-Length": str(len(payload))}

        def __enter__(self) -> "Response":
            return self

        def __exit__(self, *_args: object) -> None:
            self._stream.close()

        def read(self, size: int = -1) -> bytes:
            return self._stream.read(size)

    def fake_urlopen(request: object, timeout: float) -> Response:
        del timeout
        url = str(getattr(request, "full_url"))
        key = next(key for key in payloads if key in url)
        return Response(payloads[key])

    root = prepare_akari_data(
        bands=("90", "140", "160"),
        cache_dir=tmp_path,
        width=16,
        height=8,
        urlopen=fake_urlopen,
    )
    output = np.load(root / "akari-galactic-display.npz")
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    assert output["data"].shape == (3, 8, 16)
    assert manifest["bands_um"] == [90, 140, 160]
    assert manifest["coordinate_system"] == "galactic_plate_carree"
    assert sorted(path.name for path in root.glob("*.fits")) == [
        "akari_mollweide_160_1_4096.fits",
        "akari_mollweide_WideL_1_4096.fits",
        "akari_mollweide_WideS_1_4096.fits",
    ]
    assert sentinel.read_text(encoding="ascii") == "keep"

    prepare_akari_data(
        bands=("90", "140", "160"),
        cache_dir=tmp_path,
        width=16,
        height=8,
        delete_source=True,
        urlopen=fake_urlopen,
    )
    assert not list(root.glob("*.fits"))
