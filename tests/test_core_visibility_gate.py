from __future__ import annotations

import pytest

from zstarview.clouddisc import CloudDisc, CloudDiscConfig, VisibilityError


def test_render_now_raises_visibility_error_for_unsupported_region(tmp_path) -> None:
    clouddisc = CloudDisc(CloudDiscConfig(cache_dir=tmp_path))
    # Around lon=20 at the equator, none of G16/G18/HIMAWARI are visible enough.
    with pytest.raises(VisibilityError, match="No supported satellite for this region"):
        clouddisc.render_now(
            lat=0.0,
            lon=20.0,
            alt=90.0,
            az=180.0,
            radius_px=128,
        )
