from __future__ import annotations

import pytest

from zstarview.clouddisc import CloudDisc, CloudDiscConfig, VisibilityError


def test_fetch_source_raises_visibility_error_for_unsupported_region(tmp_path) -> None:
    clouddisc = CloudDisc(CloudDiscConfig(cache_dir=tmp_path))
    # Around lon=20 at the equator, none of G16/G18/HIMAWARI are visible enough.
    with pytest.raises(VisibilityError, match="No supported satellite for this region"):
        clouddisc.fetch_source(
            lat=0.0,
            lon=20.0,
        )
