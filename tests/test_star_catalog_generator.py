from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


def _load_module():
    root = Path(__file__).resolve().parents[1]
    mod_path = root / "src" / "zstarview" / "data" / "stars" / "generate_star_catalog.py"
    spec = importlib.util.spec_from_file_location("star_catalog_gen", mod_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_merge_catalogs_prefers_hip_for_bright_range() -> None:
    mod = _load_module()
    hip = [
        mod.StarRec(
            source_catalog="hip",
            source_id="HIP1",
            ra_hours=10.0,
            dec_deg=20.0,
            vmag=5.0,
            bv=0.4,
            hip_id=1,
        )
    ]
    # Same location/magnitude-ish but Tycho. Must be dropped in bright range.
    tyc = [
        mod.StarRec(
            source_catalog="tyc2",
            source_id="TYC1",
            ra_hours=10.00001,
            dec_deg=20.00001,
            vmag=5.1,
            bv=0.5,
            tyc_id="1-1-1",
        )
    ]
    merged = mod.merge_catalogs(
        hip,
        tyc,
        hip_priority_vmag=6.0,
        dup_sep_arcsec=10.0,
        dup_mag_diff=0.75,
    )
    assert len(merged) == 1
    assert merged[0].source_catalog == "hip"


def test_merge_catalogs_adds_non_duplicate_faint_tycho() -> None:
    mod = _load_module()
    hip = [
        mod.StarRec(
            source_catalog="hip",
            source_id="HIP10",
            ra_hours=1.0,
            dec_deg=1.0,
            vmag=6.2,
            bv=0.2,
            hip_id=10,
        )
    ]
    tyc = [
        mod.StarRec(
            source_catalog="tyc2",
            source_id="TYC2",
            ra_hours=3.0,
            dec_deg=-5.0,
            vmag=8.1,
            bv=0.9,
            tyc_id="2-2-2",
        )
    ]
    merged = mod.merge_catalogs(
        hip,
        tyc,
        hip_priority_vmag=6.0,
        dup_sep_arcsec=5.0,
        dup_mag_diff=0.75,
    )
    assert len(merged) == 2
    assert {r.source_catalog for r in merged} == {"hip", "tyc2"}
