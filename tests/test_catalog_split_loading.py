from __future__ import annotations

from pathlib import Path

from zstarview.catalog import load_star_catalog


HEADER = "HIP,HR,Name,RAh,Dec,Vmag,B-V,SourceCatalog,SourceId,TYC\n"


def _write_rows(path: Path, rows: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER + "".join(rows), encoding="utf-8")


def test_load_star_catalog_uses_split_files_when_available(tmp_path: Path) -> None:
    stars_csv = tmp_path / "stars.csv"
    split_dir = tmp_path / "stars"

    # Unified file has a sentinel star not present in split files.
    _write_rows(stars_csv, ["1,,,0,0,6.9,0.1,hip,HIP1,\n"])
    _write_rows(split_dir / "stars_base.csv", ["2,,,0,0,5.5,0.1,hip,HIP2,\n"])
    _write_rows(split_dir / "stars_extra7.csv", ["3,,,0,0,6.8,0.1,tyc2,TYC3,1-1-1\n"])
    _write_rows(split_dir / "stars_extra8.csv", [])
    _write_rows(split_dir / "stars_extra9.csv", [])

    df = load_star_catalog(str(stars_csv), vmag_threshold=7.0)
    ids = set(df.get_column("SourceId").to_list())
    assert ids == {"HIP2", "TYC3"}


def test_load_star_catalog_falls_back_to_unified_when_split_missing(tmp_path: Path) -> None:
    stars_csv = tmp_path / "stars.csv"
    split_dir = tmp_path / "stars"

    _write_rows(stars_csv, ["1,,,0,0,6.9,0.1,hip,HIP1,\n"])
    _write_rows(split_dir / "stars_base.csv", ["2,,,0,0,5.5,0.1,hip,HIP2,\n"])
    # Missing stars_extra7.csv should trigger unified fallback for threshold 7.0.

    df = load_star_catalog(str(stars_csv), vmag_threshold=7.0)
    ids = set(df.get_column("SourceId").to_list())
    assert ids == {"HIP1"}


def test_load_star_catalog_uses_extra10_when_available(tmp_path: Path) -> None:
    stars_csv = tmp_path / "stars.csv"
    split_dir = tmp_path / "stars"

    _write_rows(stars_csv, ["1,,,0,0,9.8,0.1,hip,HIP1,\n"])
    _write_rows(split_dir / "stars_base.csv", ["2,,,0,0,5.5,0.1,hip,HIP2,\n"])
    _write_rows(split_dir / "stars_extra7.csv", [])
    _write_rows(split_dir / "stars_extra8.csv", [])
    _write_rows(split_dir / "stars_extra9.csv", [])
    _write_rows(split_dir / "stars_extra10.csv", ["10,,,0,0,9.8,0.1,tyc2,TYC10,1-1-10\n"])

    df = load_star_catalog(str(stars_csv), vmag_threshold=10.0)
    ids = set(df.get_column("SourceId").to_list())
    assert ids == {"HIP2", "TYC10"}
