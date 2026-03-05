import csv
from pathlib import Path

from zstarview.asterisms import ASTERISMS, ASTERISM_KEYS_BY_SOURCE_ID, pick_rotating_asterism


def _catalog_source_ids() -> set[str]:
    csv_path = Path(__file__).resolve().parents[1] / "src" / "zstarview" / "data" / "stars" / "stars_base.csv"
    out: set[str] = set()
    with csv_path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            source_id = (row.get("SourceId") or "").strip()
            if source_id:
                out.add(source_id)
    return out


def test_asterism_definitions_present() -> None:
    assert ASTERISMS
    assert len({a.key for a in ASTERISMS}) == len(ASTERISMS)


def test_asterism_source_ids_exist_in_catalog() -> None:
    source_ids = _catalog_source_ids()
    missing: list[str] = []
    for asterism in ASTERISMS:
        for source_id in asterism.path:
            if source_id not in source_ids:
                missing.append(f"{asterism.key}:{source_id}")
    assert not missing


def test_asterism_key_lookup_by_source_id() -> None:
    assert "winter_triangle" in ASTERISM_KEYS_BY_SOURCE_ID["HIP32349"]
    assert "summer_triangle" in ASTERISM_KEYS_BY_SOURCE_ID["HIP91262"]


def test_pick_rotating_asterism_uses_slot_modulo() -> None:
    source_id = "HIP113963"  # included in multiple autumn asterisms
    keys = ASTERISM_KEYS_BY_SOURCE_ID[source_id]
    assert len(keys) >= 2
    first = pick_rotating_asterism(source_id, 0)
    second = pick_rotating_asterism(source_id, 1)
    wrap = pick_rotating_asterism(source_id, len(keys))
    assert first is not None
    assert second is not None
    assert wrap is not None
    assert first.key == keys[0]
    assert second.key == keys[1]
    assert wrap.key == keys[0]
