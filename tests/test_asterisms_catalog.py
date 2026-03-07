from pathlib import Path

from zstarview.catalog import load_star_catalog
from zstarview.asterisms import ASTERISMS, ASTERISM_KEYS_BY_SOURCE_ID, pick_rotating_asterism


def _catalog_source_ids() -> set[str]:
    csv_path = Path(__file__).resolve().parents[1] / "src" / "zstarview" / "data" / "stars" / "stars_base.csv"
    df = load_star_catalog(str(csv_path), vmag_threshold=None)
    return {str(source_id).strip() for source_id in df.get_column("SourceId").to_list() if str(source_id).strip()}


def test_asterism_definitions_present() -> None:
    assert ASTERISMS
    assert len({a.key for a in ASTERISMS}) == len(ASTERISMS)


def test_asterism_source_ids_exist_in_catalog() -> None:
    source_ids = _catalog_source_ids()
    missing: list[str] = []
    for asterism in ASTERISMS:
        for source_id in asterism.source_ids:
            if source_id not in source_ids:
                missing.append(f"{asterism.key}:{source_id}")
    assert not missing


def test_asterism_key_lookup_by_source_id() -> None:
    assert "winter_triangle" in ASTERISM_KEYS_BY_SOURCE_ID["HIP32349"]
    assert "summer_triangle" in ASTERISM_KEYS_BY_SOURCE_ID["HIP91262"]
    assert "southern_cross" in ASTERISM_KEYS_BY_SOURCE_ID["HIP60718"]
    assert "southern_pointers" in ASTERISM_KEYS_BY_SOURCE_ID["HIP71683"]
    assert "keystone" in ASTERISM_KEYS_BY_SOURCE_ID["HIP84380"]
    assert "cassiopeia_w" in ASTERISM_KEYS_BY_SOURCE_ID["HIP3179"]
    assert "jobs_coffin" in ASTERISM_KEYS_BY_SOURCE_ID["HIP101769"]


def test_pick_rotating_asterism_uses_slot_modulo() -> None:
    source_id, keys = next((source_id, keys) for source_id, keys in ASTERISM_KEYS_BY_SOURCE_ID.items() if len(keys) >= 2)
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
