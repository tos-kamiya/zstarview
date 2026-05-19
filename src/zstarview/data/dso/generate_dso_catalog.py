#!/usr/bin/env python3
"""
Build DSO (deep-sky object) catalog CSV using pyongc/OpenNGC data.

Output format (UTF-8 CSV):
Id,Name,Type,RAh,Dec,Vmag,MajorArcmin,MinorArcmin,PAdeg,SourceCatalog
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

CSV_COLUMNS = [
    "Id",
    "Name",
    "Type",
    "RAh",
    "Dec",
    "Vmag",
    "MajorArcmin",
    "MinorArcmin",
    "PAdeg",
    "SourceCatalog",
]


def _safe_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        f = float(value)
        return f if math.isfinite(f) else None
    s = str(value).strip()
    if not s:
        return None
    try:
        f = float(s)
    except ValueError:
        return None
    return f if math.isfinite(f) else None


def _split_tokens(text: str) -> List[str]:
    return [t for t in re.split(r"[\s:]+", text.strip()) if t]


def _parse_ra_hours(value: Any) -> Optional[float]:
    f = _safe_float(value)
    if f is not None:
        return f / 15.0 if f > 24.0 else f
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    toks = _split_tokens(s)
    if len(toks) < 2:
        return None
    try:
        h = float(toks[0])
        m = float(toks[1])
        sec = float(toks[2]) if len(toks) >= 3 else 0.0
    except ValueError:
        return None
    return h + (m / 60.0) + (sec / 3600.0)


def _parse_dec_deg(value: Any) -> Optional[float]:
    f = _safe_float(value)
    if f is not None:
        return f
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    sign = -1.0 if s.startswith("-") else 1.0
    toks = _split_tokens(s.lstrip("+-"))
    if len(toks) < 2:
        return None
    try:
        d = float(toks[0])
        m = float(toks[1])
        sec = float(toks[2]) if len(toks) >= 3 else 0.0
    except ValueError:
        return None
    return sign * (abs(d) + (m / 60.0) + (sec / 3600.0))


def _iter_attr_candidates(obj: Any, names: Sequence[str]) -> Iterator[Any]:
    if isinstance(obj, dict):
        lower_to_key = {str(k).lower(): k for k in obj.keys()}
        for name in names:
            key = lower_to_key.get(name.lower())
            if key is not None:
                yield obj.get(key)
        return
    for name in names:
        if hasattr(obj, name):
            yield getattr(obj, name)


def _pick_value(obj: Any, names: Sequence[str]) -> Any:
    for value in _iter_attr_candidates(obj, names):
        if value is None:
            continue
        if isinstance(value, str) and not value.strip():
            continue
        return value
    return None


def _normalize_obj_type(raw_type: Optional[str]) -> str:
    if raw_type is None:
        return ""
    s = str(raw_type).strip().lower().replace("_", " ").replace("-", " ")
    compact = s.replace(" ", "")
    mapping = {
        "g": "galaxy",
        "galaxy": "galaxy",
        "oc": "open_cluster",
        "opencluster": "open_cluster",
        "gc": "globular_cluster",
        "globularcluster": "globular_cluster",
        "en": "emission_nebula",
        "emissionnebula": "emission_nebula",
        "rn": "reflection_nebula",
        "reflectionnebula": "reflection_nebula",
        "dn": "dark_nebula",
        "darknebula": "dark_nebula",
        "pn": "planetary_nebula",
        "planetarynebula": "planetary_nebula",
        "snr": "supernova_remnant",
        "supernovaremnant": "supernova_remnant",
        "asterism": "asterism",
        "starcloud": "star_cloud",
        "partofgalaxy": "part_of_galaxy",
        "groupofgalaxies": "group_of_galaxies",
        "galaxypair": "galaxy_pair",
        "galaxytriplet": "galaxy_triplet",
        "hiiionizedregion": "hii_ionized_region",
        "starcluster+nebula": "star_cluster_nebula",
        "objectofother/unknowntype": "unknown_object",
        "associationofstars": "association_of_stars",
        "doublestar": "double_star",
        "nonexistentobject": "nonexistent_object",
        "duplicatedrecord": "duplicated_record",
        "novastar": "nova_star",
    }
    return mapping.get(compact, compact or s)


def _normalize_id(raw_id: str) -> str:
    t = raw_id.strip().upper().replace(" ", "")
    if t.startswith("MESSIER"):
        return "M" + t[7:]
    if t.startswith("NGC"):
        return "NGC" + t[3:]
    if t.startswith("IC"):
        return "IC" + t[2:]
    return t


def _infer_source_catalog(obj_id: str) -> str:
    s = obj_id.upper()
    if s.startswith("M"):
        return "messier"
    if s.startswith("NGC"):
        return "ngc"
    if s.startswith("IC"):
        return "ic"
    return "other"


def _belongs_to_allowed_catalogs(obj_id: str, allowed: set[str]) -> bool:
    src = _infer_source_catalog(obj_id)
    return src in allowed


def _fmt(v: Optional[float], digits: int = 6) -> str:
    if v is None:
        return ""
    return f"{v:.{digits}f}"


def _extract_row(obj: Any) -> Optional[Dict[str, str]]:
    # pyongc.ongc.Dso fast-path
    if hasattr(obj, "name") and hasattr(obj, "identifiers") and hasattr(obj, "magnitudes"):
        canonical = _normalize_id(str(getattr(obj, "name", "") or ""))
        if not canonical:
            return None
        identifiers = getattr(obj, "identifiers", None)
        messier_id: Optional[str] = None
        common_name = ""
        if isinstance(identifiers, tuple):
            if len(identifiers) >= 1 and identifiers[0]:
                messier_id = _normalize_id(str(identifiers[0]))
            if len(identifiers) >= 4 and isinstance(identifiers[3], list) and identifiers[3]:
                common_name = str(identifiers[3][0]).strip()
        obj_id = messier_id or canonical

        ra_h = _parse_ra_hours(getattr(obj, "ra", None))
        dec_deg = _parse_dec_deg(getattr(obj, "dec", None))
        if ra_h is None or dec_deg is None:
            return None

        vmag: Optional[float] = None
        magnitudes = getattr(obj, "magnitudes", None)
        if isinstance(magnitudes, tuple):
            for m in magnitudes:
                vmag = _safe_float(m)
                if vmag is not None:
                    break

        major = None
        minor = None
        pa = None
        dimensions = getattr(obj, "dimensions", None)
        if isinstance(dimensions, tuple):
            if len(dimensions) >= 1:
                major = _safe_float(dimensions[0])
            if len(dimensions) >= 2:
                minor = _safe_float(dimensions[1])
            if len(dimensions) >= 3:
                pa = _safe_float(dimensions[2])

        obj_type = _normalize_obj_type(str(getattr(obj, "type", "") or ""))
        name = common_name or canonical
        return {
            "Id": obj_id,
            "Name": name,
            "Type": obj_type,
            "RAh": _fmt(ra_h),
            "Dec": _fmt(dec_deg),
            "Vmag": _fmt(vmag, digits=3),
            "MajorArcmin": _fmt(major, digits=3),
            "MinorArcmin": _fmt(minor, digits=3),
            "PAdeg": _fmt(pa, digits=3),
            "SourceCatalog": _infer_source_catalog(obj_id),
        }

    raw_id = _pick_value(obj, ("id", "name", "identifier", "designation", "object_name", "main_id"))
    if raw_id is None:
        return None
    obj_id = _normalize_id(str(raw_id))
    if not obj_id:
        return None

    ra_h = _parse_ra_hours(_pick_value(obj, ("ra", "ra_h", "ra_hours", "rah", "ra2000")))
    dec_deg = _parse_dec_deg(_pick_value(obj, ("dec", "dec_deg", "de", "de2000")))
    if ra_h is None or dec_deg is None:
        return None

    vmag = _safe_float(_pick_value(obj, ("vmag", "mag", "magnitude", "visual_magnitude")))
    major = _safe_float(_pick_value(obj, ("major_axis", "major", "maj_ax", "majorarcmin", "size_major")))
    minor = _safe_float(_pick_value(obj, ("minor_axis", "minor", "min_ax", "minorarcmin", "size_minor")))
    pa = _safe_float(_pick_value(obj, ("position_angle", "pa", "pos_angle")))
    name = str(_pick_value(obj, ("common_name", "name", "names", "popular_name")) or "").strip()
    obj_type = _normalize_obj_type(str(_pick_value(obj, ("type", "obj_type", "category", "kind")) or ""))

    return {
        "Id": obj_id,
        "Name": name,
        "Type": obj_type,
        "RAh": _fmt(ra_h),
        "Dec": _fmt(dec_deg),
        "Vmag": _fmt(vmag, digits=3),
        "MajorArcmin": _fmt(major, digits=3),
        "MinorArcmin": _fmt(minor, digits=3),
        "PAdeg": _fmt(pa, digits=3),
        "SourceCatalog": _infer_source_catalog(obj_id),
    }


def _iter_pyongc_objects() -> Iterable[Any]:
    try:
        import pyongc  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "pyongc is required. Install dev deps first: uv pip install -e '.[dev]'"
        ) from exc

    # pyongc.ongc.listObjects() is the primary API in current pyongc.
    try:
        from pyongc import ongc  # type: ignore

        if hasattr(ongc, "listObjects"):
            objs = ongc.listObjects()
            if isinstance(objs, Iterable):
                return objs
    except ModuleNotFoundError:
        raise
    except Exception:
        pass

    # Support additional API variants.
    candidates: List[Tuple[Any, str]] = [
        (pyongc, "objects"),
        (pyongc, "all_objects"),
        (pyongc, "list_objects"),
        (pyongc, "get_all_objects"),
    ]
    ongc = getattr(pyongc, "ongc", None)
    if ongc is not None:
        candidates.extend(
            [
                (ongc, "objects"),
                (ongc, "all_objects"),
                (ongc, "list_objects"),
                (ongc, "get_all_objects"),
            ]
        )

    for holder, name in candidates:
        if not hasattr(holder, name):
            continue
        value = getattr(holder, name)
        if callable(value):
            value = value()
        if isinstance(value, Iterable):
            return value

    raise RuntimeError("Failed to locate iterable OpenNGC objects in pyongc API.")


def _dedupe_rows(rows: Iterable[Dict[str, str]]) -> List[Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    priority = {"messier": 0, "ngc": 1, "ic": 2, "other": 3}
    for row in rows:
        obj_id = row["Id"]
        current = out.get(obj_id)
        if current is None:
            out[obj_id] = row
            continue
        cur_p = priority.get(current.get("SourceCatalog", "other"), 9)
        new_p = priority.get(row.get("SourceCatalog", "other"), 9)
        if new_p < cur_p:
            out[obj_id] = row
    return sorted(out.values(), key=lambda r: r["Id"])


def generate(args: argparse.Namespace) -> None:
    allowed_catalogs = {v.strip().lower() for v in args.include.split(",") if v.strip()}
    if not allowed_catalogs:
        allowed_catalogs = {"messier", "ngc", "ic"}
    allowed_types = {v.strip().lower() for v in args.types.split(",") if v.strip()}
    if not allowed_types:
        allowed_types = {"galaxy", "open_cluster", "globular_cluster"}

    raw_objects = _iter_pyongc_objects()
    rows: List[Dict[str, str]] = []
    skipped_missing_coord = 0

    for obj in raw_objects:
        row = _extract_row(obj)
        if row is None:
            skipped_missing_coord += 1
            continue
        if not _belongs_to_allowed_catalogs(row["Id"], allowed_catalogs):
            continue
        if "all" not in allowed_types and row["Type"].lower() not in allowed_types:
            continue
        if args.max_vmag is not None:
            vmag = _safe_float(row["Vmag"])
            if vmag is None or vmag > float(args.max_vmag):
                continue
        rows.append(row)

    deduped = _dedupe_rows(rows)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(deduped)

    type_counts: Dict[str, int] = {}
    for row in deduped:
        t = row["Type"] or "unknown"
        type_counts[t] = type_counts.get(t, 0) + 1

    print(f"wrote: {args.output}")
    print(f"rows: {len(deduped)}")
    print(f"skipped(no usable id/coord): {skipped_missing_coord}")
    print("type counts:")
    for t in sorted(type_counts):
        print(f"  {t}: {type_counts[t]}")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description="Generate DSO CSV from OpenNGC/pyongc.")
    ap.add_argument(
        "--include",
        default="messier,ngc,ic",
        help="Comma separated catalogs to include (messier,ngc,ic,other).",
    )
    ap.add_argument(
        "--types",
        default="galaxy,open_cluster,globular_cluster",
        help="Comma separated object types to keep (default: galaxy,open_cluster,globular_cluster). Use 'all' to disable type filtering.",
    )
    ap.add_argument(
        "--max-vmag",
        type=float,
        default=12.0,
        help="Magnitude limit; keeps only objects with Vmag <= this value (default: 12.0).",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=here.parent / "dso.csv",
        help="Output CSV path.",
    )
    return ap.parse_args()


def main() -> None:
    generate(parse_args())


if __name__ == "__main__":
    main()
