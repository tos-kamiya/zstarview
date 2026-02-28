#!/usr/bin/env python3
"""
Build merged star catalogs from Hipparcos (+ optional Tycho-2 CSV).

This script is designed as the canonical catalog builder for zstarview.
It can run with Hipparcos only, then optionally merge Tycho-2 stars for
fainter ranges (e.g., 6 < Vmag <= 9).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


@dataclass(slots=True)
class StarRec:
    source_catalog: str
    source_id: str
    ra_hours: float
    dec_deg: float
    vmag: float
    bv: Optional[float]
    name: str = ""
    hip_id: Optional[int] = None
    tyc_id: Optional[str] = None

    @property
    def ra_deg(self) -> float:
        return (self.ra_hours * 15.0) % 360.0


def _fnum(text: str) -> Optional[float]:
    s = text.strip()
    if not s or s == "?":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _iter_lines(path: Path, encoding: str = "ascii") -> Iterator[str]:
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path, "r") as zf:
            names = [n for n in zf.namelist() if not n.endswith("/")]
            if not names:
                return
            # hip_main.dat.zip in this repo contains a single payload.
            with zf.open(names[0], "r") as f:
                for b in f:
                    yield b.decode(encoding, errors="ignore")
        return
    with path.open("r", encoding=encoding, errors="ignore") as f:
        yield from f


def _load_iau_name_map(path: Optional[Path]) -> Dict[int, str]:
    if path is None or not path.exists():
        return {}
    out: Dict[int, str] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            name = (row.get("proper names") or "").strip()
            hip = (row.get("HIP") or "").strip()
            if not name or not hip.isdigit():
                continue
            out[int(hip)] = name
    return out


def parse_hipparcos(
    path: Path,
    *,
    max_vmag: float,
    iau_name_by_hip: Dict[int, str],
) -> List[StarRec]:
    out: List[StarRec] = []
    for line in _iter_lines(path):
        if not line or line[0] != "H":
            continue
        hip_text = line[2:14].strip()  # cols 3-14 (1-based)
        if not hip_text.isdigit():
            continue
        hip = int(hip_text)
        vmag = _fnum(line[41:46])  # cols 42-46
        ra_deg = _fnum(line[51:63])  # cols 52-63
        dec_deg = _fnum(line[64:76])  # cols 65-76
        bmv = _fnum(line[245:251])  # cols 246-251
        if vmag is None or ra_deg is None or dec_deg is None:
            continue
        if vmag > max_vmag:
            continue
        out.append(
            StarRec(
                source_catalog="hip",
                source_id=f"HIP{hip}",
                ra_hours=ra_deg / 15.0,
                dec_deg=dec_deg,
                vmag=vmag,
                bv=bmv,
                name=iau_name_by_hip.get(hip, ""),
                hip_id=hip,
            )
        )
    return out


def _first_existing_key(row: Dict[str, str], candidates: Sequence[str]) -> Optional[str]:
    lower_to_key = {k.strip().lower(): k for k in row.keys()}
    for c in candidates:
        k = lower_to_key.get(c.lower())
        if k is not None:
            return k
    return None


def parse_tycho_csv(path: Optional[Path], *, max_vmag: float) -> List[StarRec]:
    if path is None:
        return []
    if not path.exists():
        raise FileNotFoundError(f"Tycho CSV not found: {path}")

    out: List[StarRec] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            if not row:
                continue
            id_key = _first_existing_key(row, ("TYC", "tyc", "tyc_id", "SourceId", "ID"))
            ra_key = _first_existing_key(row, ("RAh", "ra_hours", "RA_hour", "RAh2000", "RA", "ra"))
            dec_key = _first_existing_key(row, ("Dec", "dec_deg", "DEdeg", "DE", "dec"))
            vmag_key = _first_existing_key(row, ("Vmag", "VTmag", "VT", "mag", "vmag"))
            bv_key = _first_existing_key(row, ("B-V", "BV", "BT-VT", "BTminusVT", "b_v"))
            name_key = _first_existing_key(row, ("Name", "name"))

            if not ra_key or not dec_key or not vmag_key:
                continue
            ra_raw = _fnum(row.get(ra_key, ""))
            dec_deg = _fnum(row.get(dec_key, ""))
            vmag = _fnum(row.get(vmag_key, ""))
            if ra_raw is None or dec_deg is None or vmag is None:
                continue
            if vmag > max_vmag:
                continue

            # If RA value looks like degrees (> 24), convert to hours.
            ra_hours = ra_raw / 15.0 if ra_raw > 24.0 else ra_raw
            bv = _fnum(row.get(bv_key, "")) if bv_key else None
            tyc_id = (row.get(id_key, "").strip() if id_key else "")
            name = (row.get(name_key, "").strip() if name_key else "")
            out.append(
                StarRec(
                    source_catalog="tyc2",
                    source_id=tyc_id or f"TYCROW{len(out)+1}",
                    ra_hours=ra_hours,
                    dec_deg=dec_deg,
                    vmag=vmag,
                    bv=bv,
                    name=name,
                    tyc_id=tyc_id or None,
                )
            )
    return out


def parse_tycho_i259_dir(path: Optional[Path], *, max_vmag: float) -> List[StarRec]:
    """Parse Tycho-2 I/259 split files (tyc2.dat.00.gz ... tyc2.dat.19.gz)."""
    if path is None:
        return []
    if not path.exists():
        raise FileNotFoundError(f"I/259 directory not found: {path}")

    files = sorted(path.glob("tyc2.dat.*.gz"))
    if not files:
        raise FileNotFoundError(f"No tyc2.dat.*.gz files found under: {path}")

    out: List[StarRec] = []
    for fp in files:
        with gzip.open(fp, "rt", encoding="ascii", errors="ignore") as f:
            for line in f:
                if not line:
                    continue
                # Byte positions are from I/259 ReadMe (1-based).
                tyc1 = line[0:4].strip()
                tyc2 = line[5:10].strip()
                tyc3 = line[11:12].strip()
                ra_deg = _fnum(line[15:27])  # RAmdeg (J2000 mean)
                dec_deg = _fnum(line[28:40])  # DEmdeg (J2000 mean)
                bt = _fnum(line[110:116])  # BTmag
                vt = _fnum(line[123:129])  # VTmag
                hip = line[142:148].strip()
                if ra_deg is None or dec_deg is None:
                    continue
                vmag = vt if vt is not None else bt
                if vmag is None or vmag > max_vmag:
                    continue
                bv: Optional[float] = None
                if bt is not None and vt is not None:
                    # ReadMe note (7): approximate Johnson B-V.
                    bv = 0.850 * (bt - vt)
                tyc_id = "-".join(x for x in (tyc1, tyc2, tyc3) if x)
                out.append(
                    StarRec(
                        source_catalog="tyc2",
                        source_id=f"TYC{tyc_id}" if tyc_id else f"TYCROW{len(out)+1}",
                        ra_hours=ra_deg / 15.0,
                        dec_deg=dec_deg,
                        vmag=vmag,
                        bv=bv,
                        name="",
                        hip_id=int(hip) if hip.isdigit() else None,
                        tyc_id=tyc_id or None,
                    )
                )
    return out


def _ang_sep_deg(ra1_deg: float, dec1_deg: float, ra2_deg: float, dec2_deg: float) -> float:
    ra1 = math.radians(ra1_deg)
    ra2 = math.radians(ra2_deg)
    dec1 = math.radians(dec1_deg)
    dec2 = math.radians(dec2_deg)
    cos_d = math.sin(dec1) * math.sin(dec2) + math.cos(dec1) * math.cos(dec2) * math.cos(ra1 - ra2)
    cos_d = max(-1.0, min(1.0, cos_d))
    return math.degrees(math.acos(cos_d))


def _bucket_key(ra_deg: float, dec_deg: float, cell_deg: float) -> Tuple[int, int]:
    ra_bin = int((ra_deg % 360.0) // cell_deg)
    dec_bin = int((dec_deg + 90.0) // cell_deg)
    return ra_bin, dec_bin


def _neighbor_keys(key: Tuple[int, int], *, ra_bins: int) -> Iterable[Tuple[int, int]]:
    ra0, dec0 = key
    for dra in (-1, 0, 1):
        for ddec in (-1, 0, 1):
            yield ((ra0 + dra) % ra_bins, dec0 + ddec)


def _is_duplicate(
    candidate: StarRec,
    selected: List[StarRec],
    bucket_to_idx: Dict[Tuple[int, int], List[int]],
    *,
    sep_deg: float,
    max_mag_diff: float,
    cell_deg: float,
) -> bool:
    ra_bins = max(1, int(math.ceil(360.0 / cell_deg)))
    k = _bucket_key(candidate.ra_deg, candidate.dec_deg, cell_deg)
    for nk in _neighbor_keys(k, ra_bins=ra_bins):
        for idx in bucket_to_idx.get(nk, []):
            s = selected[idx]
            if abs(s.vmag - candidate.vmag) > max_mag_diff:
                continue
            if _ang_sep_deg(s.ra_deg, s.dec_deg, candidate.ra_deg, candidate.dec_deg) <= sep_deg:
                return True
    return False


def merge_catalogs(
    hip_stars: Sequence[StarRec],
    tycho_stars: Sequence[StarRec],
    *,
    hip_priority_vmag: float,
    dup_sep_arcsec: float,
    dup_mag_diff: float,
) -> List[StarRec]:
    sep_deg = dup_sep_arcsec / 3600.0
    cell_deg = max(0.02, sep_deg * 2.0)
    selected: List[StarRec] = []
    bucket_to_idx: Dict[Tuple[int, int], List[int]] = {}
    ra_bins = max(1, int(math.ceil(360.0 / cell_deg)))

    def add_star(rec: StarRec) -> None:
        idx = len(selected)
        selected.append(rec)
        k = _bucket_key(rec.ra_deg, rec.dec_deg, cell_deg)
        bucket_to_idx.setdefault(k, []).append(idx)

    # Hipparcos first: source of truth for bright stars.
    for rec in hip_stars:
        add_star(rec)

    # Tycho supplement. Skip bright side when Hipparcos should own it.
    for rec in tycho_stars:
        if rec.vmag <= hip_priority_vmag:
            continue
        if _is_duplicate(
            rec,
            selected,
            bucket_to_idx,
            sep_deg=sep_deg,
            max_mag_diff=dup_mag_diff,
            cell_deg=cell_deg,
        ):
            continue
        add_star(rec)

    # Stable output: brighter first, then RA.
    selected.sort(key=lambda r: (r.vmag, r.ra_hours))

    # Rebuild ids in case order changed.
    bucket_to_idx.clear()
    for idx, rec in enumerate(selected):
        k = _bucket_key(rec.ra_deg, rec.dec_deg, cell_deg)
        bucket_to_idx.setdefault(k, []).append(idx)
    _ = ra_bins
    return selected


def _fmt_float(v: Optional[float], digits: int) -> str:
    if v is None:
        return ""
    return f"{v:.{digits}f}"


def write_catalog(path: Path, rows: Sequence[StarRec]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "HIP",
                "HR",
                "Name",
                "RAh",
                "Dec",
                "Vmag",
                "B-V",
                "SourceCatalog",
                "SourceId",
                "TYC",
            ]
        )
        for r in rows:
            w.writerow(
                [
                    str(r.hip_id) if r.hip_id is not None else "",
                    "",
                    r.name,
                    _fmt_float(r.ra_hours, 6),
                    _fmt_float(r.dec_deg, 6),
                    _fmt_float(r.vmag, 3),
                    _fmt_float(r.bv, 3),
                    r.source_catalog,
                    r.source_id,
                    r.tyc_id or "",
                ]
            )


def _split_ranges(rows: Sequence[StarRec]) -> Dict[str, List[StarRec]]:
    out: Dict[str, List[StarRec]] = {
        "base": [],
        "extra7": [],
        "extra8": [],
        "extra9": [],
    }
    for r in rows:
        if r.vmag <= 6.0:
            out["base"].append(r)
        elif r.vmag <= 7.0:
            out["extra7"].append(r)
        elif r.vmag <= 8.0:
            out["extra8"].append(r)
        elif r.vmag <= 9.0:
            out["extra9"].append(r)
    return out


def _print_stats(label: str, rows: Sequence[StarRec]) -> None:
    n = len(rows)
    if n == 0:
        print(f"{label}: 0 rows")
        return
    bv_missing = sum(1 for r in rows if r.bv is None)
    hip_n = sum(1 for r in rows if r.source_catalog == "hip")
    tyc_n = sum(1 for r in rows if r.source_catalog == "tyc2")
    print(
        f"{label}: rows={n}, hip={hip_n}, tyc2={tyc_n}, "
        f"bv_missing={bv_missing} ({(bv_missing / n) * 100:.2f}%)"
    )


def build_catalog(args: argparse.Namespace) -> None:
    iau_map = _load_iau_name_map(args.iau_csv)
    hip = parse_hipparcos(args.hip_main, max_vmag=args.max_vmag, iau_name_by_hip=iau_map)
    tyc: List[StarRec] = []
    if args.tycho_csv is not None:
        tyc.extend(parse_tycho_csv(args.tycho_csv, max_vmag=args.max_vmag))
    if args.tycho_i259_dir is not None:
        tyc.extend(parse_tycho_i259_dir(args.tycho_i259_dir, max_vmag=args.max_vmag))
    merged = merge_catalogs(
        hip,
        tyc,
        hip_priority_vmag=args.hip_priority_vmag,
        dup_sep_arcsec=args.dup_sep_arcsec,
        dup_mag_diff=args.dup_mag_diff,
    )

    write_catalog(args.output, merged)
    split = _split_ranges(merged)
    write_catalog(args.output_dir / "stars_base.csv", split["base"])
    write_catalog(args.output_dir / "stars_extra7.csv", split["extra7"])
    write_catalog(args.output_dir / "stars_extra8.csv", split["extra8"])
    write_catalog(args.output_dir / "stars_extra9.csv", split["extra9"])

    _print_stats("hip_input", hip)
    _print_stats("tycho_input", tyc)
    _print_stats("merged_all", merged)
    _print_stats("merged_base", split["base"])
    _print_stats("merged_extra7", split["extra7"])
    _print_stats("merged_extra8", split["extra8"])
    _print_stats("merged_extra9", split["extra9"])
    print(f"wrote: {args.output}")
    print(f"wrote split files under: {args.output_dir}")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description="Generate merged star catalogs (Hipparcos + optional Tycho-2 CSV).")
    ap.add_argument(
        "--hip-main",
        type=Path,
        default=here / "hip_main.dat.zip",
        help="Path to hip_main.dat or hip_main.dat.zip",
    )
    ap.add_argument(
        "--tycho-csv",
        type=Path,
        default=None,
        help="Optional Tycho-2 normalized CSV (must include RA/Dec/Vmag-like columns).",
    )
    default_i259 = here.parent / "I-259"
    ap.add_argument(
        "--tycho-i259-dir",
        type=Path,
        default=default_i259 if default_i259.exists() else None,
        help="Optional Tycho-2 I/259 directory (tyc2.dat.*.gz).",
    )
    ap.add_argument(
        "--iau-csv",
        type=Path,
        default=here / "IAU-Catalog of Star Names (always up to date).csv",
        help="Optional IAU proper-name CSV for HIP name fill.",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=here.parent / "stars.csv",
        help="Output merged CSV path (legacy-compatible columns).",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=here,
        help="Directory to write split catalogs (base/extra7/8/9).",
    )
    ap.add_argument("--max-vmag", type=float, default=9.0)
    ap.add_argument("--hip-priority-vmag", type=float, default=6.0)
    ap.add_argument("--dup-sep-arcsec", type=float, default=5.0)
    ap.add_argument("--dup-mag-diff", type=float, default=0.75)
    return ap.parse_args()


def main() -> None:
    build_catalog(parse_args())


if __name__ == "__main__":
    main()
