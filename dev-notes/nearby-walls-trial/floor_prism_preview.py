"""Preview: extrude each floor down to the next lower floor level (global), via the app's wall renderer."""
import os, sys
from dataclasses import replace
from pathlib import Path
from statistics import mean
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtGui import QColor, QFont, QGuiApplication, QImage, QPainter
from zstarview.data.derived_tile_cache import parse_derived_tile_buildings
from zstarview.data.urban_outline_common import BuildingFootprint
from zstarview.render.terrain import draw_urban_outlines
from zstarview.types import ScreenGeometry, ViewerData
from zstarview.urban_outline_layer import _nearby_wall_faces, _project_observer_centric_urban_outline_result

HOME = Path.home()


def tokyo():
    tile = HOME / ".cache/zstarview/plateau_buildings/13103/bldg/53393599_bldg_6697_op.json"
    ground = 18.4
    tower = next(b for b in parse_derived_tile_buildings(tile) if "7aff4a51-be8b-405b-abe4-ac489697cbc8" in b.building_id)
    tower = replace(tower, ground_elevation_m=ground)
    rings = []
    for surface in tower.roof_surfaces_lonlat:
        pts = surface[:-1] if surface[0] == surface[-1] else surface
        rings.append((round(mean(p[2] for p in pts) - ground), tuple((p[0], p[1]) for p in surface)))
    levels = sorted({lv for lv, _ in rings})
    parts = []
    for lv, ring in rings:
        i = levels.index(lv)
        parts.append(BuildingFootprint("tokyo-floor", float(lv), (ring,), ground_elevation_m=ground,
                                       min_height_m=float(levels[i - 1] if i else 0)))
    print("Tokyo Tower levels:", levels, "prisms:", len(parts), flush=True)
    return (tower,), tuple(parts), (35.65691966858282, 139.74599240924002), ground, (5.0, 345.0)


def london():
    tile = next(HOME.glob(".cache/zstarview/overture_buildings/overture_building_part_lat51p5009_lonm0p1218*/bldg/overture*.json"))
    parts = tuple(b for b in parse_derived_tile_buildings(tile) if b.parent_building_id == "81f82996-243f-4ab3-85c2-d2c72d863dd6")
    levels = sorted({round(b.height_m, 1) for b in parts})
    trial = []
    for b in parts:
        i = levels.index(round(b.height_m, 1))
        trial.append(replace(b, min_height_m=levels[i - 1] if i else 0.0))
    print("London Eye parts:", len(parts), "levels:", len(levels), "top:", levels[-1], flush=True)
    return parts, tuple(trial), (51.5009, -0.1218), 1.0, (12.0, 31.0)


def render(name, data, output, radius_m=500.0):
    current, trial, (lat, lon), ground, center = data
    viewer = ViewerData(location=(lat, lon), timezone_name="UTC", city_name=name, observer_height_m=1.7,
                        ground_elevation_m=ground, view_center=center)
    image = QImage(1400, 700, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(QColor("#101821"))
    painter = QPainter(image)
    for i, (label, buildings) in enumerate((("Current walls", current), ("Trial: floor extruded to next lower floor", trial))):
        walls = _nearby_wall_faces(buildings, lat_deg=lat, lon_deg=lon, observer_elevation_m=ground + 1.7, max_distance_m=radius_m)
        print(f"  {name} / {label}: wall patches {len(walls)}", flush=True)
        projected = _project_observer_centric_urban_outline_result(SimpleNamespace(outlines=(), walls=walls),
                                                                  view_center=center, edge_fov_deg=viewer.edge_fov_deg)
        painter.save(); painter.translate(700 * i, 0)
        draw_urban_outlines(painter, ScreenGeometry(center=(350, 350), radius=350), viewer, projected or [], opacity=0.7)
        painter.setPen(QColor("white")); painter.setFont(QFont("Sans Serif", 12)); painter.drawText(12, 25, f"{name}: {label}")
        painter.restore()
    painter.end()
    if not image.save(str(output)):
        raise RuntimeError("could not save image")


if __name__ == "__main__":
    out = Path(sys.argv[1])
    app = QGuiApplication([])
    render("Tokyo Tower", tokyo(), out / "tokyo-tower-floor-prisms-preview.png")
    render("London Eye", london(), out / "london-eye-floor-prisms-preview.png")
