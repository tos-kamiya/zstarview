"""Render isolated LOD1/Overture wall comparisons; see the adjacent notes."""
import argparse
import os
import xml.etree.ElementTree as ET
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PySide6.QtGui import QColor, QFont, QGuiApplication, QImage, QPainter
from zstarview.data.derived_tile_cache import parse_derived_tile_buildings
from zstarview.data.urban_outline_common import BuildingFootprint
from zstarview.render.terrain import draw_urban_outlines
from zstarview.types import ScreenGeometry, ViewerData
from zstarview.urban_outline_layer import _nearby_wall_faces, _project_observer_centric_urban_outline_result


def tokyo(gml):
    root = ET.parse(gml).getroot()
    building = next(e for e in root.iter() if e.get('{http://www.opengis.net/gml}id') == 'bldg_7aff4a51-be8b-405b-abe4-ac489697cbc8')
    roof_rings = []
    for surface in building.iter():
        if not surface.tag.endswith('}RoofSurface'):
            continue
        for pos_list in surface.iter():
            if pos_list.tag.endswith('}posList') and pos_list.text:
                values = list(map(float, pos_list.text.split()))
                ring = tuple((values[i+1], values[i], values[i+2]) for i in range(0, len(values), 3))
                if len(ring) >= 4:
                    roof_rings.append(ring)
    roof_rings.sort(key=lambda ring: sum(p[2] for p in ring[:-1]) / max(1, len(ring)-1))
    levels = sorted({round(sum(p[2] for p in ring[:-1]) / max(1, len(ring)-1), 2) for ring in roof_rings})
    print('Tokyo LOD2 roof rings / levels:', len(roof_rings), len(levels), flush=True)
    trial_parts = []
    for ring in roof_rings:
        lower = min(p[2] for p in ring)
        upper = max(p[2] for p in ring)
        level = round(sum(p[2] for p in ring[:-1]) / max(1, len(ring)-1), 2)
        next_level = next((h for h in levels if h > level + 0.05), max(p[2] for r in roof_rings for p in r))
        if next_level - lower < 0.1:
            continue
        trial_parts.append(BuildingFootprint('roof-trial', next_level-lower, (tuple((p[0],p[1]) for p in ring),), ground_elevation_m=lower))
    tile = Path.home() / '.cache/zstarview/plateau_buildings/13103/bldg/53393599_bldg_6697_op.json'
    current = next(b for b in parse_derived_tile_buildings(tile) if '7aff4a51-be8b-405b-abe4-ac489697cbc8' in b.building_id)
    ground = 18.4
    return (replace(current, ground_elevation_m=ground),), tuple(trial_parts), (35.65691966858282, 139.74599240924002), ground, (5.0, 345.0), 200


def london():
    tile = next(Path.home().glob('.cache/zstarview/overture_buildings/overture_building_part_lat51p5009_lonm0p1218*/bldg/overture*.json'))
    parts = tuple(b for b in parse_derived_tile_buildings(tile) if b.parent_building_id == '81f82996-243f-4ab3-85c2-d2c72d863dd6')
    levels = sorted({b.min_height_m for b in parts})
    trial = tuple(replace(b, height_m=b.height_m) for b in parts if b.height_m > b.min_height_m)
    print('London parts / levels / trial parts:', len(parts), len(levels), len(trial), 'top:', max(b.height_m for b in parts), flush=True)
    return parts, trial, (51.5009, -0.1218), 1.0, (12.0, 31.0), 500


def render(name, data, output):
    current, trial, (lat, lon), ground, center, radius = data
    viewer = ViewerData(location=(lat, lon), timezone_name='UTC', city_name=name, observer_height_m=1.7, ground_elevation_m=ground, view_center=center)
    image = QImage(1400, 700, QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(QColor('#101821'))
    painter = QPainter(image)
    for i, (label, buildings) in enumerate((('Current walls', current), ('Trial: lower outline to next level', trial))):
        walls = _nearby_wall_faces(buildings, lat_deg=lat, lon_deg=lon, observer_elevation_m=ground+1.7, max_distance_m=radius)
        print(name, label, 'wall patches:', len(walls), flush=True)
        projected = _project_observer_centric_urban_outline_result(SimpleNamespace(outlines=(), walls=walls), view_center=center, edge_fov_deg=viewer.edge_fov_deg)
        painter.save()
        painter.translate(700*i, 0)
        draw_urban_outlines(painter, ScreenGeometry(center=(350, 350), radius=350), viewer, projected or [], opacity=0.7)
        painter.setPen(QColor('white'))
        painter.setFont(QFont('Sans Serif', 12))
        painter.drawText(12, 25, label)
        painter.restore()
    painter.end()
    if not image.save(str(output)):
        raise RuntimeError('Could not save image')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gml', type=Path)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    app = QGuiApplication([])
    render('Tokyo Tower', tokyo(args.gml), args.output_dir / 'tokyo-tower-roof-levels.png')
    render('London Eye', london(), args.output_dir / 'london-eye-roof-levels.png')
