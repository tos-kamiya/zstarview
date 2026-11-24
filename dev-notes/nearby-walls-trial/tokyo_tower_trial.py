"""Compare wall extrusion approaches for the cached PLATEAU Tokyo Tower.

Run with QT_QPA_PLATFORM=offscreen and the project environment. This isolated
prototype assumes a ground elevation of 18.4 m and draws only tower walls.
"""

import argparse
from pathlib import Path
from dataclasses import replace
from statistics import mean
from PySide6.QtGui import QImage, QPainter, QColor, QFont
from PySide6.QtGui import QGuiApplication
from zstarview.urban_outline_layer import parse_derived_tile_buildings, _nearby_wall_faces, _project_observer_centric_urban_outline_result
from zstarview.data.urban_outline_common import BuildingFootprint
from zstarview.types import ViewerData, ScreenGeometry
from zstarview.render.terrain import draw_urban_outlines
from types import SimpleNamespace

app=QGuiApplication([])
parser = argparse.ArgumentParser()
parser.add_argument('tile', type=Path, help='PLATEAU derived tile containing Tokyo Tower')
parser.add_argument('--output', type=Path, default=Path('tokyo-tower-comparison.png'))
args = parser.parse_args()
path = args.tile
tower=next(b for b in parse_derived_tile_buildings(path) if '7aff4a51-be8b-405b-abe4-ac489697cbc8' in b.building_id)
lat,lon=35.65691966858282,139.74599240924002
ground=18.4
observer=ground+1.7
levels=sorted({round(mean(p[2] for p in ring),1) for ring in tower.roof_surfaces_lonlat})
print('levels (absolute metres):', levels, flush=True)
tower=replace(tower, ground_elevation_m=ground)
base=_nearby_wall_faces((tower,),lat_deg=lat,lon_deg=lon,observer_elevation_m=observer)
parts=[replace(tower, height_m=levels[0]-ground, roof_surfaces_lonlat=())]
for ring in tower.roof_surfaces_lonlat:
 h=round(mean(p[2] for p in ring),1)
 lower=levels[levels.index(h)-1] if levels.index(h) else ground
 if h-lower < 0.1: continue
 parts.append(BuildingFootprint(building_id='trial', height_m=h-ground, min_height_m=max(0,lower-ground), ground_elevation_m=ground, rings_lonlat=(tuple((x,y) for x,y,_ in ring),)))
trial=_nearby_wall_faces(tuple(parts),lat_deg=lat,lon_deg=lon,observer_elevation_m=observer)
print('patches existing / trial:',len(base),len(trial),flush=True)
viewer=ViewerData(location=(lat,lon),timezone_name='Asia/Tokyo',city_name='Tokyo Tower',observer_height_m=1.7,ground_elevation_m=ground,view_center=(5.0,345.0))
image=QImage(1200,600,QImage.Format.Format_ARGB32_Premultiplied);image.fill(QColor('#101821'))
painter=QPainter(image)
for i,(label,walls) in enumerate((('Current: footprint to full height',base),('Trial: roof rings to next level',trial))):
 p=_project_observer_centric_urban_outline_result(SimpleNamespace(outlines=(),walls=walls),view_center=viewer.view_center,edge_fov_deg=viewer.edge_fov_deg)
 painter.save();painter.translate(600*i,0)
 draw_urban_outlines(painter,ScreenGeometry(center=(300,300),radius=300),viewer,p or [],opacity=0.7)
 painter.setPen(QColor('white'));painter.setFont(QFont('Sans Serif',11));painter.drawText(12,24,label)
 painter.restore()
painter.end();image.save(str(args.output))
