"""Extrude each cached Tokyo Tower roof polygon down to the highest lower floor
that lies directly beneath it (overlapping in plan), or to the ground."""
import json, math, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.patches import Polygon
from matplotlib.path import Path as MplPath
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

tile = Path.home() / ".cache/zstarview/plateau_buildings/13103/bldg/53393599_bldg_6697_op.json"
out = Path(sys.argv[1])
GROUND = 18.4
b = next(x for x in json.loads(tile.read_text())["buildings"] if "7aff4a51-be8b-405b-abe4-ac489697cbc8" in x["id"])
ring = b["rings"][0]
lon0 = sum(p[0] for p in ring[:-1]) / (len(ring) - 1); lat0 = sum(p[1] for p in ring[:-1]) / (len(ring) - 1)
kx = 111320.0 * math.cos(math.radians(lat0)); ky = 110574.0
en = lambda lon, lat: ((lon - lon0) * kx, (lat - lat0) * ky)
foot = [en(lon, lat) for lon, lat in ring]
top_h = b["height_m"]

roofs = []
for s in b["roof_surfaces"]:
    pts = [(*en(lon, lat), z - GROUND) for lon, lat, z in s]
    if pts[0] == pts[-1]:
        pts = pts[:-1]
    level = round(sum(p[2] for p in pts) / len(pts))
    roofs.append((level, [(p[0], p[1]) for p in pts]))
levels = sorted({lv for lv, _ in roofs})
MIN_OVERLAP_M2 = 0.5  # ignore polygons that only touch along an edge


def overlap_area(a, b, steps=40):
    """Estimate the plan-view intersection area of two polygons by grid sampling."""
    xs = [p[0] for p in a] + [p[0] for p in b]
    ys = [p[1] for p in a] + [p[1] for p in b]
    x0 = max(min(p[0] for p in a), min(p[0] for p in b)); x1 = min(max(p[0] for p in a), max(p[0] for p in b))
    y0 = max(min(p[1] for p in a), min(p[1] for p in b)); y1 = min(max(p[1] for p in a), max(p[1] for p in b))
    if x1 <= x0 or y1 <= y0:
        return 0.0
    dx = (x1 - x0) / steps; dy = (y1 - y0) / steps
    pts = [(x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy) for i in range(steps) for j in range(steps)]
    pa, pb = MplPath(a), MplPath(b)
    inside = sum(1 for p, ia in zip(pts, pa.contains_points(pts)) if ia and pb.contains_point(p))
    return inside * dx * dy


prisms = []  # (bottom, top, outline)
for lv, xy in roofs:
    below = [olv for olv, oxy in roofs if olv < lv and overlap_area(xy, oxy) >= MIN_OVERLAP_M2]
    prisms.append((max(below) if below else 0, lv, xy))
print("levels:", levels)
for lv in levels:
    group = [p for p in prisms if p[1] == lv]
    bottoms = {}
    for bottom, _top, _xy in group:
        bottoms[bottom] = bottoms.get(bottom, 0) + 1
    text = ", ".join(f"{n} to {b} m" for b, n in sorted(bottoms.items()))
    print(f"  {lv:4d} m: {len(group):2d} polygons -> {text}")

norm = colors.Normalize(0, max(levels)); cmap = cm.viridis
(x0, y0), (x1, y1) = foot[0], foot[1]
ux, uy = x1 - x0, y1 - y0; L = math.hypot(ux, uy); ux, uy = ux / L, uy / L
proj = lambda x, y: x * ux + y * uy

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(1, 2, 1)
s_foot = [proj(*p) for p in foot]
ax.add_patch(Polygon([(min(s_foot), 0), (max(s_foot), 0), (max(s_foot), top_h), (min(s_foot), top_h)],
                     fc="none", ec="#e4572e", lw=1, ls="--", label="current approximation (LOD1 ring x measuredHeight)"))
for bottom, top, xy in sorted(prisms, key=lambda p: p[1]):
    s = [proj(*p) for p in xy]
    ax.add_patch(Polygon([(min(s), bottom), (max(s), bottom), (max(s), top), (min(s), top)],
                         fc=cmap(norm(top)), ec="k", lw=0.3, alpha=0.35))
for lv in levels:
    ax.axhline(lv, color="gray", lw=0.4, ls=":")
    ax.text(max(s_foot) + 3, lv, f"{lv} m", va="center", fontsize=8, color="gray")
ax.set_xlim(min(s_foot) - 12, max(s_foot) + 16); ax.set_ylim(-5, top_h + 10)
ax.set_xlabel("Along footprint edge (m)"); ax.set_ylabel(f"Height above assumed ground {GROUND} m (m)")
ax.set_title("Side projection: extruded to the floor directly beneath (or ground)")
ax.legend(loc="upper right", fontsize=8); ax.grid(alpha=0.2)

ax = fig.add_subplot(1, 2, 2, projection="3d")
faces, fcols = [], []
for bottom, top, xy in prisms:
    c = cmap(norm(top))
    n = len(xy)
    for i in range(n):
        a, d = xy[i], xy[(i + 1) % n]
        faces.append([(*a, bottom), (*d, bottom), (*d, top), (*a, top)]); fcols.append(c)
    faces.append([(*p, top) for p in xy]); fcols.append(c)
ax.add_collection3d(Poly3DCollection(faces, facecolors=fcols, edgecolors=(0, 0, 0, 0.25), linewidths=0.2, alpha=0.55))
ax.set_xlim(-60, 60); ax.set_ylim(-60, 60); ax.set_zlim(0, top_h + 10)
ax.set_box_aspect((1, 1, 2.2)); ax.view_init(elev=10, azim=-60)
ax.set_title("3D view of the extruded prisms")
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=fig.axes, shrink=0.6, label="Floor (top) height (m)")
fig.savefig(out, dpi=110, bbox_inches="tight")
