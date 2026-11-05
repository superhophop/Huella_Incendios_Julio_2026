"""Build the data for the public web map in docs/ (served by GitHub Pages).

Run after notebooks 04-08, from the notebooks folder:
    ..\\.venv\\Scripts\\python build_webmap.py

It writes small web-friendly files into docs/data/:
  footprint.geojson   simplified burn footprint (WGS84)
  hotspots.json       NASA FIRMS hotspots: lon, lat, day number, local time, power
  day_01.png ...      burned area up to each day, coloured by day (for the slider)
  summary.json        key figures and the image bounds
"""
import json

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import shapely
from PIL import Image
from rasterio.warp import Resampling, calculate_default_transform, reproject, transform_bounds

from common import OUT, PROJECT, load_json

DOCS = PROJECT / "docs" / "data"
DOCS.mkdir(parents=True, exist_ok=True)

# Day colours, from yellow (first day) to dark red (last day). The page uses the same list.
COLORS = ["#ffe066", "#fec44f", "#fe9929", "#f46d2a", "#e8452c", "#d7301f", "#b81d2e", "#99003f", "#7a0177", "#5c0060"]

prog = load_json("step8_progression.json", "08_progression")
val = load_json("step7_validation.json", "07_validation_effis")
imp = load_json("step5_impact.json", "05_impact")
fire_day = pd.Timestamp(prog["daily"][0]["date"])

# 1. Footprint: simplify to ~25 m so the file stays small, then convert to latitude/longitude.
fp = gpd.read_file(OUT / "burn_footprint.gpkg", layer="footprint_dissolved")
fp["geometry"] = fp.geometry.simplify(30)
fp = fp.to_crs(4326)[["area_km2", "geometry"]]
fp["geometry"] = shapely.set_precision(fp.geometry.values, 1e-5)   # 5 decimals is about 1 m
(DOCS / "footprint.geojson").write_text(fp.to_json(drop_id=True, to_wgs84=True), encoding="utf-8")

# 2. Hotspots as compact rows [lon, lat, day number, "HH:MM", power in MW].
hs = gpd.read_file(OUT / "firms_hotspots.gpkg").to_crs(4326)
local = pd.to_datetime(hs["local"].astype(str).str[:16])
day_no = (local.dt.normalize() - fire_day).dt.days + 1
rows = [[round(p.x, 4), round(p.y, 4), int(d), t.strftime("%H:%M"), round(float(f), 1)]
        for p, d, t, f in zip(hs.geometry, day_no, local, hs["frp"])]
(DOCS / "hotspots.json").write_text(json.dumps(rows, separators=(",", ":")), encoding="utf-8")

# 3. Burn date raster -> Web Mercator (the projection web maps use), one PNG per day.
with rasterio.open(OUT / "burn_date_50m.tif") as src:
    dst_crs = "EPSG:3857"
    tr, w, h = calculate_default_transform(src.crs, dst_crs, src.width, src.height, *src.bounds)
    day = np.zeros((h, w), "uint8")
    reproject(rasterio.band(src, 1), day, src_transform=src.transform, src_crs=src.crs,
              dst_transform=tr, dst_crs=dst_crs, resampling=Resampling.nearest, src_nodata=0, dst_nodata=0)
    west, south, east, north = transform_bounds(dst_crs, "EPSG:4326", tr.c, tr.f + tr.e * h, tr.c + tr.a * w, tr.f)

rgb = np.array([[int(c[i:i + 2], 16) for i in (1, 3, 5)] for c in COLORS], "uint8")
n_days = int(day.max())
for k in range(1, n_days + 1):
    img = np.zeros((h, w, 4), "uint8")
    shown = (day > 0) & (day <= k)
    img[shown, :3] = rgb[day[shown] - 1]
    img[shown, 3] = 215
    Image.fromarray(img, "RGBA").save(DOCS / f"day_{k:02d}.png", optimize=True)

# 4. Key figures for the cards and the slider.
counts = pd.Series(day_no).value_counts().to_dict()
burned = {d["date"]: d for d in prog["daily"]}      # days with no new area are missing from this list
daily, cum = [], 0.0
for k in range(1, n_days + 1):
    date = str((fire_day + pd.Timedelta(days=k - 1)).date())
    km2 = burned.get(date, {}).get("burned_km2", 0.0)
    cum = burned.get(date, {}).get("cumulative_km2", cum)
    daily.append({"date": date, "km2": km2, "cum_km2": cum, "hotspots": int(counts.get(k, 0))})
summary = {
    "fire_day": str(fire_day.date()),
    "days": n_days,
    "colors": COLORS[:n_days],
    "bounds": [[round(south, 5), round(west, 5)], [round(north, 5), round(east, 5)]],
    "daily": daily,
    "footprint_km2": round(float(fp["area_km2"].sum()), 2),
    "provinces": {k: v["burned_km2"] for k, v in imp["provinces"].items()},
    "municipalities": imp["overlays"]["líneas límite municipales"]["features_affected"],
    "natura_pct": round(100 * imp["overlays"]["Natura 2000"]["burned_km2"] / float(fp["area_km2"].sum()), 1),
    "dice": val["dice"],
    "effis_km2": val["effis_km2"],
    "peak_day": prog["peak_day"],
    "peak_km2": prog["peak_day_km2"],
    "first_detection": prog["first_detection_local"],
    "first_lonlat": prog["first_detection_lonlat"],
    "hotspots_total": len(rows),
}
(DOCS / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding="utf-8")

size = sum(f.stat().st_size for f in DOCS.iterdir()) / 1e6
print(f"Wrote {len(list(DOCS.iterdir()))} files to docs/data ({size:.1f} MB)")
