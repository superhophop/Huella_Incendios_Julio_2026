#!/usr/bin/env python3
# =============================================================================
#  s2_burn_footprint.py
#  Burn-footprint mapping with Sentinel-2, explained step by step
#
#  Case study : July 2026 Ávila–Madrid–Toledo fire, Sierra Oeste de Madrid AOI
#  Satellite  : Sentinel-2 (ESA Copernicus), Level-2A surface reflectance
#  Method     : change detection between a pre-fire and a post-fire image
#               (Key & Benson 2006; Chuvieco 2015; Martín 1998;
#                Pinty & Verstraete 1992; Bastarrika et al. 2011)
# =============================================================================
#
#  HOW TO RUN IT
#  -------------
#  Each step can be run on its own. A step reads what the previous step saved
#  in the output folder, so you can stop, look at the results, and continue.
#
#     python s2_burn_footprint.py --step 1    # find the cloud-free images
#     python s2_burn_footprint.py --step 2    # download bands -> reflectance
#     python s2_burn_footprint.py --step 3    # spectral indices
#     python s2_burn_footprint.py --step 4    # thresholds -> burn footprint
#     python s2_burn_footprint.py --step 5    # impact (severity, land cover...)
#     python s2_burn_footprint.py --step 6    # fire weather vs climate
#     python s2_burn_footprint.py --step all  # everything in one go
#
#  WHERE THE DATA COMES FROM (all free, no account or API key)
#  ----------------------------------------------------------
#  * Microsoft Planetary Computer -> Sentinel-2 L2A images and ESA WorldCover.
#    It is searched through a STAC API ("SpatioTemporal Asset Catalog"): a
#    standard web interface where you ask "which images cover this box between
#    these dates?" and get back a list of items, each with links to its bands
#    stored as Cloud-Optimised GeoTIFFs. Only the pixels inside the AOI are
#    downloaded, not whole 100 x 100 km tiles.
#  * Open-Meteo archive API -> ERA5 reanalysis weather (temperature, humidity,
#    wind, rain) for the AOI centre, plus the 1991-2020 climate normal.
#
#  LIBRARIES (installed from PyPI with pip, not from GitHub)
#  ---------------------------------------------------------
#  pystac-client + planetary-computer : search the STAC catalogue, sign links
#  odc-stac   : read only the AOI window of each band, on a common 10 m grid
#  numpy/scipy: all the maths (indices, thresholds, filters)
#  rasterio   : read/write GeoTIFFs; turn the raster footprint into polygons
#  geopandas/shapely : vector footprint, areas, perimeters, GeoPackage
#  matplotlib : PNG maps and charts   |  requests : weather API calls
# =============================================================================

from __future__ import annotations

import argparse
import json
import sys
import warnings
from datetime import date, timedelta
from pathlib import Path

import numpy as np

# =============================================================================
#  SETTINGS - change these to study another fire or area
# =============================================================================
CFG = {
    # Area of interest (lon_min, lat_min, lon_max, lat_max), WGS84 degrees.
    # Sierra Oeste de Madrid: San Martín de Valdeiglesias, Pelayos de la Presa,
    # Navas del Rey, Chapinería, Villa del Prado, Aldea del Fresno.
    "bbox": [-4.46, 40.26, -4.17, 40.46],
    "fire_start": "2026-07-22",   # first ignition (Almorox, Toledo)
    "days_before": 28,            # pre-fire image target = ignition - 28 days
    "pre_window_days": 7,         # search +/- 7 days around that target ...
    "pre_max_window_days": 21,    # ... widening up to +/- 21 if nothing is clear
    "post_start": "2026-07-30",   # after the fire was out (Key & Benson
    "post_end": "2026-09-15",     #   "initial assessment")
    "max_tile_cloud": 40,         # quick pre-filter on the whole tile's cloud %
    "max_aoi_cloud": 1.0,         # max % of AOI pixels that may be cloud
    "min_aoi_coverage": 99.0,     # min % of AOI that must have data
    "crs": "EPSG:32630",          # WGS84 / UTM zone 30N (metres)
    "resolution": 10,             # metres
    "k_mad": 3.0,                 # threshold = unburned median +/- 3 robust sd
    "pixel_votes": 2,             # indices that must flag a pixel
    "core_votes": 3,              # indices that must agree in a patch core
    "min_patch_ha": 1.0,          # smallest patch / hole kept
    "out": "outputs_sentinel2",
}

COLLECTION = "sentinel-2-l2a"
STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
OPEN_METEO = "https://archive-api.open-meteo.com/v1/archive"

# Sentinel-2 bands we need, and why
BANDS = {
    "blue":  "B02",   # 490 nm, 10 m  - only for display
    "green": "B03",   # 560 nm, 10 m  - false-colour composite
    "red":   "B04",   # 665 nm, 10 m  - chlorophyll absorbs red
    "nir10": "B08",   # 842 nm, 10 m  - healthy leaves reflect strongly; NDVI/BAI/GEMI
    "nir":   "B8A",   # 865 nm, 20 m  - narrow NIR that pairs with B12 for NBR
    "swir2": "B12",   # 2190 nm, 20 m - burned/dry surfaces reflect strongly
}
QA_BAND = "SCL"       # Scene Classification Layer, 20 m (cloud, shadow, water...)

# SCL codes: 0 no data, 1 saturated, 3 cloud shadow, 6 water, 8/9 cloud,
# 10 thin cirrus, 11 snow. Class 3 is deliberately NOT masked: Sentinel-2's
# own classifier often labels fresh, very dark burn scars as "cloud shadow".
SCL_NODATA, SCL_CLOUD, SCL_SNOW, SCL_WATER = [0, 1], [8, 9, 10], [11], [6]

# Strong-char check for patch cores: the course's ΔBAI value (Bastarrika 2011)
STRONG_CHAR_DBAI = 144.3835

SEVERITY_BINS = [-0.25, -0.10, 0.10, 0.27, 0.44, 0.66]       # USGS dNBR classes
SEVERITY_LABELS = ["Enhanced regrowth, high", "Enhanced regrowth, low", "Unburned",
                   "Low severity", "Moderate-low severity", "Moderate-high severity",
                   "High severity"]
SEVERITY_COLORS = ["#1A9850", "#91CF60", "#D9D9D9", "#FEE08B", "#FDAE61", "#F46D43", "#A50026"]
NDVI_BINS = [0.0, 0.2, 0.4, 0.6]                              # Chuvieco (2015)
NDVI_LABELS = ["Water / shadow (<0)", "Bare soil (0-0.2)", "Sparse veg + soil (0.2-0.4)",
               "Dispersed vegetation (0.4-0.6)", "Dense vegetation (>0.6)"]
NDVI_COLORS = ["#000000", "#FFFFCC", "#D3FFBE", "#4CE600", "#267300"]
WORLDCOVER = {10: "Tree cover", 20: "Shrubland", 30: "Grassland", 40: "Cropland",
              50: "Built-up", 60: "Bare / sparse vegetation", 70: "Snow and ice",
              80: "Permanent water", 90: "Herbaceous wetland", 95: "Mangroves",
              100: "Moss and lichen"}


def say(msg: str = "") -> None:
    print(msg, flush=True)


def header(n, title: str) -> None:
    say("\n" + "=" * 72 + f"\n STEP {n} - {title}\n" + "=" * 72)


def out_dir() -> Path:
    p = Path(CFG["out"])
    p.mkdir(parents=True, exist_ok=True)
    return p


# =============================================================================
#  Small GeoTIFF helpers (rasterio)
#  A GeoTIFF = pixel values + "where on Earth" (CRS + transform). The transform
#  says: top-left corner X/Y and pixel size, so ArcGIS can place every pixel.
# =============================================================================
def write_tif(path: Path, arrays, transform, crs, nodata=None, names=None):
    import rasterio
    arrays = [arrays] if np.ndim(arrays) == 2 else list(arrays)
    h, w = arrays[0].shape
    with rasterio.open(path, "w", driver="GTiff", height=h, width=w, count=len(arrays),
                       dtype=arrays[0].dtype, crs=crs, transform=transform,
                       nodata=nodata, compress="deflate") as dst:
        for i, a in enumerate(arrays, 1):
            dst.write(a, i)
            if names:
                dst.set_band_description(i, names[i - 1])


def read_tif(path: Path):
    """Return (array or stack, transform, crs, band names)."""
    import rasterio
    with rasterio.open(path) as src:
        data = src.read()
        return (data[0] if src.count == 1 else data), src.transform, src.crs, list(src.descriptions)


def need(path: Path, step: int) -> Path:
    if not path.exists():
        sys.exit(f"Missing {path}. Run --step {step} first.")
    return path


# =============================================================================
#  STEP 1 - FIND THE CLOUD-FREE IMAGES
# =============================================================================
#  Why: a burn footprint is a comparison of two dates. Both images must be
#  free of cloud *over the AOI*. The catalogue's own "cloud %" describes the
#  whole 100 x 100 km tile, so a 5 % cloudy tile can still have its cloud right
#  on top of our area. We therefore download the small Scene Classification
#  Layer (SCL) for every candidate date and count cloud pixels inside the AOI.
#
#  Pre-fire: the clear date closest to ignition - 28 days.
#  Post-fire: the earliest clear date after the fire was out (less regrowth).
# =============================================================================
def open_catalog():
    import planetary_computer
    import pystac_client
    # sign_inplace adds a temporary access token to every band link
    return pystac_client.Client.open(STAC_URL, modifier=planetary_computer.sign_inplace)


def search_by_day(catalog, start: str, end: str) -> dict:
    """STAC search -> {date: [items]}. The AOI usually needs 2 tiles on one date."""
    items = catalog.search(collections=[COLLECTION], bbox=CFG["bbox"],
                           datetime=f"{start}/{end}",
                           query={"eo:cloud_cover": {"lt": CFG["max_tile_cloud"]}}).items()
    by_day: dict = {}
    for it in items:
        by_day.setdefault(it.datetime.date().isoformat(), []).append(it)
    return dict(sorted(by_day.items()))


def score_day(items) -> dict:
    """Load only the SCL band at 60 m over the AOI and measure real cloud cover."""
    import odc.stac
    ds = odc.stac.load(items, bands=[QA_BAND], bbox=CFG["bbox"], crs=CFG["crs"],
                       resolution=60, groupby="solar_day", resampling="nearest")
    scl = ds[QA_BAND].isel(time=0).values
    has_data = ~np.isin(scl, SCL_NODATA)
    cover = 100 * has_data.mean()
    cloud = 100 * np.isin(scl[has_data], SCL_CLOUD).mean() if has_data.any() else 100
    return {"coverage_pct": round(float(cover), 2), "aoi_cloud_pct": round(float(cloud), 3)}


def rank(by_day: dict) -> list:
    rows = []
    for d, items in by_day.items():
        try:
            s = score_day(items)
        except Exception as e:                      # a broken file, network hiccup...
            say(f"   {d}: could not be read ({e})")
            continue
        s.update(date=d, ids=[i.id for i in items],
                 satellite=sorted({i.properties.get("platform", "?") for i in items}))
        s["clear"] = (s["aoi_cloud_pct"] <= CFG["max_aoi_cloud"]
                      and s["coverage_pct"] >= CFG["min_aoi_coverage"])
        say(f"   {d}  {','.join(s['satellite']):12s} cloud {s['aoi_cloud_pct']:6.2f} %  "
            f"coverage {s['coverage_pct']:6.1f} %  {'CLEAR' if s['clear'] else ''}")
        rows.append(s)
    return rows


def step1():
    header(1, "FIND CLOUD-FREE IMAGES")
    cat = open_catalog()
    fire = date.fromisoformat(CFG["fire_start"])
    target = fire - timedelta(days=CFG["days_before"])
    say(f"Ignition {fire}; pre-fire target = {target} ({CFG['days_before']} days before)")

    pre = None
    half = CFG["pre_window_days"]
    while pre is None and half <= CFG["pre_max_window_days"]:
        start, end = target - timedelta(days=half), min(target + timedelta(days=half),
                                                        fire - timedelta(days=1))
        say(f"\nPRE-FIRE candidates {start} -> {end}:")
        clear = [r for r in rank(search_by_day(cat, str(start), str(end))) if r["clear"]]
        if clear:
            pre = min(clear, key=lambda r: abs((date.fromisoformat(r["date"]) - target).days))
        else:
            half += 7
            say(f"   none clear - widening the window to +/- {half} days")
    if pre is None:
        sys.exit("No cloud-free pre-fire image. Relax max_aoi_cloud in CFG.")

    say(f"\nPOST-FIRE candidates {CFG['post_start']} -> {CFG['post_end']}:")
    clear = [r for r in rank(search_by_day(cat, CFG["post_start"], CFG["post_end"])) if r["clear"]]
    if not clear:
        sys.exit("No cloud-free post-fire image. Extend post_end in CFG.")
    post = min(clear, key=lambda r: r["date"])

    lag_pre = (fire - date.fromisoformat(pre["date"])).days
    lag_post = (date.fromisoformat(post["date"]) - fire).days
    result = {"ignition": str(fire), "pre_target": str(target),
              "pre": {**pre, "days_before_ignition": lag_pre},
              "post": {**post, "days_after_ignition": lag_post}}
    (out_dir() / "step1_scenes.json").write_text(json.dumps(result, indent=2))
    say(f"\nSELECTED  pre-fire  {pre['date']} ({lag_pre} days before ignition)  {pre['satellite']}")
    say(f"SELECTED  post-fire {post['date']} ({lag_post} days after ignition)  {post['satellite']}")
    say(f"Saved -> {out_dir() / 'step1_scenes.json'}")


# =============================================================================
#  STEP 2 - DOWNLOAD THE BANDS AND CONVERT THEM TO REFLECTANCE
# =============================================================================
#  The images store "digital numbers" (DN), whole numbers that save space.
#  Sentinel-2 L2A (processing baseline 04.00+, i.e. since 2022) stores
#        DN = reflectance x 10000 + 1000      ->   reflectance = (DN - 1000) / 10000
#  Reflectance (0-1) = fraction of sunlight the ground reflects, which is what
#  the indices need. This is the same idea as the course's Landsat formula
#  DN x 0.0000275 - 0.2, just with Sentinel-2's own scale and offset.
#
#  All bands are put on one 10 m grid (the 20 m bands B8A, B12, SCL are
#  resampled), and the post-fire image is snapped onto exactly the same pixels
#  as the pre-fire one so the two can be subtracted pixel by pixel.
#  A "valid" mask removes cloud, snow, no-data and water (the reservoirs).
# =============================================================================
def s2_offset(items) -> int:
    vals = []
    for it in items:
        try:
            vals.append(float(it.properties.get("s2:processing_baseline")))
        except (TypeError, ValueError):
            pass
    return -1000 if (not vals or max(vals) >= 4.0) else 0


def to_reflectance(dn: np.ndarray, offset: int) -> np.ndarray:
    dn = dn.astype("float32")
    r = (dn + offset) / 10000.0
    r[dn == 0] = np.nan                  # 0 = no data
    return np.clip(r, 0, None)           # tiny negatives -> 0 (course: Con(x<0,0,x))


def load_date(items, geobox=None):
    import odc.stac
    names = list(BANDS.values()) + [QA_BAND]
    kw = dict(bands=names, groupby="solar_day",
              resampling={QA_BAND: "nearest", "*": "bilinear"})
    if geobox is None:
        kw.update(bbox=CFG["bbox"], crs=CFG["crs"], resolution=CFG["resolution"])
    else:
        kw.update(geobox=geobox)
    ds = odc.stac.load(items, **kw).isel(time=0)
    off = s2_offset(items)
    refl = {role: to_reflectance(ds[b].values, off) for role, b in BANDS.items()}
    return ds, refl, ds[QA_BAND].values


def stretch_rgb(rgb: np.ndarray, lo, hi) -> np.ndarray:
    return np.nan_to_num(np.clip((rgb - lo) / (hi - lo), 0, 1))


def step2():
    header(2, "DOWNLOAD BANDS -> REFLECTANCE")
    import matplotlib.pyplot as plt
    scenes = json.loads(need(out_dir() / "step1_scenes.json", 1).read_text())
    cat = open_catalog()
    items = {k: list(cat.search(collections=[COLLECTION], ids=scenes[k]["ids"]).items())
             for k in ("pre", "post")}

    say(f"Reading {len(BANDS) + 1} bands for {scenes['pre']['date']} (pre) ...")
    ds, pre, scl0 = load_date(items["pre"])
    gb = ds.odc.geobox
    transform = getattr(gb, "transform", None) or gb.affine
    say(f"Reading the same bands for {scenes['post']['date']} (post) on the same grid ...")
    _, post, scl1 = load_date(items["post"], geobox=gb)

    valid = ~(np.isin(scl0, SCL_NODATA + SCL_CLOUD + SCL_SNOW)
              | np.isin(scl1, SCL_NODATA + SCL_CLOUD + SCL_SNOW)
              | np.isin(scl0, SCL_WATER))
    for b in ("red", "nir10", "nir", "swir2"):
        valid &= np.isfinite(pre[b]) & np.isfinite(post[b])

    o = out_dir()
    roles = list(BANDS)
    write_tif(o / "pre_reflectance.tif", [pre[r] for r in roles], transform, CFG["crs"],
              np.nan, [f"{r} ({BANDS[r]})" for r in roles])
    write_tif(o / "post_reflectance.tif", [post[r] for r in roles], transform, CFG["crs"],
              np.nan, [f"{r} ({BANDS[r]})" for r in roles])
    write_tif(o / "valid_mask.tif", valid.astype("uint8"), transform, CFG["crs"], 255)

    # False-colour composites for ArcGIS: band 1 = NIR (B8), 2 = Red (B4), 3 = Green (B3)
    for name, b in (("pre", pre), ("post", post)):
        write_tif(o / f"{name}_fire_843.tif", [b["nir10"], b["red"], b["green"]],
                  transform, CFG["crs"], np.nan, ["NIR B8", "Red B4", "Green B3"])
    write_tif(o / "nir_change_composite.tif", [pre["nir10"], post["nir10"], pre["nir10"]],
              transform, CFG["crs"], np.nan, ["NIR pre", "NIR post", "NIR pre"])

    # One shared stretch (from the pre-fire image) so both dates are comparable
    rgb0 = np.dstack([pre["nir10"], pre["red"], pre["green"]])
    rgb1 = np.dstack([post["nir10"], post["red"], post["green"]])
    lo = np.nanpercentile(rgb0, 2, axis=(0, 1))
    hi = np.nanpercentile(rgb0, 98, axis=(0, 1))
    fig, ax = plt.subplots(1, 2, figsize=(15, 7), constrained_layout=True)
    ax[0].imshow(stretch_rgb(rgb0, lo, hi)); ax[0].set_title(f"Pre-fire {scenes['pre']['date']}")
    ax[1].imshow(stretch_rgb(rgb1, lo, hi)); ax[1].set_title(f"Post-fire {scenes['post']['date']}")
    for a in ax:
        a.set_axis_off()
    fig.suptitle("Sentinel-2 false colour 8-4-3 (NIR, Red, Green): vegetation = red, burn scar = dark")
    fig.savefig(o / "step2_false_color_843.png", dpi=130)
    plt.close(fig)

    h, w = valid.shape
    info = {"grid": {"rows": h, "cols": w, "pixel_m": CFG["resolution"], "crs": CFG["crs"],
                     "aoi_km2": round(h * w * CFG["resolution"] ** 2 / 1e6, 2),
                     "transform": list(transform)[:6]},
            "valid_pct": round(float(valid.mean() * 100), 2),
            "water_pct": round(float(np.isin(scl0, SCL_WATER).mean() * 100), 2),
            "arcgis_stretch_843": {"Band_1 NIR": [round(float(lo[0]), 3), round(float(hi[0]), 3)],
                                   "Band_2 Red": [round(float(lo[1]), 3), round(float(hi[1]), 3)],
                                   "Band_3 Green": [round(float(lo[2]), 3), round(float(hi[2]), 3)]}}
    (o / "step2_info.json").write_text(json.dumps(info, indent=2))
    say(f"Grid: {h} x {w} pixels of {CFG['resolution']} m = {info['grid']['aoi_km2']} km²")
    say(f"Valid pixels: {info['valid_pct']} %  (water {info['water_pct']} % removed)")
    say(f"ArcGIS stretch for *_fire_843.tif (Minimum Maximum, Custom): {info['arcgis_stretch_843']}")
    say(f"Saved -> pre/post_reflectance.tif, valid_mask.tif, *_843.tif, step2_false_color_843.png")


# =============================================================================
#  STEP 3 - SPECTRAL INDICES
# =============================================================================
#  Fire changes how the ground reflects light: green leaves (high NIR, low red)
#  become char and ash (low NIR, higher SWIR). Each index condenses that change
#  into one number per pixel. Two-date indices ("Δ", "d") subtract pre from post
#  (or post from pre) so only the CHANGE remains.
#
#   NDVI = (NIR - Red) / (NIR + Red)            greenness           (Rouse 1973)
#   NBR  = (NIR - SWIR2) / (NIR + SWIR2)        burn-sensitive      (Key & Benson)
#   BAI  = 1 / ((0.1 - Red)² + (0.06 - NIR)²)   closeness to char   (Martín 1998)
#   GEMI = η(1 - 0.25η) - (Red - 0.125)/(1 - Red),
#          η = (2(NIR² - Red²) + 1.5NIR + 0.5Red) / (NIR + Red + 0.5)
#                                               vegetation, robust to atmosphere
#                                               (Pinty & Verstraete 1992)
#   ΔNDVI = post - pre   (burned: negative)     dNBR = pre - post (burned: positive)
#   ΔBAI  = post - pre   (burned: positive)     ΔGEMI = post - pre (burned: negative)
#
#  NDVI, BAI and GEMI use B08 (10 m). NBR uses B8A with B12 (20 m bands).
# =============================================================================
def _nd(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        s = a + b
        return np.where(s > 0, (a - b) / s, np.nan).astype("float32")


def ndvi(nir, red):
    return _nd(nir, red)


def nbr(nir, swir2):
    return _nd(nir, swir2)


def bai(nir, red):
    with np.errstate(divide="ignore", invalid="ignore"):
        return (1.0 / ((0.1 - red) ** 2 + (0.06 - nir) ** 2)).astype("float32")


def gemi(nir, red):
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = (2 * (nir ** 2 - red ** 2) + 1.5 * nir + 0.5 * red) / (nir + red + 0.5)
        g = eta * (1 - 0.25 * eta) - (red - 0.125) / (1 - red)
        return np.where(red < 1, g, np.nan).astype("float32")


def compute_indices(pre: dict, post: dict) -> dict:
    n0, n1 = ndvi(pre["nir10"], pre["red"]), ndvi(post["nir10"], post["red"])
    b0, b1 = nbr(pre["nir"], pre["swir2"]), nbr(post["nir"], post["swir2"])
    a0, a1 = bai(pre["nir10"], pre["red"]), bai(post["nir10"], post["red"])
    g0, g1 = gemi(pre["nir10"], pre["red"]), gemi(post["nir10"], post["red"])
    return {"NDVI_pre": n0, "NDVI_post": n1, "dNDVI": n1 - n0,
            "NBR_pre": b0, "NBR_post": b1, "dNBR": b0 - b1,
            "BAI_pre": a0, "BAI_post": a1, "dBAI": a1 - a0,
            "GEMI_pre": g0, "GEMI_post": g1, "dGEMI": g1 - g0}


def read_reflectance(name: str) -> tuple:
    stack, tr, crs, _ = read_tif(need(out_dir() / f"{name}_reflectance.tif", 2))
    return {role: stack[i] for i, role in enumerate(BANDS)}, tr, crs


def step3():
    header(3, "SPECTRAL INDICES")
    import matplotlib.pyplot as plt
    o = out_dir()
    pre, tr, crs = read_reflectance("pre")
    post, _, _ = read_reflectance("post")
    valid = read_tif(need(o / "valid_mask.tif", 2))[0] == 1
    idx = compute_indices(pre, post)
    for k, v in idx.items():
        v = np.where(valid, v, np.nan).astype("float32")
        write_tif(o / f"{k}.tif", v, tr, crs, np.nan, [k])
        idx[k] = v
        say(f"   {k:10s} min {np.nanmin(v):9.3f}   median {np.nanmedian(v):9.3f}   "
            f"max {np.nanmax(v):9.3f}")

    panels = [("dNDVI", "ΔNDVI (post - pre)", "RdYlGn", -0.6, 0.3),
              ("dGEMI", "ΔGEMI (post - pre)", "RdYlGn", -0.4, 0.15),
              ("dBAI", "ΔBAI (post - pre)", "RdYlGn_r", -50, 400),
              ("dNBR", "dNBR (pre - post)", "RdYlGn_r", -0.3, 1.0)]
    fig, axs = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)
    for a, (k, t, cm, vmin, vmax) in zip(axs.ravel(), panels):
        im = a.imshow(idx[k], cmap=cm, vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=a, shrink=.8)
        a.set_title(t + "   (red = burned)")
        a.set_axis_off()
    fig.savefig(o / "step3_indices.png", dpi=120)
    plt.close(fig)
    say("Saved -> 12 index GeoTIFFs + step3_indices.png")


# =============================================================================
#  STEP 4 - THRESHOLDS -> BURN FOOTPRINT
# =============================================================================
#  4a. A THRESHOLD for each index decides "burned / not burned" per pixel.
#      The course used fixed values (Bastarrika 2011). They come from a dense
#      shrub fire and miss light grass burns, so here each threshold is learned
#      from THIS image pair: we take pixels the USGS dNBR scale calls unburned
#      (-0.1 to 0.1), measure how much each index normally varies there, and call
#      "burned" anything more than 3 robust standard deviations (MAD) beyond the
#      normal. MAD is used because a few extreme pixels (shadows, reservoir
#      edges) would inflate an ordinary standard deviation.
#  4b. VOTING (automates the course's Merge + manual clean-up):
#      - a pixel is burned if >= 2 of ΔNDVI, ΔGEMI, ΔBAI, dNBR flag it;
#      - a burned patch is kept only if it has >= 1 ha where >= 3 agree AND
#        there is strong char (ΔBAI > 144.38). Harvested cereal also loses
#        greenness, but has no char, so isolated fields are dropped.
#  4c. MAJORITY FILTER (course, 8 neighbours) removes salt-and-pepper noise;
#      patches and holes < 1 ha are removed / filled.
#  4d. VALIDATION: the M index, M = |mean_burned - mean_unburned| /
#      (sd_burned + sd_unburned). M > 1 = the classes separate well.
#  4e. RASTER TO POLYGON + DISSOLVE -> area (ha, km²) and perimeter, saved as
#      GeoPackage, GeoJSON and KML (Google Earth).
# =============================================================================
DIRECTION = {"dNDVI": "<", "dGEMI": "<", "dBAI": ">", "dNBR": ">"}


def scene_thresholds(idx: dict, valid: np.ndarray, k: float) -> dict:
    d = idx["dNBR"]
    with np.errstate(invalid="ignore"):
        ref = valid & np.isfinite(d) & (d > -0.1) & (d < 0.1)
    out = {}
    for name, op in DIRECTION.items():
        v = idx[name][ref]
        v = v[np.isfinite(v)]
        med = float(np.median(v))
        rsd = 1.4826 * float(np.median(np.abs(v - med)))     # robust standard deviation
        thr = med - k * rsd if op == "<" else med + k * rsd
        out[name] = {"op": op, "threshold": round(thr, 5), "unburned_median": round(med, 5),
                     "unburned_robust_sd": round(rsd, 5)}
    return out


def flag(values, op, thr):
    with np.errstate(invalid="ignore"):
        return (values < thr) if op == "<" else (values > thr)


def majority_filter(mask: np.ndarray) -> np.ndarray:
    """ArcGIS Majority Filter (8 neighbours, MAJORITY): a pixel takes the value
    shared by >= 5 of its 8 neighbours."""
    from scipy import ndimage
    k = np.ones((3, 3), int); k[1, 1] = 0
    n = ndimage.convolve(mask.astype(int), k, mode="nearest")
    return np.where(n >= 5, True, np.where(8 - n >= 5, False, mask))


def _label(mask):
    from scipy import ndimage
    return ndimage.label(mask, structure=np.ones((3, 3)))


def remove_small(mask: np.ndarray, min_px: int) -> np.ndarray:
    lab, n = _label(mask)
    if n == 0:
        return mask
    keep = np.bincount(lab.ravel()) >= min_px
    keep[0] = False
    return keep[lab]


def clean(mask, valid, min_px):
    m = remove_small(majority_filter(mask & valid), min_px)
    return ~remove_small(~m, min_px) & valid


def vote_and_clean(seg: dict, strong_char, valid, min_px, pixel_votes, core_votes):
    from scipy import ndimage
    votes = np.sum([m.astype(np.uint8) for m in seg.values()], axis=0)
    cand = majority_filter((votes >= pixel_votes) & valid)
    lab, n = _label(cand)
    if n == 0:
        return cand, votes
    ids = np.arange(1, n + 1)
    core = ndimage.sum((votes >= core_votes) & strong_char & valid, lab, ids)
    keep = np.zeros(n + 1, bool)
    keep[1:] = core >= min_px
    final = ~remove_small(~keep[lab], min_px)
    return final & valid, votes


def separability(values, burned, valid):
    ok = valid & np.isfinite(values)
    b, u = values[burned & ok], values[~burned & ok]
    if b.size < 2 or u.size < 2:
        return None
    return round(float(abs(b.mean() - u.mean()) / (b.std() + u.std())), 3)


def vectorise(mask, transform, crs):
    import geopandas as gpd
    from rasterio import features
    from shapely.geometry import shape
    geoms = [shape(g) for g, v in features.shapes(mask.astype("uint8"), mask=mask,
                                                   transform=transform) if v == 1]
    gdf = gpd.GeoDataFrame(geometry=geoms, crs=crs)
    if len(gdf):
        gdf = gdf.dissolve().reset_index(drop=True)
        gdf["area_ha"] = (gdf.area / 1e4).round(2)
        gdf["area_km2"] = (gdf.area / 1e6).round(4)
    return gdf


def outer_perimeter_km(gdf) -> float:
    from shapely.geometry import Polygon
    if not len(gdf):
        return 0.0
    return round(sum(Polygon(p.exterior).length
                     for p in gdf.explode(index_parts=False).geometry) / 1e3, 2)


def write_kml(gdf, path: Path, name: str):
    from shapely.geometry import Polygon
    g = gdf.to_crs(4326)
    ring = lambda cs: " ".join(f"{x:.6f},{y:.6f},0" for x, y, *_ in cs)  # noqa: E731
    parts = []
    for geom in g.geometry:
        polys = [geom] if isinstance(geom, Polygon) else list(geom.geoms)
        for p in polys:
            inner = "".join(f"<innerBoundaryIs><LinearRing><coordinates>{ring(i.coords)}"
                            "</coordinates></LinearRing></innerBoundaryIs>" for i in p.interiors)
            parts.append(f"<Polygon><outerBoundaryIs><LinearRing><coordinates>"
                         f"{ring(p.exterior.coords)}</coordinates></LinearRing>"
                         f"</outerBoundaryIs>{inner}</Polygon>")
    path.write_text('<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/'
                    f'kml/2.2"><Document><name>{name}</name><Style id="b"><LineStyle><color>'
                    'ff0000ff</color><width>2</width></LineStyle><PolyStyle><color>5a0000ff'
                    '</color></PolyStyle></Style><Placemark><name>' + name + '</name><styleUrl>'
                    '#b</styleUrl><MultiGeometry>' + "".join(parts) +
                    '</MultiGeometry></Placemark></Document></kml>', encoding="utf-8")


def read_indices(names):
    o = out_dir()
    out, tr, crs = {}, None, None
    for k in names:
        out[k], tr, crs, _ = read_tif(need(o / f"{k}.tif", 3))
    return out, tr, crs


def step4():
    header(4, "THRESHOLDS -> BURN FOOTPRINT")
    import geopandas as gpd
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from shapely.geometry import box
    o = out_dir()
    idx, tr, crs = read_indices(list(DIRECTION))
    valid = read_tif(need(o / "valid_mask.tif", 2))[0] == 1
    res = CFG["resolution"]
    px_km2 = res * res / 1e6
    min_px = max(1, int(round(CFG["min_patch_ha"] * 1e4 / res ** 2)))

    thr = scene_thresholds(idx, valid, CFG["k_mad"])
    say("4a. Thresholds learned from unburned pixels (median ± 3·robust sd):")
    seg, rows = {}, []
    for k, t in thr.items():
        seg[k] = flag(idx[k], t["op"], t["threshold"]) & valid
        c = clean(seg[k], valid, min_px)
        m = separability(idx[k], c, valid)
        rows.append({"index": k, "rule": f"{k} {t['op']} {t['threshold']}",
                     "area_km2": round(c.sum() * px_km2, 2), "M": m})
        say(f"   {k:6s} {t['op']} {t['threshold']:9.4f}   -> {rows[-1]['area_km2']:7.1f} km²"
            f"   M = {m}")

    strong = flag(idx["dBAI"], ">", STRONG_CHAR_DBAI) & valid
    final, votes = vote_and_clean(seg, strong, valid, min_px, CFG["pixel_votes"], CFG["core_votes"])
    area_km2 = round(final.sum() * px_km2, 3)
    m_final = separability(idx["dNBR"], final, valid)
    say(f"4b-c. Voting + majority filter -> FINAL footprint {area_km2} km² "
        f"({area_km2 * 100:,.0f} ha)")
    say(f"4d. Validation: dNBR separability inside vs outside footprint M = {m_final}")

    write_tif(o / "footprint_mask.tif", final.astype("uint8"), tr, crs, 255, ["burned"])
    write_tif(o / "index_votes.tif", votes.astype("uint8"), tr, crs, 255, ["votes 0-4"])

    gdf = vectorise(final, tr, crs)
    gpkg = o / "burn_footprint.gpkg"
    if gpkg.exists():
        gpkg.unlink()
    gdf.to_file(gpkg, layer="footprint_dissolved", driver="GPKG")
    h, w = final.shape
    x0, y0 = tr.c, tr.f
    aoi = box(x0, y0 - h * res, x0 + w * res, y0)
    gpd.GeoDataFrame({"area_km2": [round(aoi.area / 1e6, 2)]}, geometry=[aoi], crs=crs
                     ).to_file(gpkg, layer="aoi_extent", driver="GPKG")
    for k, m in seg.items():
        vectorise(clean(m, valid, min_px), tr, crs).to_file(gpkg, layer=f"index_{k}", driver="GPKG")
    gdf.to_crs(4326).to_file(o / "burn_footprint.geojson", driver="GeoJSON")
    write_kml(gdf, o / "burn_footprint.kml", "Burn footprint (Sentinel-2)")
    per_km = outer_perimeter_km(gdf)
    say(f"4e. Polygons saved; outer perimeter {per_km} km "
        "(pixel edges make this longer than an official perimeter)")

    # Figure: footprint over the post-fire false colour
    post, _, _ = read_reflectance("post")
    pre, _, _ = read_reflectance("pre")
    rgb0 = np.dstack([pre["nir10"], pre["red"], pre["green"]])
    lo, hi = np.nanpercentile(rgb0, 2, axis=(0, 1)), np.nanpercentile(rgb0, 98, axis=(0, 1))
    img = stretch_rgb(np.dstack([post["nir10"], post["red"], post["green"]]), lo, hi)
    fig, ax = plt.subplots(figsize=(11, 9), constrained_layout=True)
    ax.imshow(img)
    ax.contour(final, levels=[0.5], colors="yellow", linewidths=0.6)
    ax.set_title(f"Burn footprint: {area_km2:.1f} km² ({final.mean() * 100 / valid.mean():.0f} % "
                 "of valid AOI) over post-fire 8-4-3")
    ax.set_axis_off()
    fig.savefig(o / "step4_footprint.png", dpi=150)
    plt.close(fig)
    fig, axs = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)
    for a, r in zip(axs.ravel(), rows):
        a.imshow(img)
        m = clean(seg[r["index"]], valid, min_px)
        a.imshow(np.ma.masked_where(~m, m), cmap=ListedColormap(["#ff1a1a"]), alpha=.75,
                 interpolation="nearest")
        a.set_title(f"{r['rule']} -> {r['area_km2']} km²  (M = {r['M']})")
        a.set_axis_off()
    fig.savefig(o / "step4_index_segmentations.png", dpi=110)
    plt.close(fig)

    summary = {"thresholds": thr, "per_index": rows,
               "footprint": {"area_km2": area_km2, "area_ha": round(area_km2 * 100, 1),
                             "pct_of_aoi": round(100 * area_km2 / (aoi.area / 1e6), 1),
                             "outer_perimeter_km": per_km, "M_dNBR": m_final},
               "rules": {"pixel_votes": CFG["pixel_votes"], "core_votes": CFG["core_votes"],
                         "strong_char_dBAI": STRONG_CHAR_DBAI,
                         "min_patch_ha": CFG["min_patch_ha"]}}
    (o / "step4_footprint.json").write_text(json.dumps(summary, indent=2))
    say("Saved -> footprint_mask.tif, burn_footprint.gpkg/.geojson/.kml, step4_*.png/json")


# =============================================================================
#  STEP 5 - IMPACT: SEVERITY, VEGETATION, LAND COVER (+ your own layers)
# =============================================================================
#  * Severity: dNBR classes of the USGS / Key & Benson (2006) scale
#    (0.10-0.27 low, 0.27-0.44 moderate-low, 0.44-0.66 moderate-high, >0.66 high).
#  * Vegetation: NDVI zoning of Chuvieco (2015), before vs after.
#  * Land cover: ESA WorldCover 2021 (10 m, from Planetary Computer) - which
#    kind of ground burned: forest, shrub, grass, crops, built-up.
#  * Optional: any vector layer passed with --overlay NAME=PATH (municipalities,
#    Catastro buildings, Natura 2000) -> burned km² or features affected.
# =============================================================================
def classify(values, bins):
    c = np.digitize(values, bins).astype("uint8")
    c[~np.isfinite(values)] = 255
    return c


def class_table(cls, labels, px_km2, inside=None):
    rows = []
    for i, lab in enumerate(labels):
        sel = cls == i
        r = {"class": i, "label": lab, "aoi_km2": round(float(sel.sum() * px_km2), 3)}
        if inside is not None:
            r["in_footprint_km2"] = round(float((sel & inside).sum() * px_km2), 3)
        rows.append(r)
    if inside is not None:
        tot = sum(r["in_footprint_km2"] for r in rows) or 1
        for r in rows:
            r["pct_of_footprint"] = round(100 * r["in_footprint_km2"] / tot, 1)
    return rows


def landcover(final, tr, crs):
    import odc.stac
    from odc.geo.geobox import GeoBox
    cat = open_catalog()
    items = list(cat.search(collections=["esa-worldcover"], bbox=CFG["bbox"]).items())
    newest = max(i.properties.get("start_datetime", "")[:4] for i in items)
    items = [i for i in items if i.properties.get("start_datetime", "")[:4] == newest]
    gb = GeoBox(final.shape, tr, str(crs))
    lc = odc.stac.load(items, bands=["map"], geobox=gb, groupby="solar_day",
                       resampling="mode")["map"].isel(time=0).values
    return lc, newest


def step5(overlays: dict):
    header(5, "IMPACT")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch
    o = out_dir()
    idx, tr, crs = read_indices(["dNBR", "NDVI_pre", "NDVI_post"])
    final = read_tif(need(o / "footprint_mask.tif", 4))[0] == 1
    px = CFG["resolution"] ** 2 / 1e6
    result = {}

    sev = classify(idx["dNBR"], SEVERITY_BINS)
    write_tif(o / "severity_dNBR.tif", sev, tr, crs, 255, ["USGS dNBR class 0-6"])
    result["severity"] = class_table(sev, SEVERITY_LABELS, px, final)
    say("Severity inside the footprint:")
    for r in result["severity"][3:]:
        say(f"   {r['label']:24s} {r['in_footprint_km2']:8.1f} km²  {r['pct_of_footprint']:5.1f} %")

    n0, n1 = classify(idx["NDVI_pre"], NDVI_BINS), classify(idx["NDVI_post"], NDVI_BINS)
    t0, t1 = class_table(n0, NDVI_LABELS, px), class_table(n1, NDVI_LABELS, px)
    result["ndvi_classes"] = [{"label": a["label"], "pre_km2": a["aoi_km2"], "post_km2": b["aoi_km2"],
                               "change_km2": round(b["aoi_km2"] - a["aoi_km2"], 3)}
                              for a, b in zip(t0, t1)]
    say("NDVI zoning (Chuvieco), whole AOI, pre -> post:")
    for r in result["ndvi_classes"]:
        say(f"   {r['label']:32s} {r['pre_km2']:7.1f} -> {r['post_km2']:7.1f} km²  "
            f"({r['change_km2']:+.1f})")

    try:
        lc, year = landcover(final, tr, crs)
        rows = []
        for code, name in WORLDCOVER.items():
            n_aoi = int((lc == code).sum())
            if n_aoi:
                n_b = int(((lc == code) & final).sum())
                rows.append({"class": name, "aoi_km2": round(n_aoi * px, 2),
                             "burned_km2": round(n_b * px, 2),
                             "pct_of_class_burned": round(100 * n_b / n_aoi, 1)})
        tot = sum(r["burned_km2"] for r in rows) or 1
        for r in rows:
            r["pct_of_footprint"] = round(100 * r["burned_km2"] / tot, 1)
        rows.sort(key=lambda r: -r["burned_km2"])
        result["landcover"] = {"source": f"ESA WorldCover {year}", "classes": rows}
        say(f"Land cover burned (ESA WorldCover {year}):")
        for r in rows:
            say(f"   {r['class']:24s} {r['burned_km2']:8.2f} km²  = {r['pct_of_footprint']:5.1f} % of "
                f"footprint, {r['pct_of_class_burned']:5.1f} % of this class in the AOI")
    except Exception as e:
        say(f"Land-cover step skipped: {e}")

    if overlays:
        import geopandas as gpd
        fp = gpd.read_file(o / "burn_footprint.gpkg", layer="footprint_dissolved")
        fpu = fp.union_all() if hasattr(fp, "union_all") else fp.unary_union
        result["overlays"] = {}
        for name, path in overlays.items():
            lyr = gpd.read_file(path).to_crs(fp.crs)
            hit = lyr[lyr.intersects(fpu)].copy()
            r = {"features_total": len(lyr), "features_affected": len(hit)}
            if len(hit) and hit.geom_type.isin(["Polygon", "MultiPolygon"]).all():
                hit["burned_km2"] = hit.geometry.intersection(fpu).area / 1e6
                r["burned_km2"] = round(float(hit["burned_km2"].sum()), 3)
                hit.drop(columns="geometry").to_csv(o / f"overlay_{name}.csv", index=False)
            result["overlays"][name] = r
            say(f"Overlay {name}: {len(hit)} of {len(lyr)} features touched by the fire")

    fig, ax = plt.subplots(1, 2, figsize=(16, 7), constrained_layout=True)
    ax[0].imshow(np.ma.masked_equal(sev, 255), cmap=ListedColormap(SEVERITY_COLORS),
                 norm=BoundaryNorm(np.arange(-.5, 7.5), 7), interpolation="nearest")
    ax[0].contour(final, levels=[0.5], colors="black", linewidths=0.4)
    ax[0].legend(handles=[Patch(color=c, label=l) for c, l in zip(SEVERITY_COLORS, SEVERITY_LABELS)],
                 loc="lower right", fontsize=7)
    ax[0].set_title("Burn severity (dNBR, USGS)"); ax[0].set_axis_off()
    sv = result["severity"][3:]
    ax[1].barh([r["label"] for r in sv], [r["in_footprint_km2"] for r in sv],
               color=SEVERITY_COLORS[3:])
    for i, r in enumerate(sv):
        ax[1].text(r["in_footprint_km2"], i, f" {r['in_footprint_km2']:.1f} km² ({r['pct_of_footprint']}%)",
                   va="center", fontsize=9)
    ax[1].set_xlabel("km² inside the footprint"); ax[1].set_title("Area per severity class")
    fig.savefig(o / "step5_severity.png", dpi=130)
    plt.close(fig)
    (o / "step5_impact.json").write_text(json.dumps(result, indent=2))
    say("Saved -> severity_dNBR.tif, step5_severity.png, step5_impact.json")


# =============================================================================
#  STEP 6 - FIRE WEATHER vs CLIMATE (ERA5 via Open-Meteo)
# =============================================================================
#  Why: to test the "why did it get out of control" question with data.
#  We ask Open-Meteo for daily values at the AOI centre:
#    - Oct 2025 -> Aug 2026: max temperature, min humidity, wind gusts, rain
#    - 1991-2020: the climate normal to compare against
#  It shows the fuel story (wet winter -> grass growth, dry summer -> curing)
#  and the weather on the ignition days (heat, very low humidity, wind).
# =============================================================================
def om(lat, lon, start, end, daily):
    import pandas as pd
    import requests
    r = requests.get(OPEN_METEO, timeout=120, params={
        "latitude": lat, "longitude": lon, "start_date": start, "end_date": end,
        "daily": ",".join(daily), "timezone": "Europe/Madrid"})
    r.raise_for_status()
    df = pd.DataFrame(r.json()["daily"])
    df["time"] = pd.to_datetime(df["time"])
    return df.set_index("time")


def step6():
    header(6, "FIRE WEATHER vs CLIMATE")
    import matplotlib.pyplot as plt
    import pandas as pd
    lat, lon = (CFG["bbox"][1] + CFG["bbox"][3]) / 2, (CFG["bbox"][0] + CFG["bbox"][2]) / 2
    fire = pd.Timestamp(CFG["fire_start"])
    obs = om(lat, lon, f"{fire.year - 1}-10-01", (fire + pd.Timedelta(days=14)).date().isoformat(),
             ["temperature_2m_max", "relative_humidity_2m_min", "wind_gusts_10m_max",
              "precipitation_sum"])
    clim = om(lat, lon, "1991-01-01", "2020-12-31", ["temperature_2m_max", "precipitation_sum"])
    obs["tmax_normal"] = clim["temperature_2m_max"].groupby(clim.index.day_of_year).mean() \
        .reindex(obs.index.day_of_year).values
    pm = obs["precipitation_sum"].resample("MS").sum()
    pc = clim["precipitation_sum"].resample("MS").sum()
    pc = pc.groupby(pc.index.month).mean()
    rain = pd.DataFrame({"obs_mm": pm, "normal_mm": pc.reindex(pm.index.month).values})
    rain["pct_of_normal"] = (100 * rain["obs_mm"] / rain["normal_mm"]).round(0)
    rain.to_csv(out_dir() / "step6_rain_vs_normal.csv")

    wk = obs.loc[fire - pd.Timedelta(days=3): fire + pd.Timedelta(days=4)]
    hot = obs[(obs.index >= fire - pd.Timedelta(days=45)) & (obs.index <= fire)]
    res = {"point": [round(lat, 3), round(lon, 3)],
           "ignition_week": {"tmax_max_C": round(float(wk["temperature_2m_max"].max()), 1),
                             "rh_min_lowest_pct": round(float(wk["relative_humidity_2m_min"].min()), 0),
                             "gust_max_kmh": round(float(wk["wind_gusts_10m_max"].max()), 0)},
           "days_above_normal_by_5C_in_45_days_before": int(
               (hot["temperature_2m_max"] - hot["tmax_normal"] > 5).sum()),
           "rain_pct_of_normal": {d.strftime("%b %Y"): float(v)
                                  for d, v in rain["pct_of_normal"].items()}}
    (out_dir() / "step6_weather.json").write_text(json.dumps(res, indent=2))
    say(f"Ignition week: Tmax {res['ignition_week']['tmax_max_C']} °C, humidity down to "
        f"{res['ignition_week']['rh_min_lowest_pct']} %, gusts {res['ignition_week']['gust_max_kmh']} km/h")
    say(f"Days >5 °C above normal in the 45 days before: "
        f"{res['days_above_normal_by_5C_in_45_days_before']}")
    say("Rain as % of normal: " + ", ".join(f"{k} {v:.0f}%" for k, v in res["rain_pct_of_normal"].items()))

    win = obs.loc[fire - pd.Timedelta(days=40): fire + pd.Timedelta(days=14)]
    fig, ax = plt.subplots(3, 1, figsize=(10, 10), constrained_layout=True)
    ax[0].plot(win.index, win["temperature_2m_max"], c="#c0392b", label="Tmax 2026")
    ax[0].plot(win.index, win["tmax_normal"], c="grey", ls="--", label="Normal 1991-2020")
    ax[0].fill_between(win.index, win["tmax_normal"], win["temperature_2m_max"],
                       where=win["temperature_2m_max"] > win["tmax_normal"], color="#e74c3c", alpha=.2)
    ax[0].set_ylabel("°C"); ax[0].legend()
    ax[1].plot(win.index, win["relative_humidity_2m_min"], c="#2980b9", label="Min humidity %")
    a2 = ax[1].twinx()
    a2.plot(win.index, win["wind_gusts_10m_max"], c="#8e44ad", label="Max gust km/h")
    ax[1].legend(loc="upper left"); a2.legend(loc="upper right")
    for a in ax[:2]:
        a.axvline(fire, c="k", lw=1)
    x = np.arange(len(rain))
    ax[2].bar(x - .2, rain["normal_mm"], .4, color="#95a5a6", label="Normal 1991-2020")
    ax[2].bar(x + .2, rain["obs_mm"], .4, color="#16a085", label="2025-2026")
    ax[2].set_xticks(x, [d.strftime("%b %y") for d in rain.index], rotation=45)
    ax[2].set_ylabel("mm"); ax[2].legend()
    ax[0].set_title(f"Fire weather at the AOI centre ({lat:.2f}N, {lon:.2f}E), ERA5")
    ax[2].set_title("Monthly rain: wet winter (fuel grows) -> dry summer (fuel dries)")
    fig.savefig(out_dir() / "step6_fire_weather.png", dpi=140)
    plt.close(fig)
    say("Saved -> step6_fire_weather.png, step6_weather.json, step6_rain_vs_normal.csv")


# =============================================================================
def main():
    p = argparse.ArgumentParser(description="Sentinel-2 burn footprint, step by step")
    p.add_argument("--step", required=True, choices=["1", "2", "3", "4", "5", "6", "all"])
    p.add_argument("--overlay", action="append", default=[], metavar="NAME=PATH",
                   help="extra vector layer for step 5, e.g. municipios=limits.shp")
    p.add_argument("--out", default=CFG["out"])
    a = p.parse_args()
    CFG["out"] = a.out
    overlays = dict(o.split("=", 1) for o in a.overlay)
    steps = {"1": step1, "2": step2, "3": step3, "4": step4,
             "5": lambda: step5(overlays), "6": step6}
    for s in (["1", "2", "3", "4", "5", "6"] if a.step == "all" else [a.step]):
        steps[s]()
    say("\nDone.")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    main()
