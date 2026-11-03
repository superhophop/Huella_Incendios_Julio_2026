#!/usr/bin/env python3
"""
burn_footprint.py
=================
Automated pre/post-fire burn-footprint mapping with Landsat 8/9 (default) or Sentinel-2.

Case study: July 2026 Ávila–Madrid–Toledo fire complex, Sierra Oeste de Madrid AOI
(San Martín de Valdeiglesias, Pelayos de la Presa, Navas del Rey, Chapinería,
Villa del Prado, Aldea del Fresno).

The workflow automates the change-detection methodology of the course practical
"Aplicación 3 – Detección de cambios: incendio forestal en la Sierra de Gata"
(Key & Benson 2006; Chuvieco 2015; Bastarrika et al. 2011; Martín 1998;
Pinty & Verstraete 1992). What was done by hand in ArcGIS Pro is done here in code:

  Course step                         ->  This script
  ---------------------------------------------------------------------------
  3. Data acquisition (EarthExplorer) ->  STAC search (Planetary Computer), cloud
                                          cover measured INSIDE the AOI, picks the
                                          pre scene ~28 d before ignition and the
                                          first clear post scene after extinction
  4. DN -> reflectance (Raster Calc.) ->  Landsat C2 L2SP: DN*0.0000275 - 0.2,
                                          negatives -> 0 (Con(x<0, 0, x))
  5.1 5-4-3 + NIR change composite    ->  01_false_color_543.png, 02_nir_change.png
  5.2 NDVI + Chuvieco reclass, ΔNDVI  ->  03_ndvi_classes.png, ndvi_class_areas.csv
      NBR, BAI, ΔBAI, GEMI, ΔGEMI     ->  04_indices.png + GeoTIFFs
  6. Segmentation (Bastarrika 2011)   ->  ΔNDVI<-0.17767, NBRpost<-0.15006,
                                          ΔBAI>144.3835, ΔGEMI<-0.11578
     Validation: M separability index ->  separability_M.csv
     Majority filter (8 nbrs, majority) ->  applied to every segmentation
     Raster to polygon + area         ->  index_areas.csv, per-index GPKG layers
  7. Merge + manual clean + Dissolve  ->  union of all indices, polygons kept only
                                          where >= --min-votes indices agree
                                          (automates the manual clean-up), dissolve
     Layer To KML (Google Earth)      ->  burn_footprint.kml
     Compare with INFOEX / Copernicus ->  --reference perimeter -> IoU, omission,
                                          commission
  Extra: dNBR severity (USGS), impact by land cover / municipality / buildings,
         fire-weather context vs 1991-2020 climatology.

Install (conda recommended):
    conda create -n burn -c conda-forge python=3.11 pystac-client planetary-computer \
        odc-stac rioxarray rasterio geopandas scipy matplotlib pandas requests

Usage:
    python burn_footprint.py                         # Landsat 8/9, 2026 defaults
    python burn_footprint.py --sensor sentinel2      # 20 m Sentinel-2 instead
    python burn_footprint.py --method dNDVI          # one index instead of merge
    python burn_footprint.py --reference effis_perimeter.shp
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from dataclasses import asdict, dataclass, field
from datetime import date, timedelta
from pathlib import Path

import numpy as np

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
OPEN_METEO_ARCHIVE = "https://archive-api.open-meteo.com/v1/archive"

# --------------------------------------------------------------------------- #
# Sensors - every band is referred to by its role (red, nir...) from here on
# --------------------------------------------------------------------------- #
SENSORS = {
    "landsat": {
        "collection": "landsat-c2-l2",
        "resolution": 30,
        "crs": "EPSG:32630",  # WGS84 / UTM 30N - native Landsat grid over Madrid
        "bands": {"blue": "blue", "green": "green", "red": "red",   # B2 B3 B4
                  "nir": "nir08", "swir2": "swir22", "qa": "qa_pixel"},  # B5 B7
        "query": {"platform": {"in": ["landsat-8", "landsat-9"]}},  # no L7 (SLC-off)
        "band_numbers": "NIR=B5, Red=B4, Green=B3, SWIR2=B7",
    },
    "sentinel2": {
        "collection": "sentinel-2-l2a",
        "resolution": 10,
        "crs": "EPSG:32630",
        # nir10 = B08 (10 m) for the 8-4-3 composite and the red/NIR indices
        # (NDVI, BAI, GEMI); nir = B8A (20 m, narrow) pairs with B12 for NBR/dNBR.
        "bands": {"blue": "B02", "green": "B03", "red": "B04", "nir10": "B08",
                  "nir": "B8A", "swir2": "B12", "qa": "SCL"},
        "query": {},
        "band_numbers": "NIR=B8 (10 m), Red=B4, Green=B3 - S2 equivalent of Landsat 5-4-3",
    },
}


@dataclass
class Config:
    sensor: str = "landsat"
    # AOI (lon_min, lat_min, lon_max, lat_max) - Sierra Oeste de Madrid, where the
    # Almorox (Toledo) and San Martín de Valdeiglesias fronts merged 22-26 Jul 2026.
    bbox: list = field(default_factory=lambda: [-4.46, 40.26, -4.17, 40.46])
    fire_start: str = "2026-07-22"
    days_before: int = 28
    pre_window_days: int = 7
    pre_max_expand_days: int = 28
    post_start: str = "2026-07-30"     # Key & Benson initial assessment: fire out
    post_end: str = "2026-09-15"
    max_scene_cloud: float = 40.0      # tile-level pre-filter (%)
    max_aoi_cloud: float = 1.0         # % cloudy pixels allowed inside AOI
    min_aoi_coverage: float = 99.0
    crs: str | None = None             # default: sensor's
    resolution: int | None = None      # default: sensor's
    method: str = "merge"              # merge | dNDVI | NBR_post | dBAI | dGEMI | dNBR
    thresholds: str = "scene"          # scene (calibrated on this image) | bastarrika
    k_sigma: float = 3.0               # scene thresholds = unburned median ± k·MAD
    pixel_votes: int = 2               # merge: indices that must flag a pixel
    min_votes: int = 3                 # merge: indices that must agree in a patch core
    dnbr_threshold: float = 0.10       # for --method dNBR with bastarrika thresholds
    majority_neighbors: int = 8
    majority_threshold: str = "majority"   # majority | half
    min_patch_ha: float = 1.0
    pre_date: str | None = None
    post_date: str | None = None
    reference: str | None = None       # official perimeter for validation
    out: str = "outputs"
    weather: bool = True
    landcover: bool = True
    overlays: dict = field(default_factory=dict)

    def __post_init__(self):
        s = SENSORS[self.sensor]
        self.crs = self.crs or s["crs"]
        self.resolution = self.resolution or s["resolution"]


# Bastarrika et al. (2011) thresholds for Mediterranean burned areas (course, sec. 6)
BASTARRIKA = {
    "dNDVI":    ("<", -0.17767),
    "NBR_post": ("<", -0.15006),
    "dBAI":     (">", 144.3835),
    "dGEMI":    ("<", -0.11578),
}

# Direction of each index when calibrating thresholds on the scene itself
INDEX_DIRECTION = {"dNDVI": "<", "dGEMI": "<", "dBAI": ">", "dNBR": ">", "NBR_post": "<"}
# Indices voted on in the merge for each threshold set
MERGE_INDICES = {"scene": ("dNDVI", "dGEMI", "dBAI", "dNBR"),
                 "bastarrika": ("dNDVI", "NBR_post", "dBAI", "dGEMI")}
CHAR_INDEX = "dBAI"   # char check: distance to the charcoal spectral point (Martín 1998)

# Chuvieco (2015) NDVI zoning (course, sec. 5.2-A)
NDVI_BINS = [0.0, 0.2, 0.4, 0.6]
NDVI_LABELS = ["Water / shadow (<0)", "Bare soil (0-0.2)", "Sparse veg + soil (0.2-0.4)",
               "Dispersed vegetation (0.4-0.6)", "Dense vegetation (>0.6)"]
NDVI_COLORS = ["#000000", "#fffbd5", "#d4f5c4", "#4fd22a", "#1f6b12"]

# USGS / Key & Benson (2006) dNBR severity
SEVERITY_BINS = [-0.25, -0.10, 0.10, 0.27, 0.44, 0.66]
SEVERITY_LABELS = ["Enhanced regrowth, high", "Enhanced regrowth, low", "Unburned",
                   "Low severity", "Moderate-low severity", "Moderate-high severity",
                   "High severity"]
SEVERITY_COLORS = ["#1a9850", "#91cf60", "#d9d9d9", "#fee08b", "#fdae61", "#f46d43", "#a50026"]

WORLDCOVER = {10: "Tree cover", 20: "Shrubland", 30: "Grassland", 40: "Cropland",
              50: "Built-up", 60: "Bare / sparse vegetation", 70: "Snow and ice",
              80: "Permanent water", 90: "Herbaceous wetland", 95: "Mangroves",
              100: "Moss and lichen"}


def log(msg: str) -> None:
    print(f"[burn] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# Quality masks
# --------------------------------------------------------------------------- #
def qa_masks(sensor: str, qa: np.ndarray) -> dict:
    """Return boolean masks: nodata, cloud (incl. shadow/cirrus), snow, water."""
    if sensor == "landsat":
        q = qa.astype(np.uint32)
        bit = lambda b: ((q >> b) & 1) == 1  # noqa: E731
        return {"nodata": bit(0) | (q == 0),
                "cloud": bit(1) | bit(2) | bit(3) | bit(4),  # dilated, cirrus, cloud, shadow
                "snow": bit(5), "water": bit(7)}
    # Sentinel-2 SCL. Class 3 (cloud shadow) is NOT masked: Sen2Cor often labels
    # fresh burn scars as shadow, which would punch holes in the scar.
    return {"nodata": np.isin(qa, [0, 1]), "cloud": np.isin(qa, [8, 9, 10]),
            "snow": qa == 11, "water": qa == 6}


# --------------------------------------------------------------------------- #
# 1. Scene search (cloud cover measured inside the AOI)
# --------------------------------------------------------------------------- #
def open_catalog():
    import planetary_computer
    import pystac_client
    return pystac_client.Client.open(STAC_URL, modifier=planetary_computer.sign_inplace)


def search_items_by_day(catalog, cfg: Config, start: str, end: str, max_cloud: float) -> dict:
    s = SENSORS[cfg.sensor]
    query = {"eo:cloud_cover": {"lt": max_cloud}, **s["query"]}
    search = catalog.search(collections=[s["collection"]], bbox=cfg.bbox,
                            datetime=f"{start}/{end}", query=query)
    by_day: dict = {}
    for it in search.items():
        by_day.setdefault(it.datetime.date().isoformat(), []).append(it)
    return dict(sorted(by_day.items()))


def score_day(items, cfg: Config) -> dict:
    import odc.stac
    qa_name = SENSORS[cfg.sensor]["bands"]["qa"]
    ds = odc.stac.load(items, bands=[qa_name], bbox=cfg.bbox, crs=cfg.crs,
                       resolution=cfg.resolution * 3, groupby="solar_day",
                       resampling="nearest")
    m = qa_masks(cfg.sensor, ds[qa_name].isel(time=0).values)
    valid = ~m["nodata"]
    coverage = float(valid.mean() * 100)
    cloud = float(m["cloud"][valid].mean() * 100) if valid.any() else 100.0
    return {"coverage_pct": round(coverage, 2), "aoi_cloud_pct": round(cloud, 3)}


def rank_days(by_day: dict, cfg: Config) -> list[dict]:
    rows = []
    for d, items in by_day.items():
        try:
            s = score_day(items, cfg)
        except Exception as e:
            log(f"  {d}: could not score ({e})")
            continue
        s.update(date=d, n_items=len(items),
                 platforms=sorted({i.properties.get("platform", "?") for i in items}))
        log(f"  {d} {s['platforms']}: AOI cloud {s['aoi_cloud_pct']:.2f} %, "
            f"coverage {s['coverage_pct']:.1f} %")
        rows.append(s)
    return rows


def is_clear(row, cfg: Config) -> bool:
    return row["aoi_cloud_pct"] <= cfg.max_aoi_cloud and row["coverage_pct"] >= cfg.min_aoi_coverage


def forced(catalog, cfg: Config, d: str):
    by_day = search_items_by_day(catalog, cfg, d, d, 100)
    if d not in by_day:
        sys.exit(f"No {cfg.sensor} data over the AOI on {d}")
    return d, by_day[d], score_day(by_day[d], cfg)


def pick_pre(catalog, cfg: Config):
    if cfg.pre_date:
        return forced(catalog, cfg, cfg.pre_date)
    fire = date.fromisoformat(cfg.fire_start)
    target = fire - timedelta(days=cfg.days_before)
    half = cfg.pre_window_days
    while half <= cfg.pre_max_expand_days:
        start = target - timedelta(days=half)
        end = min(target + timedelta(days=half), fire - timedelta(days=1))
        log(f"PRE search {start} -> {end} (target {target} = ignition - {cfg.days_before} d)")
        by_day = search_items_by_day(catalog, cfg, start.isoformat(), end.isoformat(),
                                     cfg.max_scene_cloud)
        clear = [r for r in rank_days(by_day, cfg) if is_clear(r, cfg)]
        if clear:
            clear.sort(key=lambda r: (abs((date.fromisoformat(r["date"]) - target).days),
                                      r["aoi_cloud_pct"]))
            return clear[0]["date"], by_day[clear[0]["date"]], clear[0]
        half += 7
        log(f"  no clear pre-fire scene; widening window to +/-{half} days")
    sys.exit("No cloud-free pre-fire scene. Relax --max-aoi-cloud or set --pre-date.")


def pick_post(catalog, cfg: Config):
    if cfg.post_date:
        return forced(catalog, cfg, cfg.post_date)
    log(f"POST search {cfg.post_start} -> {cfg.post_end}")
    by_day = search_items_by_day(catalog, cfg, cfg.post_start, cfg.post_end, cfg.max_scene_cloud)
    clear = [r for r in rank_days(by_day, cfg) if is_clear(r, cfg)]
    if not clear:
        sys.exit("No cloud-free post-fire scene. Extend --post-end or relax --max-aoi-cloud.")
    clear.sort(key=lambda r: r["date"])  # earliest clear = least regrowth
    return clear[0]["date"], by_day[clear[0]["date"]], clear[0]


# --------------------------------------------------------------------------- #
# 2. Loading + DN -> reflectance
# --------------------------------------------------------------------------- #
def s2_offset(items) -> int:
    """Sentinel-2 baseline >= 04.00 stores DN = 10000*rho + 1000."""
    vals = []
    for it in items:
        try:
            vals.append(float(it.properties.get("s2:processing_baseline")))
        except (TypeError, ValueError):
            pass
    return -1000 if (not vals or max(vals) >= 4.0) else 0


def to_reflectance(dn: np.ndarray, sensor: str, s2_off: int = -1000) -> np.ndarray:
    """Course step 4: ("SR_Bx" * 0.0000275) + (-0.2), then Con(x < 0, 0, x)."""
    dn = dn.astype("float32")
    if sensor == "landsat":
        r = dn * 0.0000275 - 0.2
    else:
        r = (dn + s2_off) / 10000.0
    r[dn == 0] = np.nan                  # fill value
    return np.where(r < 0, 0, r).astype("float32")


def load_scene(items, cfg: Config, geobox=None):
    import odc.stac
    names = SENSORS[cfg.sensor]["bands"]
    kw = dict(bands=list(names.values()), groupby="solar_day",
              resampling={names["qa"]: "nearest", "*": "bilinear"})
    if geobox is None:
        kw.update(bbox=cfg.bbox, crs=cfg.crs, resolution=cfg.resolution)
    else:
        kw.update(geobox=geobox)          # post scene snaps onto the pre grid
    ds = odc.stac.load(items, **kw).isel(time=0)
    off = s2_offset(items) if cfg.sensor == "sentinel2" else 0
    bands = {role: to_reflectance(ds[asset].values, cfg.sensor, off)
             for role, asset in names.items() if role != "qa"}
    bands["qa"] = ds[names["qa"]].values
    return ds, bands


# --------------------------------------------------------------------------- #
# 3. Spectral indices (course sec. 5.2) - pure numpy, unit-testable
# --------------------------------------------------------------------------- #
def _nd(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        den = a + b
        return np.where(den > 0, (a - b) / den, np.nan).astype("float32")


def ndvi(nir, red):
    """Rouse et al. (1973): (NIR - Red) / (NIR + Red)."""
    return _nd(nir, red)


def nbr(nir, swir2):
    """Key & Benson (1996, 2006): (NIR - SWIR2) / (NIR + SWIR2)."""
    return _nd(nir, swir2)


def bai(nir, red):
    """Martín (1998): 1 / ((0.1 - Red)^2 + (0.06 - NIR)^2)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return (1.0 / ((0.1 - red) ** 2 + (0.06 - nir) ** 2)).astype("float32")


def gemi(nir, red):
    """Pinty & Verstraete (1992): eta*(1 - 0.25*eta) - (Red - 0.125)/(1 - Red)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = (2 * (nir ** 2 - red ** 2) + 1.5 * nir + 0.5 * red) / (nir + red + 0.5)
        out = eta * (1 - 0.25 * eta) - (red - 0.125) / (1 - red)
        return np.where(red < 1, out, np.nan).astype("float32")


def compute_indices(pre: dict, post: dict) -> dict:
    """Two-date indices with the course's sign conventions (Δ = post - pre),
    plus dNBR (pre - post) for USGS severity."""
    # Red/NIR indices use the sharpest NIR available (Sentinel-2 B08 at 10 m);
    # NBR always uses the NIR that pairs with SWIR2 (Landsat B5 / Sentinel-2 B8A).
    rn0, rn1 = pre.get("nir10", pre["nir"]), post.get("nir10", post["nir"])
    n0, n1 = ndvi(rn0, pre["red"]), ndvi(rn1, post["red"])
    b0, b1 = nbr(pre["nir"], pre["swir2"]), nbr(post["nir"], post["swir2"])
    a0, a1 = bai(rn0, pre["red"]), bai(rn1, post["red"])
    g0, g1 = gemi(rn0, pre["red"]), gemi(rn1, post["red"])
    with np.errstate(divide="ignore", invalid="ignore"):
        rdnbr = (b0 - b1) / np.sqrt(np.abs(b0))
    return {"NDVI_pre": n0, "NDVI_post": n1, "dNDVI": n1 - n0,
            "NBR_pre": b0, "NBR_post": b1, "dNBR": b0 - b1, "RdNBR": rdnbr.astype("float32"),
            "BAI_pre": a0, "BAI_post": a1, "dBAI": a1 - a0,
            "GEMI_pre": g0, "GEMI_post": g1, "dGEMI": g1 - g0}


def classify(values: np.ndarray, bins) -> np.ndarray:
    cls = np.digitize(values, bins).astype("uint8")
    cls[~np.isfinite(values)] = 255
    return cls


def threshold(values: np.ndarray, op: str, thr: float) -> np.ndarray:
    with np.errstate(invalid="ignore"):
        return (values < thr) if op == "<" else (values > thr)


# --------------------------------------------------------------------------- #
# 4. Segmentation post-processing (course sec. 6-7)
# --------------------------------------------------------------------------- #
def majority_filter(mask: np.ndarray, neighbors: int = 8, mode: str = "majority") -> np.ndarray:
    """ArcGIS 'Majority Filter' equivalent. A cell is replaced when MAJORITY
    (5 of 8 / 3 of 4) or HALF (4 of 8 / 2 of 4) of its neighbours share a value."""
    from scipy import ndimage
    if neighbors == 8:
        k = np.ones((3, 3), int); k[1, 1] = 0
    else:
        k = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    need = {(8, "majority"): 5, (8, "half"): 4, (4, "majority"): 3, (4, "half"): 2}[(neighbors, mode)]
    n_burn = ndimage.convolve(mask.astype(int), k, mode="nearest")
    n_unb = neighbors - n_burn
    to_burn = (n_burn >= need) & ~(n_unb >= need)
    to_unb = (n_unb >= need) & ~(n_burn >= need)
    return np.where(to_burn, True, np.where(to_unb, False, mask))


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


def clean(mask: np.ndarray, valid: np.ndarray, min_px: int, nb=8, mode="majority") -> np.ndarray:
    m = majority_filter(mask & valid, nb, mode)
    m = remove_small(m, min_px)             # isolated specks
    m = ~remove_small(~m, min_px)           # small unburned pinholes
    return m & valid


def scene_thresholds(idx: dict, valid: np.ndarray, k: float = 3.0) -> dict:
    """Thresholds calibrated on THIS image pair. Unburned reference = pixels in the
    USGS 'unburned' dNBR class (-0.1 to 0.1); threshold = median ± k·(1.4826·MAD)
    of that reference. Robust statistics are used because a few extreme pixels
    (terrain shadow when sun elevation differs between dates, reservoir edges)
    inflate an ordinary standard deviation - on Sentinel-2 10 m that pushed the
    ΔBAI threshold from ~16 to 192. Bastarrika's fixed values (calibrated on a
    dense shrub fire) miss low-severity grass and dehesa burns."""
    d = idx["dNBR"]
    with np.errstate(invalid="ignore"):
        ref = valid & np.isfinite(d) & (d > -0.1) & (d < 0.1)
    out = {}
    for name, op in INDEX_DIRECTION.items():
        v = idx[name][ref]
        v = v[np.isfinite(v)]
        med = float(np.median(v))
        s = 1.4826 * float(np.median(np.abs(v - med)))   # robust sd
        out[name] = (op, round(med - k * s if op == "<" else med + k * s, 5))
    return out


def merge_and_clean(seg: dict, valid: np.ndarray, min_votes: int, min_px: int,
                    nb=8, mode="majority", pixel_votes: int = 2, char_key: str = CHAR_INDEX,
                    strong_char: np.ndarray | None = None):
    """Course sec. 7: Merge the index segmentations, then drop polygons outside the
    real burn. The manual clean-up is automated:
      - pixel level: a pixel is burned when >= `pixel_votes` indices flag it
        (this keeps low-severity grass / dehesa burns);
      - patch level: a patch survives only if it has >= min_px 'core' pixels where
        >= `min_votes` indices agree AND there is strong char (`strong_char`,
        default: the char index, ΔBAI). Low-severity burns and freshly harvested
        cereal look alike spectrally, but a real fire patch always contains some
        truly charred ground, a harvested field does not."""
    from scipy import ndimage
    votes = np.sum([m.astype(np.uint8) for m in seg.values()], axis=0)
    if strong_char is not None:
        char = strong_char
    else:
        char = seg[char_key] if char_key in seg else np.ones_like(valid)
    union = majority_filter((votes >= pixel_votes) & valid, nb, mode)
    lab, n = _label(union)
    if n == 0:
        return union, votes
    idx = np.arange(1, n + 1)
    core = ndimage.sum((votes >= min_votes) & char & valid, lab, idx)
    size = ndimage.sum(np.ones_like(union), lab, idx)
    keep = np.zeros(n + 1, bool)
    keep[1:] = (core >= min_px) & (size >= min_px)
    final = keep[lab]
    final = ~remove_small(~final, min_px)   # Dissolve-like hole filling
    return final & valid, votes


def separability_m(values: np.ndarray, burned: np.ndarray, valid: np.ndarray):
    """Course sec. 6 'Validación': M = |mu_b - mu_u| / (sd_b + sd_u). M > 1 = good."""
    ok = valid & np.isfinite(values)
    b, u = values[burned & ok], values[~burned & ok]
    if b.size < 2 or u.size < 2:
        return None
    return {"mean_burned": float(b.mean()), "sd_burned": float(b.std()),
            "mean_unburned": float(u.mean()), "sd_unburned": float(u.std()),
            "M": float(abs(b.mean() - u.mean()) / (b.std() + u.std()))}


def area_stats(mask: np.ndarray, res_m: float) -> dict:
    n = int(mask.sum())
    return {"pixels": n, "area_ha": round(n * res_m ** 2 / 1e4, 2),
            "area_km2": round(n * res_m ** 2 / 1e6, 3)}


def class_area_table(cls: np.ndarray, labels, res_m: float, mask=None):
    import pandas as pd
    px = res_m ** 2 / 1e6
    rows = []
    for i, lab in enumerate(labels):
        sel = cls == i
        row = {"class": i, "label": lab, "aoi_km2": round(int(sel.sum()) * px, 3)}
        if mask is not None:
            row["in_footprint_km2"] = round(int((sel & mask).sum()) * px, 3)
        rows.append(row)
    df = pd.DataFrame(rows)
    if mask is not None:
        tot = df["in_footprint_km2"].sum()
        df["pct_of_footprint"] = (100 * df["in_footprint_km2"] / tot).round(1) if tot else 0.0
    return df


# --------------------------------------------------------------------------- #
# 5. Vectors: Raster to Polygon, Dissolve, KML, reference comparison
# --------------------------------------------------------------------------- #
def geobox_transform(ds):
    gb = ds.odc.geobox
    return getattr(gb, "transform", None) or gb.affine


def vectorise(mask: np.ndarray, transform, crs: str, dissolve=True):
    import geopandas as gpd
    from rasterio import features
    from shapely.geometry import shape
    geoms = [shape(g) for g, v in features.shapes(mask.astype("uint8"), mask=mask,
                                                   transform=transform) if v == 1]
    if not geoms:
        return gpd.GeoDataFrame({"area_ha": []}, geometry=[], crs=crs)
    gdf = gpd.GeoDataFrame(geometry=geoms, crs=crs)
    if dissolve:
        gdf = gdf.dissolve()
    gdf = gdf.reset_index(drop=True)
    gdf["Shape_Area_m2"] = gdf.area.round(1)
    gdf["Shape_Length_m"] = gdf.length.round(1)
    gdf["area_ha"] = (gdf.area / 1e4).round(2)
    gdf["area_km2"] = (gdf.area / 1e6).round(4)
    return gdf


def outer_perimeter_km(gdf) -> float:
    """Perimeter of the footprint's outer edges only (interior unburned holes
    excluded) - comparable to official fire-perimeter figures."""
    from shapely.geometry import Polygon
    if not len(gdf):
        return 0.0
    polys = gdf.explode(index_parts=False).geometry
    return round(sum(Polygon(p.exterior).length for p in polys) / 1e3, 2)


def write_kml(gdf, path: Path, name: str):
    """Layer To KML equivalent - opens directly in Google Earth."""
    from shapely.geometry import MultiPolygon, Polygon
    g = gdf.to_crs(4326)
    ring = lambda cs: " ".join(f"{x:.6f},{y:.6f},0" for x, y, *_ in cs)  # noqa: E731

    def poly_kml(p: Polygon):
        inner = "".join(f"<innerBoundaryIs><LinearRing><coordinates>{ring(i.coords)}"
                        f"</coordinates></LinearRing></innerBoundaryIs>" for i in p.interiors)
        return (f"<Polygon><outerBoundaryIs><LinearRing><coordinates>{ring(p.exterior.coords)}"
                f"</coordinates></LinearRing></outerBoundaryIs>{inner}</Polygon>")

    marks = []
    for _, r in g.iterrows():
        geom = r.geometry
        polys = [geom] if isinstance(geom, Polygon) else list(getattr(geom, "geoms", []))
        body = "".join(poly_kml(p) for p in polys if isinstance(p, Polygon))
        if isinstance(geom, MultiPolygon) or len(polys) > 1:
            body = f"<MultiGeometry>{body}</MultiGeometry>"
        marks.append(f"<Placemark><name>{name}</name><styleUrl>#burn</styleUrl>"
                     f"<description>{r.get('area_ha', '')} ha</description>{body}</Placemark>")
    path.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">'
        f"<Document><name>{name}</name><Style id=\"burn\"><LineStyle><color>ff0000ff</color>"
        "<width>2</width></LineStyle><PolyStyle><color>7f0000ff</color></PolyStyle></Style>"
        + "".join(marks) + "</Document></kml>", encoding="utf-8")


def compare_reference(footprint, ref_path: str, aoi_geom) -> dict:
    """Compare with an official perimeter (INFOEX in the course; EFFIS or
    Copernicus EMS for 2026), clipped to the AOI."""
    import geopandas as gpd
    ref = gpd.read_file(ref_path).to_crs(footprint.crs)
    R = _union(ref).intersection(aoi_geom)
    A = _union(footprint)
    inter = A.intersection(R).area
    return {"reference_ha_in_aoi": round(R.area / 1e4, 2),
            "footprint_ha": round(A.area / 1e4, 2),
            "overlap_ha": round(inter / 1e4, 2),
            "iou": round(inter / A.union(R).area, 3) if A.union(R).area else None,
            "commission_pct": round(100 * (A.area - inter) / A.area, 1) if A.area else None,
            "omission_pct": round(100 * (R.area - inter) / R.area, 1) if R.area else None}


def _union(gdf):
    return gdf.union_all() if hasattr(gdf, "union_all") else gdf.unary_union


# --------------------------------------------------------------------------- #
# 6. Impact by sector
# --------------------------------------------------------------------------- #
def landcover_impact(catalog, ds, mask, cfg: Config):
    import odc.stac
    import pandas as pd
    items = list(catalog.search(collections=["esa-worldcover"], bbox=cfg.bbox).items())
    if not items:
        log("WorldCover not found - skipping land-cover impact")
        return None
    newest = max(i.properties.get("start_datetime", "")[:4] for i in items)
    items = [i for i in items if i.properties.get("start_datetime", "")[:4] == newest]
    lc = odc.stac.load(items, bands=["map"], geobox=ds.odc.geobox,
                       resampling="mode")["map"].isel(time=0).values
    px = cfg.resolution ** 2 / 1e6
    rows = []
    for code, name in WORLDCOVER.items():
        in_aoi = int((lc == code).sum())
        if in_aoi:
            burned = int(((lc == code) & mask).sum())
            rows.append({"code": code, "class": name, "aoi_km2": round(in_aoi * px, 3),
                         "burned_km2": round(burned * px, 3),
                         "pct_of_class_burned": round(100 * burned / in_aoi, 1)})
    df = pd.DataFrame(rows).sort_values("burned_km2", ascending=False)
    tot = df["burned_km2"].sum()
    df["pct_of_footprint"] = (100 * df["burned_km2"] / tot).round(1) if tot else 0
    return df


def overlay_impact(name: str, path: str, footprint, out_dir: Path):
    import geopandas as gpd
    layer = gpd.read_file(path).to_crs(footprint.crs)
    fp = _union(footprint)
    hit = layer[layer.intersects(fp)].copy()
    label_col = next((c for c in layer.columns if c != "geometry" and layer[c].dtype == object), None)
    res = {"layer": name, "features_total": len(layer), "features_affected": len(hit)}
    if len(hit) and hit.geom_type.isin(["Polygon", "MultiPolygon"]).all():
        hit["burned_km2"] = hit.geometry.intersection(fp).area / 1e6
        hit["pct_burned"] = 100 * hit["burned_km2"] / (hit.area / 1e6)
        cols = ([label_col] if label_col else []) + ["burned_km2", "pct_burned"]
        tbl = hit[cols].sort_values("burned_km2", ascending=False).round(3)
        tbl.to_csv(out_dir / f"overlay_{name}.csv", index=False)
        res["burned_km2_total"] = round(float(hit["burned_km2"].sum()), 3)
        res["top"] = tbl.head(10).to_dict(orient="records")
    hit.to_file(out_dir / "burn_footprint.gpkg", layer=f"affected_{name}", driver="GPKG")
    return res


# --------------------------------------------------------------------------- #
# 7. Fire weather (ERA5 via Open-Meteo) vs 1991-2020
# --------------------------------------------------------------------------- #
def _om(lat, lon, start, end, daily):
    import pandas as pd
    import requests
    r = requests.get(OPEN_METEO_ARCHIVE, timeout=120, params={
        "latitude": lat, "longitude": lon, "start_date": start, "end_date": end,
        "daily": ",".join(daily), "timezone": "Europe/Madrid"})
    r.raise_for_status()
    df = pd.DataFrame(r.json()["daily"])
    df["time"] = pd.to_datetime(df["time"])
    return df.set_index("time")


def _pct(pm, months):
    sel = pm[pm.index.month.isin(months)]
    if sel.empty or sel["clim_mm"].sum() == 0:
        return None
    return round(float(100 * sel["obs_mm"].sum() / sel["clim_mm"].sum()), 0)


def fire_weather(cfg: Config, out_dir: Path) -> dict:
    import matplotlib.pyplot as plt
    import pandas as pd
    lat, lon = (cfg.bbox[1] + cfg.bbox[3]) / 2, (cfg.bbox[0] + cfg.bbox[2]) / 2
    fire = pd.Timestamp(cfg.fire_start)
    vars_ = ["temperature_2m_max", "relative_humidity_2m_min", "wind_speed_10m_max",
             "wind_gusts_10m_max", "precipitation_sum", "et0_fao_evapotranspiration"]
    obs = _om(lat, lon, f"{fire.year - 1}-10-01",
              (fire + pd.Timedelta(days=14)).date().isoformat(), vars_)
    clim = _om(lat, lon, "1991-01-01", "2020-12-31", ["temperature_2m_max", "precipitation_sum"])
    tclim = clim["temperature_2m_max"].groupby(clim.index.day_of_year).mean()
    obs["tmax_clim"] = tclim.reindex(obs.index.day_of_year).values
    obs["tmax_anom"] = obs["temperature_2m_max"] - obs["tmax_clim"]
    pm = obs["precipitation_sum"].resample("MS").sum()
    pclim = clim["precipitation_sum"].resample("MS").sum()
    pclim = pclim.groupby(pclim.index.month).mean()
    pm = pd.DataFrame({"obs_mm": pm, "clim_mm": pclim.reindex(pm.index.month).values})
    pm.to_csv(out_dir / "precip_monthly_vs_1991_2020.csv")

    w = obs.loc[fire - pd.Timedelta(days=7): fire + pd.Timedelta(days=5)]
    summary = {
        "point": [round(lat, 4), round(lon, 4)],
        "tmax_peak_c": round(float(w["temperature_2m_max"].max()), 1),
        "tmax_anom_mean_week_before_c": round(float(
            obs.loc[fire - pd.Timedelta(days=7): fire, "tmax_anom"].mean()), 1),
        "rh_min_lowest_pct": round(float(w["relative_humidity_2m_min"].min()), 0),
        "wind_gust_max_kmh": round(float(w["wind_gusts_10m_max"].max()), 0),
        "precip_jan_pct_of_normal": _pct(pm, [1]),
        "precip_jun_jul_pct_of_normal": _pct(pm, [6, 7]),
        "et0_sum_july_mm": round(float(obs.loc[obs.index.month == 7,
                                               "et0_fao_evapotranspiration"].sum()), 0),
    }
    win = obs.loc[fire - pd.Timedelta(days=40): fire + pd.Timedelta(days=14)]
    fig, ax = plt.subplots(3, 1, figsize=(10, 9), constrained_layout=True)
    ax[0].plot(win.index, win["temperature_2m_max"], color="#c0392b", label=f"Tmax {fire.year}")
    ax[0].plot(win.index, win["tmax_clim"], color="#7f8c8d", ls="--", label="Tmax 1991-2020")
    ax[0].fill_between(win.index, win["tmax_clim"], win["temperature_2m_max"],
                       where=win["tmax_anom"] > 0, color="#e74c3c", alpha=.2)
    ax[0].set_ylabel("°C"); ax[0].legend(loc="lower left")
    ax[1].plot(win.index, win["relative_humidity_2m_min"], color="#2980b9", label="RH min (%)")
    a1 = ax[1].twinx()
    a1.plot(win.index, win["wind_gusts_10m_max"], color="#8e44ad", label="Max gust (km/h)")
    ax[1].set_ylabel("RH %"); a1.set_ylabel("km/h")
    ax[1].legend(loc="upper left"); a1.legend(loc="upper right")
    for a in ax[:2]:
        a.axvline(fire, color="k", lw=1)
        a.text(fire, a.get_ylim()[1], " ignition", va="top", fontsize=8)
    x = np.arange(len(pm))
    ax[2].bar(x - .2, pm["clim_mm"], .4, color="#95a5a6", label="1991-2020")
    ax[2].bar(x + .2, pm["obs_mm"], .4, color="#16a085", label=f"{fire.year-1}-{fire.year}")
    ax[2].set_xticks(x, [d.strftime("%b %y") for d in pm.index], rotation=45)
    ax[2].set_ylabel("mm / month"); ax[2].legend()
    ax[2].set_title("Precipitation: wet winter -> dry summer (fuel build-up, then drying)")
    ax[0].set_title(f"Fire weather at AOI centre ({lat:.2f}N, {lon:.2f}E) - ERA5 via Open-Meteo")
    fig.savefig(out_dir / "09_fire_weather.png", dpi=150); plt.close(fig)
    return summary


# --------------------------------------------------------------------------- #
# 8. Figures
# --------------------------------------------------------------------------- #
def extent(ds):
    x, y = ds.x.values, ds.y.values
    r = abs(x[1] - x[0]) / 2
    return [x.min() - r, x.max() + r, y.min() - r, y.max() + r]


def stretch(img, ref=None):
    """2-98 % stretch per band; `ref` lets pre and post share one stretch so
    colours are comparable between dates."""
    ref = img if ref is None else ref
    out = np.empty_like(img)
    for i in range(img.shape[-1]):
        lo, hi = np.nanpercentile(ref[..., i], [2, 98])
        out[..., i] = np.clip((img[..., i] - lo) / (hi - lo + 1e-9), 0, 1)
    return np.nan_to_num(out)


def vis_nir(b):
    """Sharpest NIR band for display (Sentinel-2 B08 at 10 m, else Landsat B5)."""
    return b.get("nir10", b["nir"])


def fc543(b):
    """Course 5-4-3 false colour IR (Landsat 8/9): R=NIR, G=Red, B=Green."""
    return np.dstack([vis_nir(b), b["red"], b["green"]])


def make_figures(ds, pre, post, idx, seg_clean, final, sev, meta, cfg, out: Path):
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch
    ext = extent(ds)
    p0 = fc543(pre)
    post543 = stretch(fc543(post), p0)
    lbl = SENSORS[cfg.sensor]["band_numbers"]

    # 01 - false colour 5-4-3 pre/post
    fig, ax = plt.subplots(1, 2, figsize=(14, 6.5), constrained_layout=True)
    ax[0].imshow(stretch(p0), extent=ext); ax[0].set_title(f"Pre-fire {meta['pre_date']}")
    ax[1].imshow(post543, extent=ext); ax[1].set_title(f"Post-fire {meta['post_date']}")
    fig.suptitle(f"False-colour infrared 5-4-3 ({lbl}) - healthy vegetation = red")
    fig.savefig(out / "01_false_color_543.png", dpi=150); plt.close(fig)

    # 02 - NIR change composite: R = NIR pre, G = NIR post, B = NIR pre
    nir_comp = np.dstack([vis_nir(pre), vis_nir(post), vis_nir(pre)])
    lo, hi = np.nanpercentile(np.concatenate([vis_nir(pre).ravel(), vis_nir(post).ravel()]), [2, 98])
    nir_comp = np.nan_to_num(np.clip((nir_comp - lo) / (hi - lo), 0, 1))
    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
    ax.imshow(nir_comp, extent=ext)
    ax.legend(handles=[Patch(color="#b030b0", label="Purple = vegetation loss"),
                       Patch(color="#30b030", label="Green = vegetation gain"),
                       Patch(color="#9a9a9a", label="Grey = no change")],
              loc="lower right", fontsize=8)
    ax.set_title("NIR change composite (R: NIR pre, G: NIR post, B: NIR pre)")
    fig.savefig(out / "02_nir_change_composite.png", dpi=150); plt.close(fig)

    # 03 - NDVI Chuvieco classes pre/post
    cmap = ListedColormap(NDVI_COLORS); norm = BoundaryNorm(np.arange(-.5, 5.5), 5)
    fig, ax = plt.subplots(1, 2, figsize=(14, 6.5), constrained_layout=True)
    for a, key, d in [(ax[0], "NDVI_pre", meta["pre_date"]), (ax[1], "NDVI_post", meta["post_date"])]:
        a.imshow(np.ma.masked_equal(classify(idx[key], NDVI_BINS), 255), extent=ext,
                 cmap=cmap, norm=norm, interpolation="nearest")
        a.set_title(f"NDVI zoning (Chuvieco 2015) - {d}")
    ax[1].legend(handles=[Patch(color=c, label=l) for c, l in zip(NDVI_COLORS, NDVI_LABELS)],
                 loc="lower right", fontsize=7, framealpha=.9)
    fig.savefig(out / "03_ndvi_classes.png", dpi=150); plt.close(fig)

    # 04 - indices panel
    panels = [("dNDVI", "ΔNDVI = NDVIpost − NDVIpre", "RdYlGn", (-0.6, 0.3)),
              ("NBR_post", "NBR post", "RdYlGn", (-0.5, 0.8)),
              ("dBAI", "ΔBAI = BAIpost − BAIpre", "inferno", (-50, 400)),
              ("dGEMI", "ΔGEMI = GEMIpost − GEMIpre", "RdYlGn", (-0.4, 0.15)),
              ("dNBR", "dNBR = NBRpre − NBRpost", "RdYlGn_r", (-0.3, 1.0)),
              ("NDVI_pre", "NDVI pre", "RdYlGn", (-0.2, 0.9))]
    fig, axs = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True)
    for a, (k, t, cm, (vmin, vmax)) in zip(axs.ravel(), panels):
        im = a.imshow(idx[k], extent=ext, cmap=cm, vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=a, shrink=.8); a.set_title(t); a.set_xticks([]); a.set_yticks([])
    fig.savefig(out / "04_indices.png", dpi=130); plt.close(fig)

    # 05 - segmentations over post 5-4-3: one figure per threshold set
    titles = {"bastarrika": "Bastarrika et al. (2011) thresholds",
              "scene": f"scene-calibrated thresholds (unburned median ± {cfg.k_sigma}·MAD)"}
    for set_name, segs in seg_clean.items():
        keys = MERGE_INDICES[set_name]
        fig, axs = plt.subplots(2, 2, figsize=(13, 11), constrained_layout=True)
        for a, k in zip(axs.ravel(), keys):
            m = segs[k]
            op, t = meta["thresholds"][set_name][k]
            a.imshow(post543, extent=ext)
            a.imshow(np.ma.masked_where(~m, m), extent=ext, cmap=ListedColormap(["#ff1a1a"]),
                     alpha=.8, interpolation="nearest")
            a.set_title(f"{k} {op} {t}  -> {area_stats(m, cfg.resolution)['area_ha']:,.0f} ha")
            a.set_xticks([]); a.set_yticks([])
        fig.suptitle(f"Segmentation, {titles[set_name]}, after majority filter")
        fig.savefig(out / f"05_segmentations_{set_name}.png", dpi=130); plt.close(fig)

    # 06 - final footprint over post 5-4-3
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    ax.imshow(post543, extent=ext)
    if len(meta["gdf"]):
        meta["gdf"].boundary.plot(ax=ax, color="#ffff00", linewidth=0.9)
    ax.set_title(f"Burn footprint ({meta['method']}): {meta['area_ha']:,.1f} ha = "
                 f"{meta['area_km2']:.2f} km²  |  perimeter {meta['perimeter_km']:.1f} km")
    fig.savefig(out / "06_footprint.png", dpi=150); plt.close(fig)

    # 07 - dNBR severity
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    ax.imshow(np.ma.masked_equal(sev, 255), extent=ext, cmap=ListedColormap(SEVERITY_COLORS),
              norm=BoundaryNorm(np.arange(-.5, 7.5), 7), interpolation="nearest")
    ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(SEVERITY_COLORS, SEVERITY_LABELS)],
              loc="lower right", fontsize=7, framealpha=.9)
    ax.set_title("Burn severity (dNBR, USGS / Key & Benson 2006)")
    fig.savefig(out / "07_severity.png", dpi=150); plt.close(fig)


def write_multiband_tif(ds, bands: list, path: Path):
    """3-band float32 GeoTIFF (e.g. 5-4-3 composite) that ArcGIS/QGIS can display
    as RGB with its own stretch."""
    import rioxarray  # noqa: F401
    import xarray as xr
    arr = np.stack(bands).astype("float32")
    da = xr.DataArray(arr, coords={"band": np.arange(1, len(bands) + 1), "y": ds.y, "x": ds.x},
                      dims=("band", "y", "x"))
    da = da.rio.write_crs(ds.odc.geobox.crs.to_wkt()).rio.write_nodata(np.nan)
    da.rio.to_raster(path, compress="deflate")


def write_tif(ds, arr, path: Path, nodata):
    import rioxarray  # noqa: F401
    import xarray as xr
    da = xr.DataArray(arr, coords={"y": ds.y, "x": ds.x}, dims=("y", "x"))
    da = da.rio.write_crs(ds.odc.geobox.crs.to_wkt()).rio.write_nodata(nodata)
    da.rio.to_raster(path, compress="deflate")


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #
def run(cfg: Config) -> dict:
    import pandas as pd
    from shapely.geometry import box
    out = Path(cfg.out); out.mkdir(parents=True, exist_ok=True)
    catalog = open_catalog()
    log(f"Sensor: {cfg.sensor}  ({cfg.resolution} m, {cfg.crs})")

    # 1. Scenes
    pre_date, pre_items, pre_score = pick_pre(catalog, cfg)
    post_date, post_items, post_score = pick_post(catalog, cfg)
    log(f"Selected PRE  {pre_date} {pre_score}")
    log(f"Selected POST {post_date} {post_score}")

    # 2. Reflectance
    ds, pre = load_scene(pre_items, cfg)
    _, post = load_scene(post_items, cfg, geobox=ds.odc.geobox)
    qa0, qa1 = qa_masks(cfg.sensor, pre["qa"]), qa_masks(cfg.sensor, post["qa"])
    valid = ~(qa0["nodata"] | qa1["nodata"] | qa0["cloud"] | qa1["cloud"]
              | qa0["snow"] | qa1["snow"] | qa0["water"])
    for b in ("red", "nir", "swir2"):
        valid &= np.isfinite(pre[b]) & np.isfinite(post[b])

    # 3. Indices
    idx = compute_indices(pre, post)
    for k in idx:
        idx[k] = np.where(valid, idx[k], np.nan).astype("float32")

    # 4. Segmentation + majority filter + per-index areas
    min_px = max(1, int(round(cfg.min_patch_ha * 1e4 / cfg.resolution ** 2)))
    nb, mode = cfg.majority_neighbors, cfg.majority_threshold
    # Both threshold sets are always computed and reported, so the course's
    # Bastarrika results can be compared with the scene-calibrated ones.
    thr_sets = {"bastarrika": {**BASTARRIKA, "dNBR": (">=", cfg.dnbr_threshold)},
                "scene": scene_thresholds(idx, valid, cfg.k_sigma)}
    log("Scene-calibrated thresholds (unburned median ± "
        f"{cfg.k_sigma}·MAD): {thr_sets['scene']}")
    sep_rows, area_rows, segs_raw, segs_clean = [], [], {}, {}
    for set_name, thr in thr_sets.items():
        segs_raw[set_name] = {k: threshold(idx[k], op if op != ">=" else ">", t) & valid
                              for k, (op, t) in thr.items()}
        segs_clean[set_name] = {k: clean(m, valid, min_px, nb, mode)
                                for k, m in segs_raw[set_name].items()}
        for k, m in segs_clean[set_name].items():
            s = separability_m(idx[k], m, valid)
            if s:
                sep_rows.append({"thresholds": set_name, "index": k,
                                 **{kk: round(v, 6) for kk, v in s.items()}})
            area_rows.append({"thresholds": set_name, "index": k,
                              "rule": f"{k} {thr[k][0]} {thr[k][1]}",
                              **area_stats(m, cfg.resolution)})
    thr = thr_sets[cfg.thresholds]
    seg_raw, seg_clean = segs_raw[cfg.thresholds], segs_clean[cfg.thresholds]
    merge_keys = MERGE_INDICES[cfg.thresholds]

    # 5. Final footprint
    if cfg.method == "merge":
        # strong char for patch cores = the course's fixed ΔBAI > 144.38 (Bastarrika)
        strong = segs_raw["bastarrika"][CHAR_INDEX]
        final, votes = merge_and_clean({k: seg_raw[k] for k in merge_keys}, valid,
                                       cfg.min_votes, min_px, nb, mode, cfg.pixel_votes,
                                       strong_char=strong)
        write_tif(ds, votes.astype("uint8"), out / "index_votes.tif", 255)
        rule = (f"pixel: >= {cfg.pixel_votes} of {', '.join(merge_keys)}; patch core "
                f">= {cfg.min_patch_ha} ha with >= {cfg.min_votes} votes and "
                f"{CHAR_INDEX} {BASTARRIKA[CHAR_INDEX][0]} {BASTARRIKA[CHAR_INDEX][1]}")
    else:
        final = seg_clean[cfg.method]
        rule = f"{cfg.method} {thr[cfg.method][0]} {thr[cfg.method][1]}"
    stats = area_stats(final, cfg.resolution)
    area_rows.append({"thresholds": cfg.thresholds, "index": f"FINAL ({cfg.method})",
                      "rule": rule, **stats})
    s = separability_m(idx["dNBR"], final, valid)
    if s:
        sep_rows.append({"thresholds": cfg.thresholds, "index": "dNBR vs final footprint",
                         **{k: round(v, 6) for k, v in s.items()}})

    pd.DataFrame(sep_rows).to_csv(out / "separability_M.csv", index=False)
    pd.DataFrame(area_rows).to_csv(out / "index_areas.csv", index=False)

    # Vectors: dissolve, per-index layers, KML
    tr = geobox_transform(ds)
    gdf = vectorise(final, tr, cfg.crs)
    gpkg = out / "burn_footprint.gpkg"
    if gpkg.exists():
        gpkg.unlink()
    gdf.to_file(gpkg, layer="footprint_dissolved", driver="GPKG")
    # The exact processed extent (UTM-aligned grid). It is slightly larger than the
    # lon/lat --bbox and slightly rotated relative to it, so this is the AOI to map.
    import geopandas as gpd
    e = extent(ds)
    gpd.GeoDataFrame({"area_km2": [round((e[1] - e[0]) * (e[3] - e[2]) / 1e6, 2)]},
                     geometry=[box(e[0], e[2], e[1], e[3])], crs=cfg.crs
                     ).to_file(gpkg, layer="aoi_extent", driver="GPKG")
    for set_name, segs in segs_clean.items():
        for k, m in segs.items():
            vectorise(m, tr, cfg.crs).to_file(gpkg, layer=f"segmentation_{set_name}_{k}",
                                              driver="GPKG")
    gdf.to_crs(4326).to_file(out / "burn_footprint.geojson", driver="GeoJSON")
    write_kml(gdf, out / "burn_footprint.kml", f"Burn footprint {cfg.fire_start[:4]}")
    perimeter_km = round(float(gdf.length.sum() / 1e3), 2) if len(gdf) else 0.0  # incl. holes
    outer_km = outer_perimeter_km(gdf)

    # Rasters
    for k in ("NDVI_pre", "NDVI_post", "dNDVI", "NBR_pre", "NBR_post", "dNBR", "RdNBR",
              "BAI_pre", "BAI_post", "dBAI", "GEMI_pre", "GEMI_post", "dGEMI"):
        write_tif(ds, idx[k], out / f"{k}.tif", np.nan)
    # Reflectance composites for ArcGIS / QGIS (band 1 = NIR, 2 = Red, 3 = Green)
    write_multiband_tif(ds, [vis_nir(pre), pre["red"], pre["green"]], out / "pre_fire_543.tif")
    write_multiband_tif(ds, [vis_nir(post), post["red"], post["green"]], out / "post_fire_543.tif")
    # Course 5.1 NIR change composite: R = NIR pre, G = NIR post, B = NIR pre
    write_multiband_tif(ds, [vis_nir(pre), vis_nir(post), vis_nir(pre)],
                        out / "nir_change_composite.tif")
    sev = classify(idx["dNBR"], SEVERITY_BINS)
    write_tif(ds, sev, out / "severity_dNBR.tif", 255)
    write_tif(ds, final.astype("uint8"), out / "footprint_mask.tif", 255)

    sev_df = class_area_table(sev, SEVERITY_LABELS, cfg.resolution, final)
    sev_df.to_csv(out / "severity_areas.csv", index=False)
    ndvi_df = class_area_table(classify(idx["NDVI_pre"], NDVI_BINS), NDVI_LABELS, cfg.resolution)
    ndvi_df = ndvi_df.rename(columns={"aoi_km2": "pre_km2"})
    ndvi_df["post_km2"] = class_area_table(classify(idx["NDVI_post"], NDVI_BINS),
                                           NDVI_LABELS, cfg.resolution)["aoi_km2"]
    ndvi_df["change_km2"] = (ndvi_df["post_km2"] - ndvi_df["pre_km2"]).round(3)
    ndvi_df.to_csv(out / "ndvi_class_areas.csv", index=False)

    aoi_km2 = round(valid.size * cfg.resolution ** 2 / 1e6, 2)
    summary = {
        "config": {k: v for k, v in asdict(cfg).items() if k != "overlays"},
        "pre_scene": {"date": pre_date, **pre_score, "ids": [i.id for i in pre_items]},
        "post_scene": {"date": post_date, **post_score, "ids": [i.id for i in post_items]},
        "aoi_km2": aoi_km2, "valid_pct_of_aoi": round(float(valid.mean() * 100), 2),
        "burn_footprint": {**stats, "outer_perimeter_km": outer_km,
                           "total_edge_length_km_incl_holes": perimeter_km, "n_polygons":
                           int(len(gdf.explode(index_parts=False))) if len(gdf) else 0,
                           "pct_of_aoi": round(100 * stats["area_km2"] / aoi_km2, 1)},
        "thresholds_used": cfg.thresholds, "threshold_sets": thr_sets,
        "index_areas": area_rows, "separability_M": sep_rows,
        "severity": sev_df.to_dict(orient="records"),
        "ndvi_classes": ndvi_df.to_dict(orient="records"),
    }

    if cfg.reference:
        try:
            summary["reference_comparison"] = compare_reference(gdf, cfg.reference, box(*[
                extent(ds)[i] for i in (0, 2, 1, 3)]))
        except Exception as e:
            log(f"reference comparison failed: {e}")
    if cfg.landcover:
        try:
            lc = landcover_impact(catalog, ds, final, cfg)
            if lc is not None:
                lc.to_csv(out / "landcover_impact.csv", index=False)
                summary["landcover_impact"] = lc.to_dict(orient="records")
        except Exception as e:
            log(f"land-cover step failed: {e}")
    if cfg.overlays:
        summary["overlays"] = []
        for name, path in cfg.overlays.items():
            try:
                summary["overlays"].append(overlay_impact(name, path, gdf, out))
            except Exception as e:
                log(f"overlay {name} failed: {e}")
    if cfg.weather:
        try:
            summary["fire_weather"] = fire_weather(cfg, out)
        except Exception as e:
            log(f"weather step failed: {e}")

    meta = {"pre_date": pre_date, "post_date": post_date, "method": cfg.method, "gdf": gdf,
            "area_ha": stats["area_ha"], "area_km2": stats["area_km2"],
            "perimeter_km": outer_km}
    meta["thresholds"] = thr_sets
    make_figures(ds, pre, post, idx, segs_clean, final, sev, meta, cfg, out)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    log("-" * 64)
    log("\n" + pd.DataFrame(area_rows)[["thresholds", "index", "rule", "area_ha", "area_km2"]]
        .to_string(index=False))
    log("\n" + pd.DataFrame(sep_rows)[["thresholds", "index", "M"]].to_string(index=False))
    log(f"FINAL burned footprint: {stats['area_ha']:,} ha = {stats['area_km2']} km²")
    log(f"Outer perimeter: {outer_km} km  (all edges incl. unburned holes: {perimeter_km} km)")
    log(f"Outputs in {out.resolve()}")
    return summary


def parse_args(argv=None) -> Config:
    d = Config()
    p = argparse.ArgumentParser(description="Automated pre/post-fire burn footprint",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--sensor", choices=list(SENSORS), default=d.sensor)
    p.add_argument("--bbox", nargs=4, type=float, default=d.bbox,
                   metavar=("LON_MIN", "LAT_MIN", "LON_MAX", "LAT_MAX"))
    p.add_argument("--fire-start", default=d.fire_start)
    p.add_argument("--days-before", type=int, default=d.days_before)
    p.add_argument("--post-start", default=d.post_start)
    p.add_argument("--post-end", default=d.post_end)
    p.add_argument("--pre-date"); p.add_argument("--post-date")
    p.add_argument("--max-aoi-cloud", type=float, default=d.max_aoi_cloud)
    p.add_argument("--method", default=d.method,
                   choices=["merge", "dNDVI", "NBR_post", "dBAI", "dGEMI", "dNBR"])
    p.add_argument("--thresholds", choices=["scene", "bastarrika"], default=d.thresholds,
                   help="scene = calibrated on this image's unburned pixels; "
                        "bastarrika = fixed values from the course")
    p.add_argument("--k-sigma", type=float, default=d.k_sigma)
    p.add_argument("--pixel-votes", type=int, default=d.pixel_votes, choices=[1, 2, 3, 4])
    p.add_argument("--min-votes", type=int, default=d.min_votes, choices=[1, 2, 3, 4])
    p.add_argument("--dnbr-threshold", type=float, default=d.dnbr_threshold)
    p.add_argument("--majority-neighbors", type=int, choices=[4, 8], default=d.majority_neighbors)
    p.add_argument("--majority-threshold", choices=["majority", "half"], default=d.majority_threshold)
    p.add_argument("--min-patch-ha", type=float, default=d.min_patch_ha)
    p.add_argument("--resolution", type=int)
    p.add_argument("--reference", help="official perimeter (shp/gpkg/geojson) to validate against")
    p.add_argument("--out", default=d.out)
    p.add_argument("--no-weather", action="store_true")
    p.add_argument("--no-landcover", action="store_true")
    p.add_argument("--overlay", action="append", default=[], metavar="NAME=PATH")
    a = p.parse_args(argv)
    overlays = {}
    for o in a.overlay:
        name, _, path = o.partition("=")
        if not path:
            p.error(f"--overlay expects NAME=PATH, got {o!r}")
        overlays[name] = path
    return Config(sensor=a.sensor, bbox=a.bbox, fire_start=a.fire_start,
                  days_before=a.days_before, post_start=a.post_start, post_end=a.post_end,
                  pre_date=a.pre_date, post_date=a.post_date, max_aoi_cloud=a.max_aoi_cloud,
                  method=a.method, thresholds=a.thresholds, k_sigma=a.k_sigma,
                  pixel_votes=a.pixel_votes, min_votes=a.min_votes,
                  dnbr_threshold=a.dnbr_threshold,
                  majority_neighbors=a.majority_neighbors,
                  majority_threshold=a.majority_threshold, min_patch_ha=a.min_patch_ha,
                  resolution=a.resolution, reference=a.reference, out=a.out,
                  weather=not a.no_weather, landcover=not a.no_landcover, overlays=overlays)


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    run(parse_args())
