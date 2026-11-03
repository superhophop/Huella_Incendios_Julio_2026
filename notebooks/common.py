# =============================================================================
#  common.py - settings and small helpers shared by every notebook
#
#  Each notebook starts with   from common import *
#  so all of them use the same area, dates and output folder. Change a setting
#  here once and every notebook picks it up (restart the kernel after editing).
# =============================================================================
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)   # NaN maths in empty areas

# -----------------------------------------------------------------------------
#  SETTINGS - change these to study another fire or area
# -----------------------------------------------------------------------------
CFG = {
    # Area of interest (lon_min, lat_min, lon_max, lat_max), WGS84 degrees.
    # The whole July 2026 fire complex in Ávila, Madrid and Toledo, with a margin so
    # the burn scar doesn't touch the edge (notebook 04 checks this). The first version
    # used a smaller area, [-4.46, 40.26, -4.17, 40.46], covering only the Madrid side.
    "bbox": [-4.90, 40.20, -4.10, 40.56],
    "fire_start": "2026-07-22",   # first ignitions (Burgohondo, Ávila; Almorox, Toledo)
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
}

# Folders: the project root is the folder above /notebooks
PROJECT = Path(__file__).resolve().parent.parent
OUT = PROJECT / "outputs_madrid_avila"
OUT.mkdir(exist_ok=True)

COLLECTION = "sentinel-2-l2a"
STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"

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


# -----------------------------------------------------------------------------
#  Catalogue access
# -----------------------------------------------------------------------------
def open_catalog():
    """Connect to Planetary Computer. sign_inplace adds a temporary access
    token to every band link so the files can be downloaded."""
    import planetary_computer
    import pystac_client
    return pystac_client.Client.open(STAC_URL, modifier=planetary_computer.sign_inplace)


# -----------------------------------------------------------------------------
#  GeoTIFF helpers
#  A GeoTIFF = pixel values + "where on Earth" (CRS + transform). The transform
#  says: top-left corner X/Y and pixel size, so ArcGIS can place every pixel.
# -----------------------------------------------------------------------------
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


def need(name: str, notebook: str) -> Path:
    """Path to an output file, with a clear message if an earlier notebook
    has not been run yet."""
    p = OUT / name
    if not p.exists():
        raise FileNotFoundError(f"{p.name} not found. Run notebook {notebook} first.")
    return p


def read_reflectance(which: str) -> tuple:
    """Load pre_ or post_reflectance.tif as {band role: 2-D array}."""
    stack, tr, crs, _ = read_tif(need(f"{which}_reflectance.tif", "02"))
    return {role: stack[i] for i, role in enumerate(BANDS)}, tr, crs


def aoi_polygon(crs):
    """The AOI rectangle in another CRS. Extra points along its sides keep them
    accurate, because a straight lon/lat line is slightly curved in UTM."""
    import geopandas as gpd
    from shapely import segmentize
    from shapely.geometry import box
    return gpd.GeoSeries([segmentize(box(*CFG["bbox"]), 0.001)], crs=4326).to_crs(crs).iloc[0]


def aoi_mask(shape, transform, crs) -> np.ndarray:
    """True for pixels fully inside the AOI. The image grid is a UTM rectangle
    slightly larger than the AOI, so its corners stick out past the AOI."""
    from rasterio import features
    aoi = aoi_polygon(crs)
    inside = features.geometry_mask([aoi], shape, transform, invert=True)
    on_edge = features.geometry_mask([aoi.boundary], shape, transform, invert=True,
                                     all_touched=True)
    return inside & ~on_edge


def read_valid():
    """Usable pixels: no cloud or water (from notebook 02), and inside the AOI."""
    mask, transform, crs, _ = read_tif(need("valid_mask.tif", "02"))
    return (mask == 1) & aoi_mask(mask.shape, transform, crs)


def secret(name: str) -> str:
    """Read a private value (like an API key) from the .env file in the project
    folder. .env is listed in .gitignore, so it never gets published."""
    env = PROJECT / ".env"
    if env.exists():
        for line in env.read_text(encoding="utf-8").splitlines():
            key, _, value = line.partition("=")
            if key.strip() == name:
                return value.strip()
    raise KeyError(f"Add a line {name}=... to {env}")


def save_json(obj, name: str):
    (OUT / name).write_text(json.dumps(obj, indent=2))


def load_json(name: str, notebook: str):
    return json.loads(need(name, notebook).read_text())


# -----------------------------------------------------------------------------
#  Display helpers
# -----------------------------------------------------------------------------
def stretch_rgb(rgb: np.ndarray, lo, hi) -> np.ndarray:
    """Scale reflectance to 0-1 brightness for display (like ArcGIS stretch)."""
    return np.nan_to_num(np.clip((rgb - lo) / (hi - lo), 0, 1))


def false_color(pre: dict, post: dict):
    """8-4-3 false colour (NIR, Red, Green) for both dates, using ONE stretch
    taken from the pre-fire image so the two pictures are comparable."""
    rgb0 = np.dstack([pre["nir10"], pre["red"], pre["green"]])
    rgb1 = np.dstack([post["nir10"], post["red"], post["green"]])
    lo = np.nanpercentile(rgb0, 2, axis=(0, 1))
    hi = np.nanpercentile(rgb0, 98, axis=(0, 1))
    return stretch_rgb(rgb0, lo, hi), stretch_rgb(rgb1, lo, hi), lo, hi


def aoi_map(zoom: int = 10):
    """Interactive web map (folium) with the AOI rectangle drawn on it."""
    import folium
    x0, y0, x1, y1 = CFG["bbox"]
    m = folium.Map(location=[(y0 + y1) / 2, (x0 + x1) / 2], zoom_start=zoom,
                   tiles="OpenStreetMap")
    folium.TileLayer(
        "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        attr="Esri World Imagery", name="Satellite").add_to(m)
    folium.Rectangle([[y0, x0], [y1, x1]], color="red", fill=False, tooltip="AOI").add_to(m)
    return m
