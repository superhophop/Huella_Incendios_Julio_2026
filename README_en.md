# Burn footprint of the July 2026 fires in Ávila, Madrid and Toledo

Pre-fire and post-fire burn mapping of the fire complex in central Spain, using Sentinel-2 images
(24 June and 8 August 2026). The method automates the change-detection workflow of the course
practical *Aplicación 3 – Detección de cambios (Sierra de Gata)*: what was done step by step in
ArcGIS Pro (Raster Calculator, Composite Bands, Reclassify, Majority Filter, Raster to Polygon,
Merge, Dissolve, Layer To KML) runs here in Python.

## Results

| Measure | Result |
|---|---|
| Study area (AOI) | lon −4.90…−4.10, lat 40.20…40.56, 2,714.64 km² |
| Mapped burn footprint | **618.03 km²** (61,803 ha), with no burned pixels at the AOI edge |
| By province | Ávila 386.07 km², Madrid 219.41 km², Toledo 12.54 km² |
| Municipalities touched | 30 of 78 |
| Separability | dNBR separability M = 1.513 (above 1 means good separation) |
| Official perimeter (EFFIS) | 641.02 km² in 7 fires; 97.9% of the footprint lies inside it, 94.3% of it is captured, Dice 0.961 |
| Natura 2000 | 502.73 km² (81.3% of the footprint) inside protected sites; EFFIS gives 81.7% |
| Spread (VIIRS hotspots) | First detection 22 July 15:35 local; peak day 25 July, about 200 km² |

The footprint agrees closely with the official EFFIS perimeters (notebook 07). Press figures (more than
500 km² in Ávila and 268.6 km² in Madrid) are higher than both maps and come from other sources and dates.

## 0. Notebooks (recommended way to work)

The workflow is in `notebooks/` as step-by-step teaching notebooks (`00_start_here` → `08_progression`),
with shared settings in `notebooks/common.py` and results written to `outputs_madrid_avila/`.
`notebooks/fetch_overlays.py` downloads the layers that notebook 05 overlays.

Notebook 08 needs a free NASA FIRMS key (https://firms.modaps.eosdis.nasa.gov/api/map_key/). Put it in
a file called `.env` in the project folder, as `FIRMS_MAP_KEY=your_key`. `.env` is listed in `.gitignore`,
so the key is never published.

Open the `incendio_2026` folder in VS Code and the notebooks use the project's `.venv` automatically.
To rebuild that environment:

```powershell
& "$env:USERPROFILE\anaconda3\python.exe" -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
```

The standalone scripts described below are the earlier version of this project. They still default
to the smaller Madrid-only area.

## 1. Install (scripts)

```bash
conda create -n burn -c conda-forge python=3.11 pystac-client planetary-computer \
    odc-stac rioxarray rasterio geopandas scipy matplotlib pandas requests
conda activate burn
```

No account needed. Landsat 8/9 Collection 2 Level-2 (the same L2SP product as in the course)
and Sentinel-2 L2A are read straight from Microsoft Planetary Computer, so there is no
EarthExplorer download or .tar unpacking.

## 2. Run

```bash
python burn_footprint.py                     # Landsat 8/9, 30 m (default)
python burn_footprint.py --sensor sentinel2  # Sentinel-2, 20 m, more frequent revisit
```

Defaults for the 2026 fire:

| Setting | Value |
|---|---|
| AOI | lon −4.46…−4.17, lat 40.26…40.46 (San Martín de Valdeiglesias, Pelayos de la Presa, Navas del Rey, Chapinería, Villa del Prado, Aldea del Fresno) |
| Ignition | 22 Jul 2026 (Almorox; San Martín de Valdeiglesias on 23 Jul) |
| Pre-fire scene | cloud-free date closest to **24 Jun 2026** (ignition − 28 days), searched ±7 days, widening up to ±28 |
| Post-fire scene | earliest cloud-free date from 30 Jul to 15 Sep, i.e. Key & Benson's *initial assessment* once the fire is out |
| Cloud-free | ≤ 1 % cloud/shadow pixels **inside the AOI** (Landsat QA_PIXEL bits 1–4) and ≥ 99 % data coverage |
| Projection | WGS 84 / UTM 30N (EPSG:32630), the Madrid equivalent of the course's UTM 29N |

Landsat 8 + 9 together revisit every 8 days. If the ±7-day window has no clear scene, the
script widens it automatically. The candidate dates and their AOI cloud % are printed, and
the chosen scene IDs are stored in `summary.json`.

## 3. How each course step maps to the script

| Course step | Script | Output |
|---|---|---|
| 4. DN → reflectance: `(SR_Bx × 0.0000275) − 0.2`, then `Con(x<0, 0, x)` | `to_reflectance()` | – |
| 5.1 5-4-3 false-colour IR | `fc543()` (R=B5 NIR, G=B4 Red, B=B3 Green) | `01_false_color_543.png` |
| 5.1 NIR change composite (R: NIR pre, G: NIR post, B: NIR pre) | – | `02_nir_change_composite.png` (purple = loss) |
| 5.2-A NDVI + Chuvieco (2015) reclassification | `ndvi()`, `classify()` | `03_ndvi_classes.png`, `ndvi_class_areas.csv` |
| 5.2-B ΔNDVI = NDVIpost − NDVIpre | `compute_indices()` | `dNDVI.tif` |
| 5.2-C NBR = (B5 − B7)/(B5 + B7) | `nbr()` | `NBR_pre/post.tif` |
| 5.2-D/E BAI (Martín 1998), ΔBAI | `bai()` | `BAI_*.tif`, `dBAI.tif` |
| 5.2-F/G GEMI (Pinty & Verstraete 1992), ΔGEMI | `gemi()` | `GEMI_*.tif`, `dGEMI.tif` |
| 6. Segmentation, Bastarrika et al. (2011): ΔNDVI < −0.17767, NBRpost < −0.15006, ΔBAI > 144.3835, ΔGEMI < −0.11578 | `BASTARRIKA`, `threshold()` | `05_segmentations.png` |
| 6. Validation: M = \|μb − μu\| / (σb + σu) | `separability_m()` | `separability_M.csv` |
| 6. Majority Filter (8 neighbours, MAJORITY) | `majority_filter()` | applied to every segmentation |
| 6. Raster to Polygon + Sum of area per index | `vectorise()` | `index_areas.csv`, GPKG layers `segmentation_*` |
| 7. Merge + manual clean-up + Dissolve | `merge_and_clean()` | `burn_footprint.gpkg` → `footprint_dissolved` |
| 7. Layer To KML (Google Earth) | `write_kml()` | `burn_footprint.kml` |
| Compare with the official perimeter (INFOEX / Copernicus EMS in the course) | `--reference` | IoU, omission %, commission % in `summary.json` |

Extras beyond the course:

- dNBR burn severity (USGS / Key & Benson 2006): `07_severity.png` and `severity_areas.csv`
- burned area by land cover: `landcover_impact.csv`
- fire weather against the 1991–2020 climatology: `09_fire_weather.png`

### Threshold sets: course (Bastarrika) vs scene-calibrated (default)

Both sets are always computed and reported side by side (`index_areas.csv`,
`separability_M.csv`, `05_segmentations_bastarrika.png`, `05_segmentations_scene.png`).

- **`bastarrika`**: the course's fixed values. They were calibrated on a dense shrub fire,
  and in the 2026 Sierra Oeste run they caught only the intense core: 119 km².
- **`scene`** (default): each index's threshold = median ± 3·MAD (robust standard deviation)
  of this image pair's own unburned pixels (the USGS "unburned" dNBR class). This catches
  low-severity grass and dehesa burns. Robust statistics stop a few extreme pixels (terrain
  shadow, reservoir edges) from inflating the thresholds. In the 2026 run it gave 279.5 km²
  (Landsat 8, 30 m) and 259.8 km² (Sentinel-2, 10 m), with < 1 % false positives in unburned
  control zones (SW mountains, SE farmland).

Use `--thresholds bastarrika` to reproduce the course exactly.

### How the manual clean-up is automated

In the course, polygons outside the real burn (red, healthy patches in the 5-4-3 composite)
were deleted by hand after the Merge. The script decides this with the indices themselves:

- **Vegetation-loss indices:** ΔNDVI and ΔGEMI use only red and NIR. A cereal field
  harvested in July looks burned to them.
- **Char indices:** NBR post (SWIR) and ΔBAI (distance to the charcoal spectral point)
  respond to char, not to harvest.
- **The rule:** a pixel enters the footprint when ≥ `--pixel-votes` (default 2) of ΔNDVI,
  ΔGEMI, ΔBAI and dNBR flag it. A patch is kept only if it has a core of ≥ 1 ha where
  ≥ `--min-votes` (default 3) agree **and** there is strong char (ΔBAI > 144.38, the
  course value).
- **Why it works:** low-severity grass burns and fresh stubble look alike spectrally, but a
  real fire patch always contains some truly charred ground and a harvested field does not.

In the synthetic test this kept 100 % of a low-severity burn ring and removed 100 % of an
isolated harvested field. A harvested field that *touches* the burn can still be included,
so check the 5-4-3 composite along farmland edges.

## 4. Options

```bash
--method merge|dNDVI|NBR_post|dBAI|dGEMI|dNBR   # final footprint from one index or the merge
--min-votes 3                                   # stricter merge
--majority-neighbors 4 --majority-threshold half
--pre-date 2026-06-24 --post-date 2026-08-07    # force scenes
--reference effis_or_ems_perimeter.shp          # validation
--overlay municipios=recintos_municipales.shp   # burned km² per municipality
--overlay edificios=catastro_buildings.gpkg     # buildings inside the footprint
--overlay natura=natura2000.shp                 # protected habitat burned
--bbox -4.50 40.25 -4.10 40.50                  # different AOI
```

## 5. Data sources for overlays and validation

- **Official perimeter:** EFFIS burnt areas (Copernicus), or the Copernicus EMS Rapid
  Mapping delineation, the same type of product shown at the end of the course
- **Municipalities:** CNIG *Líneas límite municipales*
- **Buildings:** Catastro INSPIRE Buildings (ATOM feed per municipality)
- **Protected areas:** Natura 2000 (MITECO / EEA)
- **National fire statistics:** MITECO *Estadísticas de incendios forestales* and its *avances informativos*

## Caveats

- The scripts' default AOI covers only the Madrid side. The notebooks cover the whole complex.
- The Bastarrika thresholds were calibrated on Landsat TM in Mediterranean areas. Check the
  M values: if an index has M < 1 in your scene, its threshold is not separating the classes
  well, so drop it (`--method`) or adjust it.
- Key & Benson's *extended assessment* (post scene in the next growing season, June–July 2027)
  cannot be made yet. Rerun the script with `--post-start 2027-06-01 --post-end 2027-07-31`
  next summer to see delayed mortality and regrowth.
- Burned houses are only partly visible at 30 m. Use the Catastro overlay to count buildings.
