"""Download the overlay layers notebook 05 uses into ../data, for the AOI in common.py.

Run after notebook 04, from the notebooks folder: python fetch_overlays.py
Buildings are only fetched for municipalities the footprint touches, so it needs the footprint first.
"""
import re
import subprocess
import sys
import time
import unicodedata
import urllib.parse
import zipfile

import geopandas as gpd
import pandas as pd
from shapely import segmentize
from shapely.geometry import box

from common import CFG, OUT, PROJECT

sys.stdout.reconfigure(encoding="ascii", errors="replace")
ROOT = PROJECT / "data_downloads"
ROOT.mkdir(exist_ok=True)
W, S, E, N = CFG["bbox"]
DATA = PROJECT / "data"
DATA.mkdir(exist_ok=True)
UA = "Mozilla/5.0 (research script)"
BB84 = f"{S},{W},{N},{E},urn:ogc:def:crs:EPSG::4326"


def curl(url, dest=None, timeout=240):
    cmd = ["curl", "-sS", "-L", "-g", "-m", str(timeout), "-A", UA, "-w", "%{http_code}"]
    if dest:
        cmd += ["-o", str(dest)]
    return subprocess.run(cmd + [url.replace(" ", "%20")], capture_output=True).stdout


def plain(s):
    s = "".join(c for c in unicodedata.normalize("NFD", s) if unicodedata.category(c) != "Mn").upper()
    return re.sub(r"[^A-Z ]", " ", s).split()


def paged(base, typename, step, name):
    frames, start = [], 0
    page = ROOT / f"{name}_page.gml"
    while True:
        q = urllib.parse.urlencode({"service": "WFS", "version": "2.0.0", "request": "GetFeature", "typenames": typename, "count": str(step),
                                    "startIndex": str(start), "srsName": "urn:ogc:def:crs:EPSG::4326", "bbox": BB84})
        for attempt in range(6):
            curl(f"{base}?{q}", page, 200)
            try:
                g = gpd.read_file(page)
                break
            except Exception:
                print(" retry", name, start, attempt, flush=True)
                time.sleep(15)
        else:
            raise SystemExit(f"{name} gave up at {start}")
        frames.append(g)
        print(" ", name, start, len(g), flush=True)
        if len(g) < step:
            return gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), geometry="geometry", crs=frames[0].crs).drop_duplicates(subset="gml_id")
        start += step


fp = gpd.read_file(OUT / "burn_footprint.gpkg", layer="footprint_dissolved")
fpu = fp.union_all()
aoi = gpd.GeoSeries([segmentize(box(W, S, E, N), 0.001)], crs=4326).to_crs(fp.crs).iloc[0]

# municipalities
q = urllib.parse.urlencode({"service": "WFS", "version": "2.0.0", "request": "GetFeature", "typenames": "au:AdministrativeUnit", "count": "300",
                            "srsName": "urn:ogc:def:crs:EPSG::4326", "bbox": BB84, "outputFormat": "application/geo+json"})
raw = ROOT / "munis.geojson"
curl("https://www.ign.es/wfs-inspire/unidades-administrativas?" + q, raw, 400)
m = gpd.read_file(raw)
m = m[m["nationalLevel"].astype(str).str.contains("4thOrder")].copy()
m["nombre"] = m["name"].map(lambda v: (re.search(r"'text': '([^']+)'", str(v)) or [None, None])[1])
m["codigo"] = m["nationalCode"].astype(str)
m["provincia"] = m["codigo"].str[4:6].map({"05": "Ávila", "28": "Madrid", "45": "Toledo"})
m = m[["codigo", "nombre", "provincia", "geometry"]].sort_values("nombre")
print("municipalities", len(m), m["provincia"].value_counts(dropna=False).to_dict(), "null names", int(m["nombre"].isna().sum()), flush=True)
m.to_file(DATA / "líneas_límite_municipales.geojson", driver="GeoJSON")
mu = m.to_crs(fp.crs)
touched = mu[mu.intersects(fpu)]
print("touched by fire", len(touched), touched["provincia"].value_counts().to_dict(), flush=True)

# Natura 2000
q = urllib.parse.urlencode({"where": "1=1", "geometry": f"{W},{S},{E},{N}", "geometryType": "esriGeometryEnvelope", "inSR": "4326",
                            "spatialRel": "esriSpatialRelIntersects", "outFields": "*", "outSR": "4326", "f": "geojson"})
raw = ROOT / "natura.geojson"
curl("https://bio.discomap.eea.europa.eu/arcgis/rest/services/ProtectedSites/Natura2000Sites/MapServer/2/query?" + q, raw, 300)
n = gpd.read_file(raw)[["SITECODE", "SITENAME", "SITETYPE", "RELEASE_DATE", "Area_km2", "geometry"]]
print("natura", len(n), flush=True)
n.to_file(DATA / "natura2000.geojson", driver="GeoJSON")

# rivers
r = paged("https://servicios.idee.es/wfs-inspire/hidrografia", "hy-p:Watercourse", 500, "rios")
r = r.rename(columns={"text": "nombre", "inspireId|Identifier|localId": "id_ign", "LocalisedCharacterString": "tipo"})
r = r[[c for c in ("id_ign", "nombre", "tipo", "origin", "geometry") if c in r]]
r = r[r.is_valid & ~r.is_empty]
print("rivers", len(r), flush=True)
r.to_file(DATA / "rios.geojson", driver="GeoJSON")

# population centres
t = paged("https://www.ign.es/wfs-inspire/ngbe", "gn:NamedPlace", 1000, "ngbe")
kinds = ["Núcleos de población", "Capital de municipio", "Entidad singular", "Entidad menor de población"]
t = t[t["LocalisedCharacterString"].isin(kinds)].to_crs(fp.crs)
t = t[t.within(aoi)].copy()
t["tipo"] = t["LocalisedCharacterString"]
t["nombre"] = t["text"].map(lambda v: (re.findall(r"'([^']+)'", str(v)) or [str(v)])[0])
t["rank"] = t["tipo"].map({k: i for i, k in enumerate(kinds)})
t["cell"] = (t.geometry.x // 500).astype(int).astype(str) + "_" + (t.geometry.y // 500).astype(int).astype(str)
t = t.sort_values("rank").drop_duplicates(subset=["nombre", "cell"])
print("settlements", len(t), flush=True)
t[["localId", "nombre", "tipo", "geometry"]].to_crs(4326).to_file(DATA / "nucleos_poblacion.geojson", driver="GeoJSON")

# Catastro buildings
links = {}
for pr in ("05", "28", "45"):
    body = curl(f"https://www.catastro.hacienda.gob.es/INSPIRE/buildings/{pr}/ES.SDGC.bu.atom_{pr}.xml")
    for z in re.findall(r'href="([^"]+\.zip)"', body.decode("utf-8", "replace")):
        links[(pr, frozenset(plain(z.split("/")[-2].split("-", 1)[1])))] = z
frames, missing = [], []
for name, prov in zip(touched["nombre"], touched["codigo"].str[4:6]):
    z = links.get((prov, frozenset(plain(name))))
    if z is None:
        missing.append(name)
        print("NO CATASTRO MATCH", name, prov, flush=True)
        continue
    d = ROOT / "catastro" / z.split("/")[-1].removesuffix(".zip")
    d.mkdir(parents=True, exist_ok=True)
    zp = d / "data.zip"
    for attempt in range(4):
        curl(z, zp)
        try:
            with zipfile.ZipFile(zp) as zf:
                member = [x for x in zf.namelist() if x.lower().endswith(".building.gml")][0]
                (d / "building.gml").write_bytes(zf.read(member))
            break
        except Exception as e:
            print(" retry", name, type(e).__name__, flush=True)
            time.sleep(10)
    g = gpd.read_file(d / "building.gml").to_crs(fp.crs)
    g["municipio"] = name
    frames.append(g)
    print(" ", name, len(g), flush=True)
b = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=fp.crs)
b = b[b.intersects(aoi)].copy()
b["numberOfDwellings"] = pd.to_numeric(b["numberOfDwellings"], errors="coerce").fillna(0).astype(int)
b = b[["localId", "municipio", "currentUse", "conditionOfConstruction", "numberOfDwellings", "beginning", "geometry"]]
b["beginning"] = b["beginning"].astype(str)
out = DATA / "edificios_catastro.gpkg"
if out.exists():
    out.unlink()
b.to_file(out, layer="edificios", driver="GPKG")
print("buildings", len(b), "missing municipalities", missing, flush=True)
print("DONE", flush=True)
