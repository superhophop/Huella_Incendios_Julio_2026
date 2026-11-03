"""Bottom-up smoke and CO estimate: emission = area x fuel consumed x emission factor
(Seiler & Crutzen 1980). Fuel: van Leeuwen et al. 2014, Table 3. Factors: Andreae 2019, Table 1.
Range: 5th-95th percentile from sampling each mean and SD as a lognormal.

Run after notebook 05, from the notebooks folder:  python estimate_emissions.py
Writes emissions_by_cover.json and emissions_total.json to the output folder.
These are order-of-magnitude estimates, not measurements."""
import json

import numpy as np

from common import OUT, load_json

impact = load_json("step5_impact.json", "05")
D = OUT
area = {c["class"]: c["burned_km2"] for c in impact["landcover"]["classes"]}
# cover: (biome used, FC mean, FC sd, CO mean, CO sd, PM2.5 mean, PM2.5 sd, TPM mean, TPM sd)
SAV = (69, 20, 6.7, 3.3, 8.7, 3.1)
TF = (113, 50, 18.5, 14.4, 18.4, 8.3)
AG = (76, 55, 8.2, 4.4, 12.9, 7.2)
P = {
    "Tree cover": ("Temperate forest", 58, 72, *TF),
    "Grassland": ("Grassland savanna", 4.3, 2.2, *SAV),
    "Shrubland": ("Chaparral (fuel), savanna and grassland (factors)", 27, 19, *SAV),
    "Cropland": ("Crop residue (fuel), agricultural residues (factors)", 6.5, 9.0, *AG),
}
rng = np.random.default_rng(2026)
N = 200_000


def logn(mean, sd):
    s2 = np.log(1 + (sd / mean) ** 2)
    return rng.lognormal(np.log(mean) - s2 / 2, np.sqrt(s2), N)


totals = {"co": 0, "pm25": 0, "tpm": 0, "dm": 0}
rows = []
for cover, (biome, fc, fcs, co, cos, pm, pms, tp, tps) in P.items():
    ha = area[cover] * 100
    fc_s = logn(fc, fcs)
    dm = ha * fc_s
    # same fuel draw for all species in a cover, tonnes
    s = {"dm": dm, "co": dm * logn(co, cos) / 1000, "pm25": dm * logn(pm, pms) / 1000, "tpm": dm * logn(tp, tps) / 1000}
    for k in totals:
        totals[k] = totals[k] + s[k]
    rows.append({"land_cover": cover, "burned_km2": area[cover], "biome_used": biome, "fuel_t_ha": fc, "fuel_sd": fcs,
                 "ef_co_g_kg": co, "ef_pm25_g_kg": pm, "ef_tpm_g_kg": tp,
                 "fuel_burned_t": round(ha * fc), "co_t": round(ha * fc * co / 1000),
                 "pm25_t": round(ha * fc * pm / 1000), "tpm_t": round(ha * fc * tp / 1000)})
tot = {k: sum(r[c] for r in rows) for k, c in (("co", "co_t"), ("pm25", "pm25_t"), ("tpm", "tpm_t"), ("dm", "fuel_burned_t"))}
summ = []
for k, name, unit in (("co", "Carbon monoxide (CO)", "t"), ("pm25", "Fine smoke particles (PM2.5)", "t"),
                      ("tpm", "Total smoke particles (TPM)", "t"), ("dm", "Vegetation burned (dry matter)", "t")):
    lo, hi = np.percentile(totals[k], [5, 95])
    summ.append({"id": k, "species": name, "central_t": int(tot[k]), "low_t": int(round(lo, -2)), "high_t": int(round(hi, -2))})
print(json.dumps(rows, indent=0)[:1800])
print(json.dumps(summ))
print("uncounted km2", {k: v for k, v in area.items() if k not in P})
print("tree share of CO", round(rows[0]["co_t"] / tot["co"], 3))
(D / "emissions_by_cover.json").write_text(json.dumps({"rows": rows}, ensure_ascii=False), encoding="utf-8")
(D / "emissions_total.json").write_text(json.dumps({"rows": summ}, ensure_ascii=False), encoding="utf-8")
