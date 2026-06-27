"""
find_wettest_cell_per_gcam_region.py
--------------------------------------
For each of GCAM's 32 regions, find the grid cell with the highest
long-term mean P - E (precipitation minus evaporation) from ERA5.

ERA5 convention: evaporation is negative (land→atmosphere flux),
so  P - E  =  total_precipitation + evaporation.

Approach:
  1. Compute long-term mean P - E on the native ERA5 global grid.
  2. Assign each ERA5 grid-cell center to a GCAM region via a
     point-in-polygon spatial join against the GCAM shapefile.
  3. Return the highest-P-E cell in each region.

Output: CSV with one row per region — reg_id, region name, lat, lon,
        and long-term mean annual P - E in mm/yr.
"""
# %%
import coiled
import icechunk
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from dask.distributed import Client

# ──────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ──────────────────────────────────────────────────────────────────────
# *******************************************************
# Cropland filtering — set to 0.0 to disable (no filtering, original behavior).
# e.g. 0.10 means only consider ERA5 cells where cropland fraction ≥ 10 %.
crop_fraction_threshold = 0.
cropland_zarr_path = (
    "s3://carbonplan-carbon-removal/ew-workflows-data/gcam/"
    "cropland+gcamRegions/ds_cropland+gcam_1deg.zarr/"
)

_crop_tag = f"_crop{int(crop_fraction_threshold * 100)}pct" if crop_fraction_threshold > 0 else ""
save_path = (
    "s3://carbonplan-carbon-removal/ew-workflows-data/gcam/"
    f"wettest-cell-per-region/wettest_cell_per_gcam_region{_crop_tag}.csv"
)
# *******************************************************

year_range = [2000, 2024]

era5_bucket = "carbonplan-carbon-removal"
era5_dir    = "era5/preprocessed_icechunk"

gcam_shapefile = (
    "s3://carbonplan-carbon-removal/ew-workflows-data/gcam/"
    "gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/main_outputs/"
    "region_boundaries_moirai_landcells_3p1_0p5arcmin.shp"
)


# %%
# ──────────────────────────────────────────────────────────────────────
# 0.  Coiled cluster — same region as the S3 bucket so data stays local
# ──────────────────────────────────────────────────────────────────────
cluster = coiled.Cluster(
    n_workers=20,
    region="us-west-2",           # match carbonplan-carbon-removal bucket
    package_sync=True,            # mirrors your local env on workers
)
client = Client(cluster)
print(client.dashboard_link)

# ──────────────────────────────────────────────────────────────────────
# 1.  ERA5 full-grid data via icechunk  (lazy, ~859 GB)
# ──────────────────────────────────────────────────────────────────────
print("Opening ERA5 icechunk store (lazy) ...")
storage = icechunk.s3_storage(bucket=era5_bucket, prefix=era5_dir, from_env=True)
repo    = icechunk.Repository.open(storage)
session = repo.readonly_session("main")
ds_era5 = xr.open_zarr(session.store, consolidated=False)

# subset to year range
ds_era5 = ds_era5.sel(
    time=slice(f"{year_range[0]}-01-01", f"{year_range[1]}-12-31")
)

# %%
# ──────────────────────────────────────────────────────────────────────
# 2.  Compute long-term mean P - E on the native ERA5 grid
#     ERA5 evaporation is negative (land→atm), so P - E = tp + e
#     Time-mean preserves the relative ranking across cells.
# ──────────────────────────────────────────────────────────────────────
print("Computing long-term mean P - E (triggers ERA5 download) ...")
pme = (
    ds_era5["total_precipitation"].mean("time")
    + ds_era5["evaporation"].mean("time")
).compute()
# pme: dims (latitude, longitude); units: m hr-1 (mean hourly accumulation)

# %%
# ──────────────────────────────────────────────────────────────────────
# 3.  Build a flat dataframe of ERA5 grid-cell centers + P-E values
# ──────────────────────────────────────────────────────────────────────
lats_era5 = pme.latitude.values   # 1-D (n_lat,)
lons_era5 = pme.longitude.values  # 1-D (n_lon,), -179.6 to 179.4

lat_grid, lon_grid = np.meshgrid(lats_era5, lons_era5, indexing="ij")

df_cells = pd.DataFrame({
    "lat": lat_grid.flatten(),
    "lon": lon_grid.flatten(),
    "pme": pme.values.flatten(),
})

# %%
# ──────────────────────────────────────────────────────────────────────
# 4.  Assign each ERA5 cell to a GCAM region via spatial join
#     The shapefile defines the 32 GCAM land regions; ocean cells
#     won't match any polygon and are dropped automatically.
# ──────────────────────────────────────────────────────────────────────
print("Assigning ERA5 cells to GCAM regions (spatial join) ...")
gdf_regions = gpd.read_file(gcam_shapefile)[["reg_id", "reg_nm", "geometry"]]
gdf_regions = gdf_regions.to_crs("EPSG:4326")
gdf_regions = gdf_regions.dissolve(by="reg_id").reset_index()   # one row per region

gdf_cells = gpd.GeoDataFrame(
    df_cells,
    geometry=gpd.points_from_xy(df_cells["lon"], df_cells["lat"]),
    crs="EPSG:4326",
)

gdf_joined = gdf_cells.sjoin(
    gdf_regions[["reg_id", "reg_nm", "geometry"]],
    how="left",
    predicate="within",
)

# drop ocean / unmatched cells
gdf_land = gdf_joined.dropna(subset=["reg_id"]).copy()
gdf_land["reg_id"] = gdf_land["reg_id"].astype(int)

# %%
# ──────────────────────────────────────────────────────────────────────
# 4b. Load cropland fractions onto all land cells (needed for fallback)
# ──────────────────────────────────────────────────────────────────────
if crop_fraction_threshold > 0:
    print("Loading cropland fractions ...")
    da_crop = xr.open_zarr(cropland_zarr_path)["cropland_frac"].rename({"y": "latitude", "x": "longitude"})
    da_crop_regrid = da_crop.interp(
        latitude=pme.latitude, longitude=pme.longitude, method="nearest"
    ).compute()
    crop_vals = pd.DataFrame({
        "lat":       lat_grid.flatten(),
        "lon":       lon_grid.flatten(),
        "crop_frac": da_crop_regrid.values.flatten(),
    })
    gdf_land = gdf_land.merge(crop_vals, on=["lat", "lon"], how="left")
    gdf_land["crop_frac"] = gdf_land["crop_frac"].fillna(0.0)
    n_crop_pass = (gdf_land["crop_frac"] >= crop_fraction_threshold).sum()
    print(f"  {n_crop_pass:,} / {len(gdf_land):,} land cells have cropland fraction ≥ {crop_fraction_threshold:.0%}")

# %%
# ──────────────────────────────────────────────────────────────────────
# 5.  Find the highest-P-E cell per region, with cropland threshold fallback.
#     For each region, start at crop_fraction_threshold and step down by 5%
#     until the wettest qualifying cell reaches PME_MIN_MM_YR.
#     If even threshold=0 doesn't reach it, use that cell's location but
#     clamp pme to PME_MIN_MM_YR.
# ──────────────────────────────────────────────────────────────────────
PME_MIN_MM_YR = 10.0
PME_MIN_RAW   = PME_MIN_MM_YR / (1000 * 8760)

# Descending threshold sequence: crop_fraction_threshold → 0 in 5% steps
if crop_fraction_threshold > 0:
    n_steps = int(round(crop_fraction_threshold / 0.05)) + 1
    threshold_steps = [round(crop_fraction_threshold - i * 0.05, 2) for i in range(n_steps)]
    if threshold_steps[-1] != 0.0:
        threshold_steps.append(0.0)
else:
    threshold_steps = [0.0]

records = []
for reg_id, grp in gdf_land.groupby("reg_id"):
    chosen_row  = None
    used_thresh = None

    for thresh in threshold_steps:
        candidates = grp[grp["crop_frac"] >= thresh] if thresh > 0 else grp
        if candidates.empty:
            continue
        idx = candidates["pme"].idxmax()
        row = candidates.loc[idx]
        chosen_row  = row    # always update so we end with the best-at-threshold-0 result
        used_thresh = thresh
        if row["pme"] >= PME_MIN_RAW:
            break            # criterion met — stop lowering threshold

    if chosen_row is None:
        continue

    clamped = chosen_row["pme"] < PME_MIN_RAW
    records.append({
        "reg_id":              int(chosen_row["reg_id"]),
        "reg_nm":              chosen_row["reg_nm"],
        "lat":                 chosen_row["lat"],
        "lon":                 chosen_row["lon"],
        "pme":                 PME_MIN_RAW if clamped else float(chosen_row["pme"]),
        "used_crop_threshold": used_thresh,
    })
    if clamped:
        print(f"  WARNING: {chosen_row['reg_nm']!r} — best P-E is "
              f"{chosen_row['pme']*1000*8760:.1f} mm/yr even at crop≥{used_thresh:.0%}; "
              f"clamped to {PME_MIN_MM_YR:.0f} mm/yr")

df_wettest = pd.DataFrame(records).rename(columns={"reg_nm": "region"})

# ──────────────────────────────────────────────────────────────────────
# 6.  Convert units and tidy
#     mean [m hr-1] × 1000 [mm m-1] × 8760 [hr yr-1] = mm yr-1
# ──────────────────────────────────────────────────────────────────────
df_wettest["region_nospaces"] = df_wettest["region"].str.replace(" ", "_", regex=False)
df_wettest["pme_mm_yr"]       = df_wettest["pme"] * 1000 * 8760

df_wettest = (
    df_wettest[["reg_id", "region", "region_nospaces", "lat", "lon", "pme_mm_yr", "used_crop_threshold"]]
    .sort_values("reg_id")
    .reset_index(drop=True)
)

# %%
# ──────────────────────────────────────────────────────────────────────
# 7.  Print and save
# ──────────────────────────────────────────────────────────────────────
print("\nWettest cell per GCAM region (highest P - E):")
print(df_wettest.to_string(index=False))

df_wettest.to_csv(save_path, index=False)
print(f"\nSaved: {save_path}")

client.close()
cluster.close()

# %%
