"""
ERA5-Land Infiltration Proxy Extractor — wettest site per GCAM region
=======================================================================
Same logic as era5_land_edh_infiltration.py, but instead of sampling the
representative sites from the validation parquet, this script samples the
grid cell with the highest long-term mean P - E in each GCAM region.

Sites are read from the pre-computed CSV:
  s3://carbonplan-carbon-removal/ew-workflows-data/gcam/
      wettest-cell-per-region/wettest_cell_per_gcam_region.csv

Output zarrs are written directly to S3:
  s3://carbonplan-carbon-removal/ew-workflows-data/gcam/
      gcam-era5-wettestSite-zarr/era5_sites_{year}.zarr

NOTE on "infiltration":
  Soil water infiltration is not a native ERA5-Land output variable.
  This script uses surface runoff (sro) and sub-surface runoff (ssro)
  as proxies. If you need a true infiltration estimate, a common
  approximation is:  infiltration ≈ total_precip - surface_runoff
  which can be computed from `tp` and `sro`. Set VARIABLES below
  to match your needs.

Requirements:
    pip install xarray zarr dask aiohttp s3fs pandas

Auth setup (one-time):
    Add the following to ~/.netrc  (create if it doesn't exist):

        machine data.earthdatahub.destine.eu
            login edh
            password <your EDH personal access token>

    Then restrict permissions:
        chmod 600 ~/.netrc

    AWS credentials must also be configured (e.g. via ~/.aws/credentials
    or environment variables) with write access to carbonplan-carbon-removal.
"""
# %%
import os

import coiled
from dask.distributed import Client
import pandas as pd
import xarray as xr
import numpy as np
import s3fs

# ─────────────────────────────────────────────
# 1.  CONFIGURATION
# ─────────────────────────────────────────────

VARIABLES_accum_nonE = ["sro", "ssro", "ro", "tp"]
VARIABLES_accum_E    = ["e", "pev"]
VARIABLES_accum      = VARIABLES_accum_nonE + VARIABLES_accum_E
VARIABLES_other      = ["swvl1", "swvl2", "swvl3", "stl1"]
VARIABLES_all        = VARIABLES_accum_nonE + VARIABLES_accum_E + VARIABLES_other

# *******************************************************
# Cropland filtering — must match the value used in find_wettest_cell_per_gcam_region.py.
# Set to 0.0 to use the unfiltered CSV (original behavior).
# e.g. 0.25 means the wettest-cell CSV was built from cells where cropland fraction ≥ 25 %.
crop_fraction_threshold = 0.
_crop_tag = f"_crop{int(crop_fraction_threshold * 100)}pct" if crop_fraction_threshold > 0 else ""
# *******************************************************

WETTEST_CELL_CSV = (
    "s3://carbonplan-carbon-removal/ew-workflows-data/gcam/"
    f"wettest-cell-per-region/wettest_cell_per_gcam_region{_crop_tag}.csv"
)

# Output: one zarr per year, written directly to S3
OUT_S3_DIR = (
    "s3://carbonplan-carbon-removal/ew-workflows-data/gcam/"
    f"gcam-era5-wettestSite{_crop_tag}-zarr"
)

# Earth Data Hub ERA5-Land Zarr URL
atoken = "edh_pat_5006183c5b85f3a21566968ca1710fe2f600be9e71841e2f750adeb4ee818ae35b2595d1d470d641df0560f317133278"
EDH_URL = (
    f"https://edh:{atoken}@data.earthdatahub.destine.eu"
    "/era5/reanalysis-era5-land-no-antartica-v0.zarr"
)

# %%
# ─────────────────────────────────────────────
# 2.  LOAD WETTEST-CELL SITES
# ─────────────────────────────────────────────
# CSV columns: reg_id, region, region_nospaces, lat, lon, pme_mm_yr

print("Reading wettest-cell CSV ...")
fs = s3fs.S3FileSystem(anon=False)
with fs.open(WETTEST_CELL_CSV) as f:
    sites = pd.read_csv(f)

print(f"  {len(sites)} GCAM regions loaded")
print(sites[["reg_id", "region", "lat", "lon", "pme_mm_yr"]].to_string(index=False))

lats    = xr.DataArray(sites["lat"].values, dims="site")
lons_180 = xr.DataArray(sites["lon"].values, dims="site")
lons    = lons_180 % 360   # EDH uses 0–360 longitude convention

# %%
# ─────────────────────────────────────────────
# 3.  OPEN THE ZARR STORE (lazy — nothing downloaded yet)
# ─────────────────────────────────────────────

print("Opening ERA5-Land Zarr store (lazy) ...")
ds_full = xr.open_dataset(
    EDH_URL,
    storage_options={"client_kwargs": {"trust_env": True}},  # uses ~/.netrc
    chunks={},   # enables dask lazy loading
    engine="zarr",
)
print(f"Dataset variables available: {list(ds_full.data_vars)}")

# %%
# ─────────────────────────────────────────────
# 4.  SUBSET: variables + time + sites
# ─────────────────────────────────────────────

ds = ds_full[VARIABLES_all]

ds = ds.sel(
    latitude=lats,
    longitude=lons,
    method="nearest",
).sel(valid_time=slice("2000-01-01", "2024-12-31"))

# Attach human-readable coordinates
ds = ds.assign_coords(
    site_id=("site", sites["reg_id"].values),
    target_lat=("site", sites["lat"].values),
    target_lon=("site", sites["lon"].values),
)

# ERA5 convention: evaporation is negative (land→atmosphere), flip sign
for var in VARIABLES_accum_E:
    ds[var] = ds[var] * -1

# %%
# ─────────────────────────────────────────────
# 5.  COMPUTE ANNUAL STATS AND SAVE PER YEAR
# ─────────────────────────────────────────────
# Accumulation vars: daily max → annual sum
# State vars: annual mean
# Each year is saved as a separate zarr directly on S3.

EXPECTED_VARS = (
    [f"{v}_annual_sum" for v in VARIABLES_accum]
    + [f"{v}_annual_mean" for v in VARIABLES_other]
)

def year_is_complete(path, s3):
    """Return True if the zarr at path exists on S3 and has all expected variables."""
    try:
        if not s3.exists(path):
            return False
        ds_check = xr.open_zarr(path, storage_options={"anon": False})
        return all(v in ds_check for v in EXPECTED_VARS)
    except Exception:
        return False

cluster = coiled.Cluster(n_workers=4, region="us-east-1")
client = Client(cluster)

for year, ds_year in ds.groupby("valid_time.year"):
    out_path = os.path.join(OUT_S3_DIR, f"era5_sites_{year}.zarr")
    print(f"\nYear {year} ...")

    if year_is_complete(out_path, fs):
        print("  complete on S3 — skipping")
        continue

    exists = fs.exists(out_path)
    if exists:
        print("  incomplete zarr found on S3 — overwriting")

    # accumulation variables: daily max → annual sum
    ds_daily_max = ds_year[VARIABLES_accum].resample(valid_time="1D").max()
    ds_accum = ds_daily_max.sum(dim="valid_time").rename(
        {v: f"{v}_annual_sum" for v in VARIABLES_accum}
    )

    # state variables: annual mean
    ds_mean = ds_year[VARIABLES_other].mean(dim="valid_time").rename(
        {v: f"{v}_annual_mean" for v in VARIABLES_other}
    )

    ds_annual = xr.merge([ds_accum, ds_mean]).expand_dims(year=[year])

    print(f"  Downloading {year} ...")
    ds_annual = ds_annual.compute()

    mode = "w" if exists else "w-"
    ds_annual.to_zarr(out_path, mode=mode, storage_options={"anon": False})
    print(f"  Saved → {out_path}")

client.close()
cluster.close()
print("\nDone.")

# %%
# ─────────────────────────────────────────────
# 6.  OPTIONAL: READ BACK AND COMBINE
# ─────────────────────────────────────────────

# annual_datasets = []
# for year in range(2000, 2025):
#     zarr_path = os.path.join(OUT_S3_DIR, f"era5_sites_{year}.zarr")
#     tmpds = xr.open_zarr(zarr_path, storage_options={"anon": False})
#     annual_datasets.append(tmpds)

# ds_combined = xr.concat(annual_datasets, dim="year")
# print(ds_combined)

# %%
