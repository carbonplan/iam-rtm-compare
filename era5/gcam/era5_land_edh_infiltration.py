"""
ERA5-Land Infiltration Proxy Extractor — Earth Data Hub
=========================================================
Streams ERA5-Land data directly from the Earth Data Hub Zarr store,
extracts annual means for 32 specific grid cells, and saves as
NetCDF + Zarr. No file downloads required — only the data you need
is transferred.

NOTE on "infiltration":
  Soil water infiltration is not a native ERA5-Land output variable.
  This script uses surface runoff (sro) and sub-surface runoff (ssro)
  as proxies. If you need a true infiltration estimate, a common
  approximation is:  infiltration ≈ total_precip - surface_runoff
  which can be computed from `tp` and `sro`. Set VARIABLES below
  to match your needs.

Requirements:
    pip install xarray zarr dask aiohttp netcdf4

Auth setup (one-time):
    Add the following to ~/.netrc  (create if it doesn't exist):

        machine data.earthdatahub.destine.eu
            login edh
            password <your EDH personal access token>

    Then restrict permissions:
        chmod 600 ~/.netrc
"""
# %%
import os 

import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path

# ─────────────────────────────────────────────
# 1.  CONFIGURATION
# ─────────────────────────────────────────────

# Variables to extract. Options relevant to infiltration:
#   "sro"   - Surface runoff (m)
#   "ssro"  - Sub-surface runoff (m)
#   "ro"    - Total runoff (m)
#   "tp"    - Total precipitation (m)  -- useful for infiltration = tp - sro
# Add/remove as needed.
VARIABLES_accum_nonE= ["sro", "ssro", "ro", "tp"]
VARIABLES_accum_E = ["e", "pev"]
VARIABLES_accum = VARIABLES_accum_nonE + VARIABLES_accum_E
VARIABLES_other = ["swvl1", "swvl2", "swvl3"]
VARIABLES_all = VARIABLES_accum_nonE + VARIABLES_accum_E + VARIABLES_other

# Output directory
OUT_DIR = Path("era5_infiltration_output")
OUT_DIR.mkdir(exist_ok=True)

# Earth Data Hub ERA5-Land Zarr URL
atoken = "edh_pat_5006183c5b85f3a21566968ca1710fe2f600be9e71841e2f750adeb4ee818ae35b2595d1d470d641df0560f317133278"
EDH_URL = (
    f"https://edh:{atoken}@data.earthdatahub.destine.eu"
    "/era5/reanalysis-era5-land-no-antartica-v0.zarr"
)

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
# 4.  SUBSET: variables + years + bounding box
#     This defines what will actually be downloaded.
# ─────────────────────────────────────────────
import geopandas as gpd
# load the sites we want to look at 
spin_parquet = 's3://carbonplan-carbon-removal/ew-workflows-data/valdation-data/sampled_locations_updated.parquet'
rtgdf = gpd.read_parquet(spin_parquet)
lats = xr.DataArray(rtgdf.geometry.y.values, dims="site")
lons_180 = xr.DataArray(rtgdf.geometry.x.values, dims="site")
# convert lons to 0-360 
lons = lons_180 % 360

# Keep only the variables we need
ds = ds_full[VARIABLES_all]

# EDH uses 'latitude'/'longitude' dimension names
ds = ds.sel(
    latitude=lats,
    longitude=lons,
    method="nearest",
).sel(valid_time = slice("2000-01-01", "2024-12-31"))

# Attach human-readable coordinates
ds = ds.assign_coords(
    site_id=("site", rtgdf["reg_id"].values),
    target_lat=("site", rtgdf.geometry.y.values),
    target_lon=("site", rtgdf.geometry.x.values),
)

# multiply E by -1 (makes E positive)
for var in VARIABLES_accum_E:
    ds[var] = ds[var] * -1

# %% 
# --- annual means
last_complete_year = 2007
annual_datasets = []
OUT_DIR = "/Users/tylerkukla/Downloads/tmp-era5"
for year, ds_year in ds.groupby("valid_time.year"):
    print(f"Computing annual stats for {year} ...")

    if year <= last_complete_year:
        print("already saved... moving to next year ")
        continue

    # --- accumulation variables: daily max → annual sum
    ds_daily_max = ds_year[VARIABLES_accum].resample(valid_time="1D").max()
    ds_accum = ds_daily_max.sum(dim="valid_time").rename(
        {v: f"{v}_annual_sum" for v in VARIABLES_accum}
    )

    # --- other variables: annual mean
    ds_mean = ds_year[VARIABLES_other].mean(dim="valid_time").rename(
        {v: f"{v}_annual_mean" for v in VARIABLES_other}
    )

    # --- merge and tag with year
    ds_annual = xr.merge([ds_accum, ds_mean]).expand_dims(year=[year])
    annual_datasets.append(ds_annual)

    # --- compute and save each year individually
    print(f"  Downloading {year} ...")
    ds_annual = ds_annual.compute()
    # save 
    out_path = os.path.join(OUT_DIR, f"era5_sites_{year}.zarr")
    ds_annual.to_zarr(out_path)
    print(f"  Saved {out_path}")

# %%
# --- read in data to combine
annual_datasets = []
OUT_DIR = "/Users/tylerkukla/Downloads/tmp-era5"
for year, ds_year in ds.groupby("valid_time.year"):
   tmpds = xr.open_zarr(os.path.join(OUT_DIR, f"era5_sites_{year}.zarr"))
   annual_datasets.append(tmpds)

ds_combined = xr.concat(annual_datasets, dim="year")

# %% 
# --- 
# # Trigger computation (this is where data actually downloads)
# print("Downloading selected data chunks (this may take a moment) ...")
# ds_combined = ds_combined.compute()
# print("Download complete.")
# print(ds_combined)

