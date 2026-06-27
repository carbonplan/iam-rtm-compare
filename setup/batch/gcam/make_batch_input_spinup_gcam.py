# ---------------------------------------------
# 
# generate spinup inputs for gcam runs
# # 
# data validation: https://github.com/carbonplan/SCEPTER-data-validation
# # 
# soilGrids docs: https://docs.isric.org/globaldata/soilgrids/SoilGrids_faqs_01.html
# 
# ---------------------------------------------
# %% 
import os

import geopandas as gpd
import numpy as np
import pandas as pd 
import xarray as xr
# [ for optional plotting at bottom of script ]
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import matplotlib.pyplot as plt

# --- read in batch helper functions  
from ew_workflows import batch_helperFxns as bhf
# ---

# %%
# *************************************************************** #
# --- set version (tracks soil/model inputs; independent of site strategy)
this_version = "v1"           # PmE runoff; pw pH 15-30cm;
# this_version = "v2"             # as v1 but max tPh = 7.0;
# this_version = "v3"           # as v1 but soilCO2 is doubled;
# this_version = "v1colombia"   # all sites have colombia spinup inputs

# --- set site strategy
# SITE_STRATEGY = "largest_cropland"      # original densest-ag-region sites
# SITE_STRATEGY = "wettestSite"           # wettest P-E cell per GCAM region
SITE_STRATEGY = "wettestSite_crop25pct"  # wettest P-E cell with ≥25% cropland

SITE_DATASETS = {
    "largest_cropland": {
        "parquet": "s3://carbonplan-carbon-removal/ew-workflows-data/valdation-data/sampled_locations_updated.parquet",
        "zarr_dir": "s3://carbonplan-carbon-removal/ew-workflows-data/gcam/gcam-era5-zarr/",
    },
    "wettestSite": {
        "parquet": "s3://carbonplan-carbon-removal/ew-workflows-data/valdation-data/sampled_locations_wettestSite.parquet",
        "zarr_dir": "s3://carbonplan-carbon-removal/ew-workflows-data/gcam/gcam-era5-wettestSite-zarr/",
    },
    "wettestSite_crop25pct": {
        "parquet": "s3://carbonplan-carbon-removal/ew-workflows-data/valdation-data/sampled_locations_wettestSite_crop25pct.parquet",
        "zarr_dir": "s3://carbonplan-carbon-removal/ew-workflows-data/gcam/gcam-era5-wettestSite_crop25pct-zarr/",
    },
}
# *************************************************************** #
#
# --- read in parquet data
spin_parquet = SITE_DATASETS[SITE_STRATEGY]["parquet"]
rtgdf = gpd.read_parquet(spin_parquet)
# rtgdf.columns

# --- read in zarr
datpath = SITE_DATASETS[SITE_STRATEGY]["zarr_dir"]
year_range = [2000, 2024]
ann_datasets = []
for year in range(year_range[0], year_range[1]+1):
    tds = xr.open_zarr(os.path.join(datpath, f"era5_sites_{year}.zarr"))
    ann_datasets.append(tds)
ds = xr.concat(ann_datasets, dim="year")

# %% 
# --- compute target soilCO2
# using Gwiazda and Broecker 1994 method as adopted by Godderis et al. 2010 (Time evolution of the mineralogical composition of Mississippi Valley loess over the last 10 kyr);
# Gaillardet et al. 2019 (Global climate control on carbonate weathering intensity); and Zeng et al. 2022 (Natural and anthropogenic driving forces of carbonate weathering and the related carbon sink flux: a model comparison study at global scale)
# 
# Note NPP appears to be in KgC/m2/yr / 0.0001 (scale factor) ; equation requires gC/m2/yr
# (NPP user guide: https://www.umt.edu/numerical-terradynamic-simulation-group/project/modis/user-guides/mod17c61usersguidev11mar112021.pdf)
# 
# [1] convert NPP to gC/m2/yr
rtgdf['MODIS_NPP_gC_m2_yr'] = rtgdf['MODIS_NPP'] * 0.0001 * 1e3
# [2] define atmospheric CO2 ppmv
co2_atm_ppmv = 420
# [3] calculate soil pco2 (eq. 22 of Gaillardet et al. 2019; https://doi.org/10.1016/j.chemgeo.2018.05.009)
A_CONST = 1.03e6    # unit conversion constant
fractNPP = 0.75     # [] fraction of NPP set to soil CO2 prod in deepest roots
rtgdf['soil_pco2_GB94'] = co2_atm_ppmv + ((A_CONST * fractNPP * rtgdf['MODIS_NPP_gC_m2_yr']) / (rtgdf['ERA5_stl1']**2))


# %% 
# [ OPTIONAL : Explore inputs ]
# # --- compare like vars
# rtgdf['diff_ph_0-20'] = rtgdf['SoilGrids_phh2o_15-30cm'] - rtgdf['ISRIC_PHAQ_40-60cm']*10
# rtgdf['diff_cec_0-20'] = rtgdf['SoilGrids_cec_5-15cm'] - rtgdf['ISRIC_ECEC_0-20cm']*10
# rtgdf['diff_soc_0-20'] = rtgdf['SoilGrids_soc_5-15cm']/10 - rtgdf['ISRIC_ORGC_0-20cm']

# pcol = "diff_ph_0-20"
# pcol = "diff_cec_0-20"
# pcol = 'diff_soc_0-20'

# fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
# ax.add_feature(cfeature.COASTLINE)
# rtgdf.plot(column=pcol, legend=True, ax=ax, transform=ccrs.PlateCarree(), markersize=200, cmap=cmocean.cm.balance, vmin=-20, vmax=20)
# plt.show() 
# plt.close()



# %% 
# columns we care about: 
# PORO: 'soil_porosity_GLDAS' 
# pH: 'SoilGrids_phh2o_0-5cm' ; 'SoilGrids_phh2o_5-15cm' ;  'SoilGrids_phh2o_15-30cm'
# CEC: 'SoilGrids_cec_0-5cm' ; 'SoilGrids_cec_5-15cm' ; 'SoilGrids_cec_15-30cm'
# TEMP: 'ERA5_stl1'
# pcol = "ISRIC_ORGC_0-20cm"
# pcol = "SoilGrids_soc_0-5cm"
# pcol = "soil_pco2_GB94"

# fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
# ax.add_feature(cfeature.COASTLINE)
# rtgdf.plot(column=pcol, legend=True, ax=ax, transform=ccrs.PlateCarree(), markersize=200,vmin=100, vmax=10e3)
# plt.show() 
# plt.close()

# %% 
# --- add era5 data to dataframe
# [1] calculate infiltration (need units = [m yr-1])
dfqrun = (ds['tp_annual_sum'].mean('year') - ds['e_annual_sum'].mean('year')).to_dataframe(name='qrun').reset_index()
dfqrun = dfqrun.rename(columns={'site_id': 'reg_id'})
# [2] calculate soil moisture (need units [m3 m-3])
#.    weighted mean (layer1=0-7; layer2=7-28; layer3=28-100)
dfsoilmoisture = ((ds['swvl1_annual_mean'].mean('year')*7 + ds['swvl2_annual_mean'].mean('year')*21 + ds['swvl3_annual_mean'].mean('year')*72)/100).to_dataframe(name='soilmoisture').reset_index()
dfsoilmoisture = dfsoilmoisture.rename(columns={'site_id': 'reg_id'})

# merge
rtgdf = rtgdf.merge(dfqrun[['reg_id', 'qrun']], on='reg_id', how='left')
rtgdf = rtgdf.merge(dfsoilmoisture[['reg_id', 'soilmoisture']], on='reg_id', how='left')


# %% 
# --- create a new dataframe with just the vars we want 
dfx = rtgdf[['reg_id', 'region', 'region_nospaces']].copy()

dfx = dfx.assign(
    # [1] target soil organic matter (varname = 'tsom'; need units = [wt%]) 
    # (note, taking weighted mean by depth; weights sum to one; divide by 10 to convert from dg/kg to g/kg, then 10 more to convert g/kg to wt %)
    tsom=(
        rtgdf['SoilGrids_soc_0-5cm'] * (1/6)
        + rtgdf['SoilGrids_soc_5-15cm'] * (2/6)
        + rtgdf['SoilGrids_soc_15-30cm'] * (3/6)
    ) / 100,

    # [2] target porewater pH (varname = 'tph'; need units = [])
    tph=rtgdf['SoilGrids_phh2o_15-30cm']/10,

    # [3] target exchangeable acidity (varname = 'tec'; need units = %CEC)
    tec=101 - rtgdf['ISRIC_BSAT_0-20cm'],    # 101 to avoid negative values (EC is 1-baseSat)

    # [4] target soil CO2 (varname = 'tsoilco2'; need units = log10 atm)
    tsoilco2=np.log10(rtgdf['soil_pco2_GB94'] / 1e6),   # convert ppmv to log10 atm

    # [5] mean annual temp (varname = 'mat'; need units = degC)
    mat=rtgdf['ERA5_stl1'] - 273.15,

    # [6] porosity (varname = 'poro'; need units = [m3 m-3])
    poro=rtgdf['soil_porosity_GLDAS'],

    # [7] CEC (varname = 'cec'; need units = cmol kg-1)
    # (note, assumes CEC7 is approximately the total potential CEC -- soilGrids CEC is CEC7)
    cec=(   # mean of two depths ; divide by 10 to convert units (note ISRIC CEC data do not require the divide by 10 conversion)
        (rtgdf['SoilGrids_cec_5-15cm'] + rtgdf['SoilGrids_cec_15-30cm']) / 2
    ) / 10,

    # [8] runoff (varname = 'qrun'; need units [m yr-1])
    qrun=rtgdf['qrun'],

    # [9] soil moisture (varname = 'soilmoisture'; need units in [m3 m-3])
    soilmoisture=rtgdf['soilmoisture'],
)

# %% 
# --- [version specific updates]
if this_version == "v2":
    dfx['tph'] = dfx['tph'].clip(upper=7.0)
if this_version == "v3": 
    dfx['tsoilco2'] = dfx['tsoilco2'] + np.log10(2)   # add log10(2) to double soil CO2
# --- add constants and non-numerics
dfx['aws_bucket'] = f"s3://carbonplan-carbon-removal/SCEPTER/scepter_output/spinups/gcam/{this_version}_{SITE_STRATEGY}/"
dfx['aws_save'] = "move"  # ["move", "copy", None] whether to "move" file to aws, just copy it, or nothing at all
dfx['tunespin_timeout'] = 60*60*48   # [s] time before we give up on a tunespin run (7200 is 120*60 = 2 hours)
dfx['scepter_exec_name'] = "scepter" # (note scepter_richards maybe not yet working w/ tunespin script --- need to update frame.in)
dfx['skip_lab_run'] = True  # whether to skip the lab run (True means we skip it and just do the field tunespin)
dfx['spinname'] = (
                    "spinup_gcam_"
                    + rtgdf["region_nospaces"].astype(str)
                    + "_id_"
                    + rtgdf["reg_id"].astype(str)
                )

# %% 
# --- save result -----------------------------------------------------------------------
savepath_batch = "s3://carbonplan-carbon-removal/ew-workflows-data/scepter/batch/"
fn = f"spinup_gcam_{this_version}_{SITE_STRATEGY}.csv"
multi_run_split = False   # whether to split the csv into multiple files
max_iters_per_set = 20    # [int] the number of runs per csv (only used if multi_run_split is True)
bhf.save_df(dfx, savepath_batch, fn, multi_run_split, max_iters_per_set)
# ---------------------------------------------------------------------------------------
print(fn)
# %%
