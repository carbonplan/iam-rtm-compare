# ---------------------------------------------------
# 
# make batch file for a set of generic experimental
# conditions. 
# 
# Define the experimental conditions and the sites, 
# and simulate *some* combinations of inputs across 
# ( some inputs are explored one-at-a-time, and 
#   others, like dustrate, are explored for all 
#   values across all cases )
# 
# T Kukla ; nov 2025
# 
# ---------------------------------------------------
# %% 
import itertools
import os
import sys

import numpy as np
import pandas as pd

# --- read in batch helper functions  
from ew_workflows import batch_helperFxns as bhf
# ---

# %%
#  ----------------------------------------------------------------
# [ RUN NAMES ]
# exe_central1p5 ; exe_central2p0 ; exe_limBio1p5 ; exe_lowCCS1p5 ; 
# exe_lowCDR1p5 ; exe_noBiochar1p5 ; exe_noWeathering1p5 ; 
# exe_ref1p5 ; exe_techBehaviorChange1p5 ; exe_techBehaviorChange2p0 ;
# [ FOCUS CASES ]
# exe_central1p5 ; exe_techBehaviorChange1p5 ; exe_lowCDR1p5 ; exe_noWeathering1p5
# ----------------------------------------------------------------
# [ UNIVERSAL ]
RUNNAME = f"exe_central1p5"
EXP_GROUP = "gcam_opr2_psdfull_10xrough"   # (opr = one(cell)-per-region)
SUFF = ""
EXP_SET = f"{RUNNAME}{SUFF}"
dustsp = "gbas"
runtype = "yearly_ltm"
extra_tag = "v0"
pref = f"{EXP_GROUP}_{EXP_SET}_{runtype}_{extra_tag}"
fn = pref + ".csv"
spinup_tag = "_no-amnt"  # "_no-amnt"

# save preferences
savepath_batch = "s3://carbonplan-carbon-removal/ew-workflows-data/scepter/batch/"
multi_run_split = False   # whether to split the csv into multiple files
max_iters_per_set = 20    # [int] the number of runs per csv (only used if multi_run_split is True)

# %% 
# ==========================================================================
# [ read in site pickle ]
celltypes = "one-cell-per-region"
deg_resolution = 1
dfsites = pd.read_pickle(f's3://carbonplan-carbon-removal/ew-workflows-data/gcam/apprate_calculated/{celltypes}/{RUNNAME}/annual_rockflx_densest-ag-region_{deg_resolution}deg.pkl')
# 
# %% 
# [ CONSTANTS ]
# ==========================================================================
# common adjustments:
const_dict = {
    # --- commonly changed ---
    "duration": 80,    # [yr] duration of run (starts from earliest year)
    "dustsp": dustsp,   
    "dustsp_2nd": "amnt",
    "spinrun": f'390_southeast{spinup_tag}', # f'390_southeast{spinup_tag}',
    "dustrate": 1.,   # OVERWRITTEN BY GCAM-DERIVED DUST FLUXES (but need a non-zero number for identifying controls later)

    # --- MULTI-YEAR SPECIFIC ---
    # "dust_ts_dir": f"s3://carbonplan-carbon-removal/ew-workflows-data/scepter/dust/",
    # "dust_ts_fn": f"{dustsp}_150yr_100on-50off_001.csv", 
    # ---------------------------

    # --- PH-REACT SPECIFIC ----
    # "max_time": 100,    # [yr] total amount of time to simulate
    # "duration": 1,    # [yr] duration of each individual run (starts from earliest year)
    # "tph": 7.0,       # [] target pH for each iteration
    # "phnorm_pw": True,  # metric for target pH (true = porewater pH; false = soil pH)
    # ---------------------------

    # --- CONTROL VALUES ---
    "dustrad": 5,   # [micron]
    "secondary_min_rule": "sld_track",
    "poro_updated": "None", # 0.35,
    # "dustrate_2nd": 0.,    # [g/m2/yr] (=t/ha/yr * 100)

    # --- sometimes changed ---
    "imix": 3,    # mixing style (1=fickian; 2=homogeneous; 3=tilling)
    "cec_adsorption_on": True, 

    # --- surface area / poro / psd rules
    "include_psd_full": True, 
    "include_psd_bulk": False,
    'psdrain_log10_sd': 0.05, # [] log 10 standard deviation for psd
    'psdrain_wt': 1.0,       # [] weight for the psd
    'use_psdrain_datfile': False,  # False means we construct the PSD based on inputs, rather than copy an existing data file
    'poro_iter_field': False,      # [default='false'] porosity iteration
    'poro_evol': False,            # [default='false'] porosity evolves with solid phase dissolution
    'sa_rule1': False,       # [True, False, "spinvalue"] SA decreases as porosity increases
    'sa_rule2': True,       # [True, False, "spinvalue"] SA increases as porosity increases
    'include_roughness_sa': False,  # [True, False] whether to apply roughness calculation from Navarre-Sitchler & Brantley, 2007 
                                   # (roughness factor = (beta / a)^0.33) where beta is particle radius in m and a is BET measurement resolution in m, taken to be 10^-10 
    
    # --- other
    "singlerun_seasonality": True,
    'climatedir': f's3://carbonplan-carbon-removal/ew-workflows-data/scepter/clim/one-cell-per-region_1deg/{EXP_SET}/',   # climate input main directory
    "skip_lab_run": True,
    'spindir': "s3://carbonplan-carbon-removal/SCEPTER/scepter_output/spinups/", 

    # # --- update CEC?
    # 'cec_update_from_spinup': True,   # [True, False] whether to update CEC and alpha vars relative to the spinup value (False means no change to cec.in is made)
    # 'cec': 45,                   # [cmol kg-1] [only used if cec_update_from_spinup == True] cation exchange capacity
    # 'alpha': 3.4,                   # [only used if cec_update_from_spinup == True]
    
    # --- compute specific
    'aws_save': "move",              # ["move", "copy", None] whether to "move" file to aws, just copy it, or nothing at all
    'aws_bucket': f"s3://carbonplan-carbon-removal/SCEPTER/scepter_output/{EXP_GROUP}/{EXP_SET}",  # where to save at AWS (only used if 'aws_save'=True)
    
    # --- which executable to use
    'scepter_exec_name': 'scepter_richards'  # ['scepter', 'scepter_rateA', ...]
}
# ==========================================================================


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# %% 
# [ SIMULATION GROUPS ]
# ==========================================================================
individual_cases = { 
    # [ Each element is one (and only) tweak from ctrl ]
    # "dustrad": [10, 300],    # [micron]
    # "secondary_min_rule": ["remove"],
    # "poro_updated": [0.25, 0.45],    # initial porosity
}
# ==========================================================================


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# %% 
# [ BY SITE ]
# ==========================================================================
# sites = [f'4_central_valley', f'264_cornbelt', f'390_southeast', f'468_northeast']
sites = list(dfsites['region_nospaces'].unique())
climate_timestep = runtype   # ( suffix for dirs here: s3://carbonplan-carbon-removal/ew-workflows-data/scepter/clim/era5_2006-01-01_2020-12-31/ )
by_site = {   # values must have same order as 'sites' var
    # "spinrun": [f'{site}{spinup_tag}' for site in sites],
    "climatefiles": [
        f'{site}_{deg_resolution}_{runtype}' for site in sites
    ]  # these serve as the site name when there is no cliamte file to use
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# %% 
# [ ALL COMBINATIONS OF THESE ]
# ==========================================================================
all_combinations = { # [ include ctrl values (!!) ]
    # "dustrate_2nd": [0, 5, 35], # [g/m2/yr] (=t/ha/yr * 100)
    # "dustrate": [x*100 for x in [0.1, 1, 10]],
}
# ==========================================================================

# %% 
# --- BUILD DATAFRAME
df = bhf.build_df_one_at_a_time(
    pref, 
    const_dict, 
    sites, 
    by_site, 
    all_combinations, 
    individual_cases,
    add_ctrl=True,
    add_index=True,
)

# [ account for psd tracking ]
if const_dict['include_psd_full'] or const_dict['include_psd_bulk']:
    df['psdrain_meanRad'] = df['dustrad'] / 1e6  # (convert um to m)


# %%
# --- SAVE
bhf.save_df(df, savepath_batch, fn, multi_run_split, max_iters_per_set)

# %%
pref
# %%

