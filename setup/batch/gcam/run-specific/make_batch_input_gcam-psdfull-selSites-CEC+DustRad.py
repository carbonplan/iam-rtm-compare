# ---------------------------------------------------
#
#      !! BUILDS CSV FOR "2AXIS" RUNS !!
#
# make batch file for a selected subset of sites,
# varying over all combinations of cec_factors and
# dustrads (full factorial grid).
#
# Apprate expansion is optional: set repeat_clim_for_apprates
# to True and apprates to a single-element array to use one
# apprate without breaking the climatedir/aws_bucket logic.
#
# Site selection mirrors the DGSA script: list regions
# in SELECTED_REGIONS; the script warns on mismatches.
#
# T Kukla ; july 2026
#
# ---------------------------------------------------
# %%
import numpy as np
import pandas as pd

# --- read in batch helper functions
from ew_workflows import batch_helperFxns as bhf
# ---

# %%
#
# ---------------------------------------------------
# [ test cases ]
# exe_central1p5 ; exe_noWeathering1p5 ; exe_techBehaviorChange1p5
## ---------------------------------------------------
# [ UNIVERSAL ]
RUNNAME = "exe_central1p5"
EXP_GROUP = "gcam_opr_cecdrapp1"   # !! UPDATE AS NEEDED !!
EXP_SET = f"{RUNNAME}"
CLIM_SELECTION_LABEL = "repcells1"  # see labels in `create_clim_inputs_one-cell-per-region.ipynb`
dustsp = "gbas"
spinversion = "v4_repcells1"
apprate_suffix = "-tHaYr-appRate-FIXED"
runtype = "yearly_ltm"
extra_tag = spinversion
pref = f"{EXP_GROUP}_{EXP_SET}_{dustsp}_{runtype}_{extra_tag}"
fn = pref + ".csv"

# save preferences
savepath_batch = "s3://carbonplan-carbon-removal/ew-workflows-data/scepter/batch/"
multi_run_split = False
max_iters_per_set = 20


# %%
# ====================================================================
# [ SITE SELECTION ]
# Specify the region names matching 'region_nospaces' in dfsites.
# The script will warn if any name is not found in dfsites.
# ====================================================================
SELECTED_REGIONS = [
    "China",
    "India",
    "USA",
    "Southeast_Asia",
    "Brazil",
    "EU-15",
    "South_Asia",
]
# ====================================================================


# %%
# ====================================================================
# [ PARAMETER ARRAYS ]
# Full factorial grid: every combination of cec_factors × dustrads.
# Set apprates to a single-element array to skip multi-apprate expansion.
# ====================================================================
cec_factors = np.array([0.1,1,3]) # np.linspace(0.1,4,10)    # [] multipliers on spinup CEC

dustrads    = np.array([0.25, 0.5, 1, 10, 30]) # np.linspace(0.25,100,10)     # [micron] dust particle radii (specified directly,
                                           # 0.25-100 captures roughness 1-30 and dustrad 10-100)

# apprate config — RULES: values must come from
# [0.05, 0.1, 0.25, 0.5, 0.75, 1., 2.5, 5., 7.5, 10., 15., 20., 30., 40., 50.]
# and must have decimals to match file naming conventions.
repeat_clim_for_apprates = True
apprates = np.array([0.05, 0.1, 0.5, 1, 2.5])  # np.array([30.0])        # [t/ha/yr]
# ====================================================================


# %%
# ====================================================================
# [ CONSTANTS ]
# ====================================================================
const_dict = {
    # --- commonly changed ---
    "duration": 80,
    "dustsp": dustsp,
    "dustsp_2nd": "amnt",
    "dustrate": 1.,   # OVERWRITTEN BY GCAM-DERIVED DUST FLUXES

    # --- secondary mineral rule
    "secondary_min_rule": "sld_track",
    "poro_updated": "None",

    # --- mixing
    "imix": 3,
    "cec_adsorption_on": True,

    # --- surface area / psd rules
    "include_psd_full": True,
    "include_psd_bulk": False,
    "psdrain_log10_sd": 0.05,
    "psdrain_wt": 1.0,
    "use_psdrain_datfile": False,
    "poro_iter_field": False,
    "poro_evol": False,
    "sa_rule1": False,
    "sa_rule2": True,
    "include_roughness_sa": False,

    # --- other
    "singlerun_seasonality": True,
    "skip_lab_run": True,
    "spindir": f"s3://carbonplan-carbon-removal/SCEPTER/scepter_output/spinups/gcam/{spinversion}/",

    # --- CEC update from spinup
    "cec_update_from_spinup": True,

    # --- compute
    "aws_save": "move",
    "scepter_exec_name": "scepter_richards",
}
# ====================================================================


# %%
# ====================================================================
# [ READ SITE DATA + VALIDATE SELECTED REGIONS ]
# ====================================================================
celltypes = "one-cell-per-region"
deg_resolution = 1
dfsites = pd.read_pickle(
    f"s3://carbonplan-carbon-removal/ew-workflows-data/gcam/apprate_calculated/"
    f"{celltypes}/{RUNNAME}/1.0-cropFraction/"
    f"annual_rockflx_densest-ag-region_{deg_resolution}deg.pkl"
)

all_region_names = dfsites["region_nospaces"].unique()
missing = [r for r in SELECTED_REGIONS if r not in all_region_names]
if missing:
    print(
        f"\nWARNING: the following regions in SELECTED_REGIONS were not found in dfsites:\n"
        f"  {missing}\n"
        f"Available regions:\n  {sorted(all_region_names)}\n"
    )

site_rows = (
    dfsites[dfsites["region_nospaces"].isin(SELECTED_REGIONS)][["region_nospaces", "reg_id"]]
    .drop_duplicates(subset="region_nospaces")
)
sites = list(site_rows["region_nospaces"])
reg_ids = list(site_rows["reg_id"])
# ====================================================================


# %%
# ====================================================================
# [ BUILD BASE DATAFRAME ]
# ====================================================================
by_site = {
    "region_nospaces": sites,
    "climatefiles": [f"{site}_{deg_resolution}_{runtype}" for site in sites],
    "spinrun": [
        f"spinup_gcam_{site}_id_{reg_id}_spintuneup4"
        for site, reg_id in zip(sites, reg_ids)
    ],
}

df = bhf.build_df_one_at_a_time(
    pref,
    const_dict,
    sites,
    by_site,
    all_combinations={},
    individual_cases={},
    add_ctrl=True,
    add_index=True,
)
# ====================================================================


# %%
# ====================================================================
# [ READ SPINUP CEC ]
# ====================================================================
spin_batch_path = (
    f"s3://carbonplan-carbon-removal/ew-workflows-data/scepter/batch/"
    f"spinup_gcam_{spinversion}.csv"
)
df_spin_cec = (
    pd.read_csv(spin_batch_path)[["region_nospaces", "cec"]]
    .rename(columns={"cec": "spinup_cec"})
)
# ====================================================================


# %%
# ====================================================================
# [ EXPAND OVER CEC FACTORS × DUSTRADS (full factorial) ]
# ====================================================================
df_ctrl  = df[df["dustrate"] == 0].copy()
df_treat = df[df["dustrate"] != 0].copy()

df_treat = df_treat.merge(df_spin_cec, on="region_nospaces", how="left")
df_ctrl  = df_ctrl.merge(df_spin_cec, on="region_nospaces", how="left")

# --- treatments: all combinations of cec_factor × dustrad
treat_frames = []
for cec_f in cec_factors:
    for dustrad in dustrads:
        frame = df_treat.copy()
        frame["cec"] = frame["spinup_cec"] * cec_f
        frame["cec_factor"] = cec_f
        frame["dustrad"] = dustrad
        treat_frames.append(frame)
df_treat = pd.concat(treat_frames, ignore_index=True)
df_treat = df_treat.drop(columns="spinup_cec")

# --- controls: one per site × cec_factor only (dustrad=dustrads[0] is a placeholder;
#     it is not used in runs with dustrate=0, but must be a real number for psd setup)
ctrl_frames = []
for cec_f in cec_factors:
    frame = df_ctrl.copy()
    frame["cec"] = frame["spinup_cec"] * cec_f
    frame["cec_factor"] = cec_f
    frame["dustrad"] = dustrads[0]
    ctrl_frames.append(frame)
df_ctrl = pd.concat(ctrl_frames, ignore_index=True)
df_ctrl = df_ctrl.drop(columns="spinup_cec")
# ====================================================================


# %%
# --- reunite and compute psdrain_meanRad
df = pd.concat([df_treat, df_ctrl], ignore_index=True)

if const_dict["include_psd_full"] or const_dict["include_psd_bulk"]:
    df["psdrain_meanRad"] = df["dustrad"] / 1e6   # um → m


# %%
# ====================================================================
# [ CLIMATEDIR + AWS_BUCKET ]
# ====================================================================
if not repeat_clim_for_apprates:
    df["climatedir"] = (
        f"s3://carbonplan-carbon-removal/ew-workflows-data/scepter/clim/"
        f"one-cell-per-region_1deg/{EXP_SET}/1.0-cropFraction/"
    )
    df["aws_bucket"] = (
        f"s3://carbonplan-carbon-removal/SCEPTER/scepter_output/{EXP_GROUP}/{EXP_SET}/"
    )
else:
    df_treat = df[df["dustrate"] != 0].copy()
    df_ctrl  = df[df["dustrate"] == 0].copy()

    treat_app_frames = []
    for apprate in apprates:
        frame = df_treat.copy()
        frame["climatedir"] = (
            f"s3://carbonplan-carbon-removal/ew-workflows-data/scepter/clim/"
            f"one-cell-per-region_1deg/{CLIM_SELECTION_LABEL}/{EXP_SET}/{apprate}{apprate_suffix}/"
        )
        frame["aws_bucket"] = (
            f"s3://carbonplan-carbon-removal/SCEPTER/scepter_output/"
            f"{EXP_GROUP}/{EXP_SET}/{apprate}{apprate_suffix}"
        )
        frame["apprate_base"] = apprate
        treat_app_frames.append(frame)
    df_treat = pd.concat(treat_app_frames, ignore_index=True)

    # controls: use first apprate's climatedir; one aws_bucket per cec_factor
    df_ctrl["climatedir"] = (
        f"s3://carbonplan-carbon-removal/ew-workflows-data/scepter/clim/"
        f"one-cell-per-region_1deg/{CLIM_SELECTION_LABEL}/{EXP_SET}/{apprates[0]}{apprate_suffix}/"
    )
    df_ctrl["aws_bucket"] = df_ctrl["cec_factor"].apply(
        lambda f: (
            f"s3://carbonplan-carbon-removal/SCEPTER/scepter_output/"
            f"{EXP_GROUP}/{EXP_SET}/control_cec{str(f).replace('.', 'p')}"
        )
    )
    df_ctrl["cec_update_from_spinup"] = True

    df = pd.concat([df_treat, df_ctrl], ignore_index=True)
# ====================================================================


# %%
# ====================================================================
# [ UPDATE newrun_id ]
# ====================================================================
def _fmt(val):
    """Format a float for embedding in a run ID (replaces '.' with 'p')."""
    return str(val).replace(".", "p")


treat_mask = df["dustrate"] != 0
ctrl_mask  = df["dustrate"] == 0

if repeat_clim_for_apprates:
    df.loc[treat_mask, "newrun_id"] = (
        df.loc[treat_mask, "newrun_id"]
        + "_dustrad" + df.loc[treat_mask, "dustrad"].apply(_fmt)
        + "_cec"     + df.loc[treat_mask, "cec_factor"].apply(_fmt)
        + "_ar"      + df.loc[treat_mask, "apprate_base"].apply(_fmt)
    )
else:
    df.loc[treat_mask, "newrun_id"] = (
        df.loc[treat_mask, "newrun_id"]
        + "_dustrad" + df.loc[treat_mask, "dustrad"].apply(_fmt)
        + "_cec"     + df.loc[treat_mask, "cec_factor"].apply(_fmt)
    )

df.loc[ctrl_mask, "newrun_id"] = (
    df.loc[ctrl_mask, "newrun_id"]
    + "_cec" + df.loc[ctrl_mask, "cec_factor"].apply(_fmt)
)
# ====================================================================


# %%
# --- SAVE
bhf.save_df(df, savepath_batch, fn, multi_run_split, max_iters_per_set)

# %%
pref
# %%
