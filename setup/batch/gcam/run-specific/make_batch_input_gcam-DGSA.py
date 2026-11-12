# ---------------------------------------------------
#
#      !! BUILDS CSV FOR "DGSA" RUNS !!
#
# Generate batch inputs for a Distance-based Generalized
# Sensitivity Analysis (DGSA). Uses Latin Hypercube
# Sampling (LHS) to draw from parameter distributions
# for apprate, cec_factor, roughness_factor, dustrad_init,
# and qrun_factor across 7 user-selected sites.
#
# For every (site, LHS sample) pair there is:
#   - one treatment run (rock applied)
#   - one control run  (dustrate=0, same cec_factor)
#
# T Kukla ; july 2026
#
# ---------------------------------------------------
# %%
import numpy as np
import pandas as pd
from scipy.stats import norm, qmc, truncnorm

# --- read in batch helper functions
from ew_workflows import batch_helperFxns as bhf
# ---


# %%
#
#
# [ UNIVERSAL ]
RUNNAME = "exe_central1p5"
DGSA_VERSION = "v1.2"               # bump to avoid overwriting prior DGSA climatedir sets
EXP_GROUP = f"gcam_dgsa_{DGSA_VERSION}"  
EXP_SET = f"{RUNNAME}"
CLIM_SELECTION_LABEL = "repcells1"  # see labels in `create_clim_inputs_one-cell-per-region.ipynb`
dustsp = "gbas"
spinversion = "v4_repcells1"
runtype = "yearly_ltm"
extra_tag = spinversion
pref = f"{EXP_GROUP}_{EXP_SET}_{dustsp}_{runtype}_{extra_tag}"
fn = pref + ".csv"

# climatedir config
apprate_suffix = "-tHaYr-appRate-FIXED"
ref_apprate = 5.0   # [t/ha/yr] existing source dir to copy+scale from
                    # must be one of: [0.05, 0.1, 0.5, 1., 2.5, 5., 7.5, 10., 15., 20., 30., 40., 50.]
dryrun_climatedirs = False   # True = skip S3 writes (test path logic only)

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
    # [ 75% of rock application and, usually, cdr_adv ]
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
# [ SAMPLING CONFIG ]
# ====================================================================
# (scroll to bottom for visualization of uniform vs log-uniform)
N_SAMPLES = 250   # !! number of LHS draws (= simulations per site) !!
RANDOM_SEED = 42  # for reproducibility

# When True:  sample dustrad as a single LHS dimension (1 degree of freedom).
# When False: sample roughness_factor and dustrad_init independently (2 degrees
#             of freedom) and derive dustrad = dustrad_init / roughness_factor.
SAMPLE_DUSTRAD_COMBINED = True

# Distribution config per parameter.
# Supported dist types: "uniform", "log-uniform", "normal", "truncated_normal"
#   uniform:           params = {"low": ..., "high": ...}
#   log-uniform:       params = {"low": ..., "high": ...}   (bounds in original units)
#   normal:            params = {"loc": ..., "scale": ...}
#   truncated_normal:  params = {"loc": ..., "scale": ..., "low": ..., "high": ...}
PARAM_DISTS = {
    "apprate": {                         # [t/ha/yr] application rate
        "dist": "uniform",
        "params": {"low": 0.2, "high": 20.},
    },
    "cec_factor": {                      # [] multiplier on spinup CEC
        "dist": "uniform",
        "params": {"low": 0.1, "high": 4.},
    },
    "roughness_factor": {                # [] dustrad = dustrad_init / roughness_factor
        "dist": "uniform",
        "params": {"low": 2., "high": 30.},  # Strefler et al. has RF ~30 for dustrad ~10 um.
    },
    "dustrad_init": {                    # [micron] initial dust particle radius
        "dist": "uniform",
        "params": {"low": 10., "high": 50.},  # GCAM default is 10; set as lower bound
    },
    "qrun_factor": {                     # [] multiplier on q_temp.in runoff values
        "dist": "log-uniform",
        "params": {"low": 0.25, "high": 4.},
    },
}

# When SAMPLE_DUSTRAD_COMBINED is True, swap roughness_factor + dustrad_init for a
# single dustrad dimension. Bounds are derived from the limits already defined above,
# so dustrad spans exactly the range achievable by dustrad_init / roughness_factor.
# roughness_factor and dustrad_init will be set to NaN in the dataframe (not sampled).
if SAMPLE_DUSTRAD_COMBINED:
    _rf = PARAM_DISTS.pop("roughness_factor")["params"]
    _di = PARAM_DISTS.pop("dustrad_init")["params"]
    PARAM_DISTS["dustrad"] = {
        "dist": "log-uniform",
        "params": {
            "low":  _di["low"]  / _rf["high"],  # dustrad_init.low  / roughness_factor.high
            "high": _di["high"] / _rf["low"],   # dustrad_init.high / roughness_factor.low
        },
    }
# ====================================================================


# %%
# ====================================================================
# [ CONSTANTS ]
# ====================================================================
const_dict = {
    # --- commonly changed ---
    "duration": 80,          # [yr]
    "dustsp": dustsp,
    "dustsp_2nd": "amnt",
    "dustrate": 1.,          # OVERWRITTEN BY GCAM-DERIVED DUST FLUXES

    # --- secondary mineral rule
    "secondary_min_rule": "sld_track",
    "poro_updated": "None",

    # --- mixing
    "imix": 3,               # 1=fickian; 2=homogeneous; 3=tilling
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
# [ LHS SAMPLING ]
# ====================================================================
def _sample_from_dist(u, cfg):
    """Transform uniform [0,1] samples u to the distribution defined by cfg."""
    dist = cfg["dist"]
    p = cfg["params"]
    if dist == "uniform":
        return p["low"] + (p["high"] - p["low"]) * u
    elif dist == "log-uniform":
        log_low = np.log10(p["low"])
        log_high = np.log10(p["high"])
        return 10 ** (log_low + (log_high - log_low) * u)
    elif dist == "normal":
        return norm.ppf(u, loc=p["loc"], scale=p["scale"])
    elif dist == "truncated_normal":
        a = (p["low"] - p["loc"]) / p["scale"]
        b = (p["high"] - p["loc"]) / p["scale"]
        return truncnorm.ppf(u, a, b, loc=p["loc"], scale=p["scale"])
    else:
        raise ValueError(f"Unknown distribution type: '{dist}'")


param_names = list(PARAM_DISTS.keys())
sampler = qmc.LatinHypercube(d=len(param_names), seed=RANDOM_SEED)
unit_samples = sampler.random(n=N_SAMPLES)   # shape (N_SAMPLES, len(param_names))

lhs_df = pd.DataFrame(
    {name: _sample_from_dist(unit_samples[:, i], PARAM_DISTS[name])
     for i, name in enumerate(param_names)}
)
lhs_df.index.name = "lhs_sample_id"
lhs_df = lhs_df.reset_index()   # lhs_sample_id column: 0 … N_SAMPLES-1
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
# [ BUILD TREATMENT DATAFRAME ]
# One row per (site, LHS sample) → N_SAMPLES × len(sites) rows
# ====================================================================
treat_frames = []
for site, reg_id in zip(sites, reg_ids):
    frame = lhs_df.copy()
    frame["region_nospaces"] = site
    frame["reg_id"] = reg_id
    treat_frames.append(frame)
df_treat = pd.concat(treat_frames, ignore_index=True)
df_treat["site"] = df_treat["region_nospaces"]

# --- CEC and dustrad from LHS samples
df_treat = df_treat.merge(df_spin_cec, on="region_nospaces", how="left")
df_treat["cec"] = df_treat["spinup_cec"] * df_treat["cec_factor"]
if SAMPLE_DUSTRAD_COMBINED:
    # dustrad sampled directly; roughness_factor and dustrad_init were not LHS dims
    df_treat["roughness_factor"] = np.nan
    df_treat["dustrad_init"] = np.nan
else:
    df_treat["dustrad"] = df_treat["dustrad_init"] / df_treat["roughness_factor"]
df_treat = df_treat.drop(columns="spinup_cec")

# --- constants
for k, v in const_dict.items():
    df_treat[k] = v

# --- by-site columns
df_treat["climatefiles"] = df_treat["region_nospaces"] + f"_{deg_resolution}_{runtype}"
df_treat["spinrun"] = (
    "spinup_gcam_"
    + df_treat["region_nospaces"]
    + "_id_"
    + df_treat["reg_id"].astype(str)
    + "_spintuneup4"
)

# --- apprate stored as apprate_base (actual dustrate overwritten later by GCAM-derived fluxes)
df_treat = df_treat.rename(columns={"apprate": "apprate_base"})

# --- assign integer paramSet ID (1-indexed) from lhs_sample_id so all sites sharing
#     the same LHS draw group under the same paramSet_* directory.
df_treat["paramSet_id"] = df_treat["lhs_sample_id"] + 1

# --- climatedir: copy + scale from reference apprate dir; one dir per row (paramSet_id).
df_treat["climatedir"] = bhf.build_dgsa_climatedirs(
    df_treat=df_treat,
    base_clim_s3="s3://carbonplan-carbon-removal/ew-workflows-data/scepter/clim/one-cell-per-region_1deg",
    clim_selection_label=CLIM_SELECTION_LABEL,
    exp_set=EXP_SET,
    dgsa_version=DGSA_VERSION,
    ref_apprate=ref_apprate,
    apprate_suffix=apprate_suffix,
    runtype=runtype,
    deg_resolution=deg_resolution,
    dryrun=dryrun_climatedirs,
    save_validation_plots = True,
)
df_treat["aws_bucket"] = (
    f"s3://carbonplan-carbon-removal/SCEPTER/scepter_output/{EXP_GROUP}/{EXP_SET}/"
)

# --- psd tracking (must come after dustrad is set)
if const_dict["include_psd_full"] or const_dict["include_psd_bulk"]:
    df_treat["psdrain_meanRad"] = df_treat["dustrad"] / 1e6   # um → m

# --- unique run ID
df_treat["newrun_id"] = (
    pref + "_"
    + df_treat["region_nospaces"]
    + "_lhs" + df_treat["lhs_sample_id"].astype(str)
)
# ====================================================================


# %%
# ====================================================================
# [ BUILD CONTROL DATAFRAME ]
# One control per (site, LHS sample) — matches each treatment in CEC
# but with dustrate=0 (no rock application).
# ====================================================================
df_ctrl = df_treat.copy()
df_ctrl["dustrate"] = 0.
if SAMPLE_DUSTRAD_COMBINED:
    # dustrad_init is NaN; keep dustrad as inherited from df_treat (placeholder, not used)
    pass
else:
    df_ctrl["dustrad"] = df_ctrl["dustrad_init"]   # placeholder; not used when dustrate=0
if const_dict["include_psd_full"] or const_dict["include_psd_bulk"]:
    df_ctrl["psdrain_meanRad"] = df_ctrl["dustrad"] / 1e6
df_ctrl["newrun_id"] = (
    pref + "_ctrl_"
    + df_ctrl["region_nospaces"]
    + "_lhs" + df_ctrl["lhs_sample_id"].astype(str)
)
# ====================================================================


# %%
# --- COMBINE AND SAVE
df = pd.concat([df_treat, df_ctrl], ignore_index=True)
# %% 
bhf.save_df(df, savepath_batch, fn, multi_run_split, max_iters_per_set)

# %%
pref
# %%
#
#
#
# %%
#
#
#
#
#
#
#
# ————————————————————————————————————————————————————————————————————————
# --- SCRATCH quickplot uniform vs log-uniform dists
import matplotlib.pyplot as plt
A, B = 0.5, 2
n_samples = 1000

uniform_samples = np.random.uniform(A, B, n_samples)
log_uniform_samples = np.exp(np.random.uniform(np.log(A), np.log(B), n_samples))

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].hist(uniform_samples, bins=50, color='steelblue', edgecolor='white', alpha=0.8)
axes[0].set_title(f'Uniform({A}, {B})')
axes[0].set_xlabel('Value')
axes[0].set_ylabel('Count')

axes[1].hist(log_uniform_samples, bins=50, color='coral', edgecolor='white', alpha=0.8)
axes[1].set_title(f'Log-Uniform({A}, {B})')
axes[1].set_xlabel('Value')
# axes[1].set_xscale('log')

plt.suptitle('Uniform vs Log-Uniform Sampling', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.show()
# %%
