"""
find_representative_cropland_cell_all_regions.py
─────────────────────────────────────────────────
For each GCAM region, identify the ERA5 grid cell whose spinup variables
are closest to the cropland-area-weighted regional mean.

Saves per-region results to:
  {S3_OUT_ROOT}/{region_nospaces}/
      spinup_vars.csv    — one-row CSV; pd.concat() across dirs for a full table
      cell_map.png       — 9-panel map of each spinup variable
      distance_map.png   — normalized distance from cropland-weighted centroid
      gridded_data.zarr  — all ERA5 cells in the region with every derived variable
                           (ERA5 climate, SoilGrids, ISRIC, MODIS NPP, GLDAS porosity,
                            spinup vars, cropland fraction, rep_distance)
"""
# %%
import io
import traceback

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import s3fs
import icechunk
import rioxarray  # noqa: F401
from pyproj import Transformer
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import obstore
try:
    from obspec_utils.registry import ObjectStoreRegistry
except ImportError:
    from virtualizarr.registry import ObjectStoreRegistry
from virtual_tiff import VirtualTIFF
import rasterix
from affine import Affine

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# Water infiltration proxy for qrun:
#   "pme" — P − E  (can be negative in arid/cold regions)
#   "tp"  — total precipitation  (always positive; use when P−E is unreasonably low)
RUNOFF_PROXY = "tp"

# If True, skip any region whose output dir already contains all three expected
# files (spinup_vars.csv, cell_map.png, distance_map.png).  Set to False to
# recompute and overwrite everything.
SKIP_EXISTING = True

YEAR_RANGE   = [2000, 2024]
MODIS_YEAR   = 2000
CO2_ATM_PPMV = 420
A_CONST      = 1.03e6
FRACT_NPP    = 0.75
REP_VARS     = ["tsom", "tph", "tec", "tsoilco2", "cec", "qrun", "mat", "soilmoisture"]

# Optional hard constraints applied before the representative-cell search.
# Keys must be in REP_VARS; omit a variable (or leave its dict empty) for no constraint.
# Units match derive_spinup_vars() output:
#   qrun         → m yr⁻¹        (100 mm yr⁻¹ = 0.1)
#   tph          → pH units
#   mat          → °C
#   soilmoisture → m³ m⁻³
#   tsom         → wt%
#   tec          → % CEC
#   cec          → cmol kg⁻¹
#   tsoilco2     → log₁₀(atm)
VAR_CONSTRAINTS = {
    "qrun": {"min": 0.1},            # > 100 mm yr⁻¹
    "tph":  {"min": 5.4, "max": 7.0},
}

S3_OUT_ROOT          = f"s3://carbonplan-carbon-removal/ew-workflows-data/gcam/spinup-conditions/cropwtd-mean/{RUNOFF_PROXY}-qrun"
ERA5_BUCKET          = "carbonplan-carbon-removal"
ERA5_ICECHUNK_DIR    = "era5/preprocessed_icechunk"
VALIDATION_PREFIX    = "ew-workflows-data/valdation-data"
ISRIC_BASE           = f"s3://carbonplan-carbon-removal/{VALIDATION_PREFIX}/ISRIC"
S3_CROPLAND_ZARR     = (
    "s3://carbonplan-carbon-removal/ew-workflows-data/gcam/"
    "cropland+gcamRegions/ds_cropland+gcam_1deg.zarr/"
)
S3_GCAM_SHP = (
    "s3://carbonplan-carbon-removal/ew-workflows-data/gcam/"
    "gcam_boundaries_moirai_3p1_0p5arcmin_wgs84/main_outputs/"
    "region_boundaries_moirai_landcells_3p1_0p5arcmin.shp"
)

SG_VARS   = ["phh2o", "soc", "bdod", "cec"]
SG_DEPTHS = ["0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm"]
SG_BASE   = "https://files.isric.org/"
SG_RES    = 5000

ISRIC_VARS   = ["BSAT", "CECS", "CECc", "ECEC", "PHAQ", "ORGC"]
ISRIC_DEPTHS = ["D1", "D2", "D3", "D4", "D5", "D6", "D7"]
ISRIC_LABELS = ["0-20cm", "20-40cm", "40-60cm", "60-80cm", "80-100cm", "100-150cm", "150-200cm"]

# ERA5 variable names (long CDS names used in the carbonplan icechunk store)
ERA5_VARS = {
    "tp":    "total_precipitation",
    "e":     "evaporation",
    "swvl1": "volumetric_soil_water_layer_1",
    "swvl2": "volumetric_soil_water_layer_2",
    "swvl3": "volumetric_soil_water_layer_3",
    "stl1":  "soil_temperature_level_1",
}
# ──────────────────────────────────────────────────────────────────────────────


# ── GLOBAL DATASET LOADERS ────────────────────────────────────────────────────

def open_era5(era5_vars):
    """Open ERA5 icechunk store; fall back short→long variable names."""
    storage = icechunk.s3_storage(bucket=ERA5_BUCKET, prefix=ERA5_ICECHUNK_DIR, from_env=True)
    repo    = icechunk.Repository.open(storage)
    session = repo.readonly_session("main")
    ds      = xr.open_zarr(session.store, consolidated=False)
    available = set(ds.data_vars)
    resolved  = {}
    for short, long in era5_vars.items():
        if long in available:
            resolved[short] = long
        elif short in available:
            resolved[short] = short
        else:
            print(f"  ERA5 WARNING: '{long}' / '{short}' not found")
    return ds, resolved


def open_soilgrids():
    """Return lazy SoilGrids dataset (ESRI:54052) and Affine transform."""
    sg_store    = obstore.store.from_url(SG_BASE)
    sg_registry = ObjectStoreRegistry({SG_BASE: sg_store})
    parser      = VirtualTIFF(ifd=0)

    def _manifest(var, depth):
        url = (
            f"{SG_BASE}soilgrids/latest/data_aggregated/{SG_RES}m/"
            f"{var}/{var}_{depth}_mean_{SG_RES}.tif"
        )
        return xr.open_zarr(
            parser(url=url, registry=sg_registry), zarr_format=3, consolidated=False
        )["0"]

    da_ref    = _manifest("phh2o", "0-5cm")
    mps       = da_ref.attrs["model_pixel_scale"]
    mtp       = da_ref.attrs["model_tiepoint"]
    transform = Affine(mps[0], 0.0, mtp[3], 0.0, -mps[1], mtp[4])

    depth_slices = []
    for depth in SG_DEPTHS:
        das = []
        for var in SG_VARS:
            da = _manifest(var, depth)
            da = da.proj.assign_crs(spatial_ref="ESRI:54052")
            idx = rasterix.RasterIndex.from_transform(
                transform, width=da.sizes["x"], height=da.sizes["y"]
            )
            da = da.assign_coords(xr.Coordinates.from_xindex(idx))
            das.append(da.rename(var))
        depth_slices.append(xr.merge(das))

    sg_ds = xr.concat(depth_slices, dim=pd.Index(SG_DEPTHS, name="depth"))
    return sg_ds, transform


def open_isric():
    """Open ISRIC raster and depth property tables."""
    isric_raster     = xr.open_dataset(
        f"{ISRIC_BASE}/wise_30sec_v1.tif", engine="rasterio", chunks="auto",
    )["band_data"].squeeze("band")
    tsv              = pd.read_csv(
        f"{ISRIC_BASE}/wise_30sec_v1.tsv", sep="\t"
    ).rename(columns={"pixel_vaue": "pixel_value"})
    pixel_to_newsuid = tsv.set_index("pixel_value")["description"]
    depth_tables     = {}
    for d, label in zip(ISRIC_DEPTHS, ISRIC_LABELS):
        df_d = pd.read_csv(f"{ISRIC_BASE}/HW30s_w{d}.txt", quotechar='"')
        depth_tables[label] = df_d.set_index("NEWSUID")[ISRIC_VARS]
    return isric_raster, pixel_to_newsuid, depth_tables


# ── PER-REGION PROCESSING ─────────────────────────────────────────────────────

def era5_ltm_for_region(ds_era5, resolved_vars, bounds, year_range):
    """Compute ERA5 long-term means within the region bounding box."""
    lat_dim  = "latitude" if "latitude" in ds_era5.dims else "lat"
    lon_dim  = "longitude" if "longitude" in ds_era5.dims else "lon"
    minx, miny, maxx, maxy = bounds
    buf = 1.0
    ds_sub = ds_era5.sel({
        "time":   slice(f"{year_range[0]}-01-01", f"{year_range[1]}-12-31"),
        lat_dim:  slice(maxy + buf, miny - buf),
        lon_dim:  slice(minx - buf, maxx + buf),
    })
    ltm = {short: ds_sub[varname].mean("time").compute()
           for short, varname in resolved_vars.items()
           if varname in ds_sub.data_vars}
    lats = ltm["tp"][lat_dim].values
    lons = ltm["tp"][lon_dim].values
    return ltm, lats, lons


def build_cell_gdf(ltm, lats, lons, region_geom):
    """Flatten ERA5 LTM arrays to a GeoDataFrame filtered to the region polygon."""
    lat_grid, lon_grid = np.meshgrid(lats, lons, indexing="ij")
    df = pd.DataFrame({
        "lat":      lat_grid.flatten(),
        "lon":      lon_grid.flatten(),
        "pme_m_hr": (ltm["tp"].values + ltm["e"].values).flatten(),
        "tp_m_hr":  ltm["tp"].values.flatten(),
        "swvl1":    ltm["swvl1"].values.flatten() if "swvl1" in ltm else np.nan,
        "swvl2":    ltm["swvl2"].values.flatten() if "swvl2" in ltm else np.nan,
        "swvl3":    ltm["swvl3"].values.flatten() if "swvl3" in ltm else np.nan,
        "stl1_K":   ltm["stl1"].values.flatten() if "stl1" in ltm else np.nan,
    })
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["lon"], df["lat"]), crs="EPSG:4326")
    return gdf[gdf.within(region_geom)].reset_index(drop=True)


def sample_cropland(da_crop, gdf):
    vals = da_crop.interp(
        y=xr.DataArray(gdf["lat"].values, dims="points"),
        x=xr.DataArray(gdf["lon"].values, dims="points"),
        method="nearest",
    ).compute().values
    return np.where(np.isfinite(vals), vals, 0.0).clip(min=0.0)


def sample_modis(modis_da, gdf):
    raw = (
        modis_da
        .sel(
            x=xr.DataArray(gdf["lon"].values, dims="points"),
            y=xr.DataArray(gdf["lat"].values, dims="points"),
            method="nearest",
        )
        .compute().values
    )
    return raw * 0.0001 * 1e3  # → gC m⁻² yr⁻¹


def sample_gldas(gldas_da, gdf):
    return (
        gldas_da
        .sel(
            lon=xr.DataArray(gdf["lon"].values, dims="points"),
            lat=xr.DataArray(gdf["lat"].values, dims="points"),
            method="nearest",
        )
        .compute().values
    )


def sample_soilgrids(sg_ds, gdf):
    """Sample SoilGrids at ERA5 cell centres (EPSG:4326 → ESRI:54052 conversion)."""
    transformer = Transformer.from_crs("EPSG:4326", "ESRI:54052", always_xy=True)
    xs, ys      = transformer.transform(gdf["lon"].values, gdf["lat"].values)
    sg_sampled  = sg_ds.sel(
        x=xr.DataArray(xs, dims="points"),
        y=xr.DataArray(ys, dims="points"),
        method="nearest",
    )
    sg_df = (
        sg_sampled.compute()
        .drop_vars(["x", "y", "spatial_ref"], errors="ignore")
        .to_dataframe()
        .unstack("depth")
    )
    sg_df.columns = [f"SoilGrids_{var}_{depth}" for var, depth in sg_df.columns]
    return sg_df.reset_index(drop=True)


def _nearest_newsuid(isric_raster, pixel_to_newsuid, valid_pixel_set, lon, lat):
    for half_w in [0.5, 1.0, 2.0, 5.0, 10.0, 15.0]:
        xi = np.where(
            (isric_raster.x.values >= lon - half_w) & (isric_raster.x.values <= lon + half_w)
        )[0]
        yi = np.where(
            (isric_raster.y.values >= lat - half_w) & (isric_raster.y.values <= lat + half_w)
        )[0]
        if xi.size == 0 or yi.size == 0:
            continue
        window = isric_raster.isel(x=xi, y=yi).compute()
        pvals  = window.values.astype(int)
        valid_mask = np.isin(pvals, list(valid_pixel_set))
        if not valid_mask.any():
            continue
        ys_w, xs_w = np.meshgrid(window.y.values, window.x.values, indexing="ij")
        dists = (xs_w - lon) ** 2 + (ys_w - lat) ** 2
        dists[~valid_mask] = np.inf
        iy, ix = np.unravel_index(dists.argmin(), dists.shape)
        return pixel_to_newsuid[int(pvals[iy, ix])]
    return np.nan


def sample_isric(isric_raster, pixel_to_newsuid, depth_tables, gdf):
    valid_pixel_set = set(pixel_to_newsuid.index)
    pixel_vals = (
        isric_raster
        .sel(
            x=xr.DataArray(gdf["lon"].values, dims="points"),
            y=xr.DataArray(gdf["lat"].values, dims="points"),
            method="nearest",
        )
        .compute().values.astype(int)
    )
    newsuids = pixel_to_newsuid.reindex(pixel_vals).values
    nan_idxs = np.where(pd.isna(newsuids))[0]
    if len(nan_idxs):
        for i in nan_idxs:
            newsuids[i] = _nearest_newsuid(
                isric_raster, pixel_to_newsuid, valid_pixel_set,
                gdf["lon"].iloc[i], gdf["lat"].iloc[i],
            )
    isric_df = pd.DataFrame(index=range(len(gdf)))
    for label, props_df in depth_tables.items():
        joined = props_df.reindex(newsuids).reset_index(drop=True)
        for var in ISRIC_VARS:
            isric_df[f"ISRIC_{var}_{label}"] = joined[var].values
    return isric_df


def derive_spinup_vars(gdf, runoff_proxy):
    _col = "pme_m_hr" if runoff_proxy == "pme" else "tp_m_hr"
    gdf = gdf.copy()
    gdf["qrun"]         = gdf[_col] * 8760
    gdf["soilmoisture"] = (gdf["swvl1"] * 7 + gdf["swvl2"] * 21 + gdf["swvl3"] * 72) / 100
    gdf["mat"]          = gdf["stl1_K"] - 273.15
    gdf["soil_pco2_GB94"] = (
        CO2_ATM_PPMV
        + (A_CONST * FRACT_NPP * gdf["MODIS_NPP_gC_m2_yr"]) / (gdf["stl1_K"] ** 2)
    )
    gdf["tsoilco2"] = np.log10(gdf["soil_pco2_GB94"] / 1e6)
    gdf["tsom"] = (
        gdf["SoilGrids_soc_0-5cm"]   * (1 / 6)
        + gdf["SoilGrids_soc_5-15cm"]  * (2 / 6)
        + gdf["SoilGrids_soc_15-30cm"] * (3 / 6)
    ) / 100
    gdf["tph"]  = gdf["SoilGrids_phh2o_15-30cm"] / 10
    gdf["tec"]  = 101 - gdf["ISRIC_BSAT_0-20cm"]
    gdf["cec"]  = ((gdf["SoilGrids_cec_5-15cm"] + gdf["SoilGrids_cec_15-30cm"]) / 2) / 10
    return gdf


def find_representative_cell(gdf, rep_vars, var_constraints=None):
    # Build mask of cells that satisfy all hard constraints
    constraint_mask = np.ones(len(gdf), dtype=bool)
    if var_constraints:
        for v, bounds in var_constraints.items():
            if v not in gdf.columns:
                continue
            vals = gdf[v].values.astype(float)
            if "min" in bounds:
                constraint_mask &= np.isfinite(vals) & (vals >= bounds["min"])
            if "max" in bounds:
                constraint_mask &= np.isfinite(vals) & (vals <= bounds["max"])
    gdf_cand = gdf[constraint_mask].reset_index(drop=True)
    if len(gdf_cand) == 0:
        raise ValueError(f"No cells satisfy var_constraints={var_constraints}")

    # Weighted mean/std computed over constrained candidates only
    weights  = gdf_cand["cropland_frac"].values.copy()
    wt_mean  = {}
    wt_std   = {}
    for v in rep_vars:
        vals = gdf_cand[v].values.astype(float)
        mask = np.isfinite(vals) & np.isfinite(weights)
        w, x = weights[mask], vals[mask]
        wt_mean[v] = np.average(x, weights=w)
        wt_std[v]  = np.sqrt(np.average((x - wt_mean[v]) ** 2, weights=w))

    # Compute distance for all cells so maps show the full picture
    dist2 = np.zeros(len(gdf))
    for v in rep_vars:
        std = wt_std[v]
        if std < 1e-12:
            continue
        dist2 += ((gdf[v].values.astype(float) - wt_mean[v]) / std) ** 2
    gdf = gdf.copy()
    gdf["rep_distance"] = np.sqrt(dist2)

    # Best cell chosen only from constrained candidates
    cand_indices = gdf.index[constraint_mask].tolist()
    best_idx = int(gdf.loc[cand_indices, "rep_distance"].idxmin())
    return gdf, gdf.loc[best_idx], wt_mean, wt_std


# ── FIGURES ───────────────────────────────────────────────────────────────────

def _map_extent(bounds, pad=2.0):
    minx, miny, maxx, maxy = bounds
    return [minx - pad, maxx + pad, miny - pad, maxy + pad]


def make_cell_map(gdf, rep, bounds, region_name, region_id, runoff_proxy):
    qrun_label = "P−E [m yr⁻¹]" if runoff_proxy == "pme" else "Precip [m yr⁻¹]"
    plot_vars  = [
        ("cropland_frac", "Cropland fraction",     cmocean.cm.algae,   None, None),
        ("tsom",          "SOM [wt%]",             cmocean.cm.matter,  None, None),
        ("tph",           "Porewater pH",          cmocean.cm.haline,  4.5,  8.5),
        ("tec",           "Exch. acidity [%CEC]",  cmocean.cm.solar,   0,    100),
        ("tsoilco2",      "Soil CO₂ [log₁₀ atm]", cmocean.cm.deep,    None, None),
        ("cec",           "CEC [cmol kg⁻¹]",       cmocean.cm.tempo,   None, None),
        ("qrun",          qrun_label,              cmocean.cm.rain,    None, None),
        ("mat",           "Mean ann. T [°C]",      cmocean.cm.thermal, None, None),
        ("soilmoisture",  "Soil moisture [m³ m⁻³]", cmocean.cm.rain,  0,    0.6),
    ]
    ncols, nrows = 3, 3
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5.5 * ncols, 4 * nrows),
        subplot_kw={"projection": ccrs.PlateCarree()}, constrained_layout=True,
    )
    extent = _map_extent(bounds)
    for ax, (vname, title, cmap, vmin, vmax) in zip(axes.flat, plot_vars):
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.BORDERS,   linewidth=0.3, alpha=0.5)
        sc = ax.scatter(
            gdf["lon"].values, gdf["lat"].values,
            c=gdf[vname].values, cmap=cmap, vmin=vmin, vmax=vmax,
            s=60, transform=ccrs.PlateCarree(), zorder=3,
        )
        ax.scatter(
            [rep["lon"]], [rep["lat"]], s=220, marker="*",
            color="red", edgecolors="white", linewidths=0.8,
            transform=ccrs.PlateCarree(), zorder=5,
        )
        plt.colorbar(sc, ax=ax, orientation="vertical", pad=0.02, shrink=0.8)
        ax.set_title(title, fontsize=10)
    fig.suptitle(
        f"Spinup variables — {region_name} (GCAM id={region_id})\n"
        f"★ = representative cell  (lat={rep['lat']:.3f}, lon={rep['lon']:.3f})",
        fontsize=12,
    )
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def make_distance_map(gdf, rep, bounds, region_name):
    fig, ax = plt.subplots(1, 1, figsize=(7, 5), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent(_map_extent(bounds), crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeature.BORDERS,   linewidth=0.3, alpha=0.5)
    sc = ax.scatter(
        gdf["lon"].values, gdf["lat"].values,
        c=gdf["rep_distance"].values, cmap="viridis_r",
        s=60, transform=ccrs.PlateCarree(), zorder=3,
    )
    ax.scatter(
        [rep["lon"]], [rep["lat"]], s=220, marker="*",
        color="red", edgecolors="white", linewidths=0.8,
        transform=ccrs.PlateCarree(), zorder=5,
    )
    plt.colorbar(sc, ax=ax, label="Normalized distance from cropland-weighted mean")
    ax.set_title(
        f"Distance from cropland-weighted centroid — {region_name}\n"
        f"★ = representative cell  (lat={rep['lat']:.3f}, lon={rep['lon']:.3f})",
        fontsize=10,
    )
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


# ── SAVE TO S3 ────────────────────────────────────────────────────────────────

def save_region(fs, reg_id, region_name, region_nospaces, gdf, rep, wt_mean, wt_std, bounds):
    s3_dir = f"{S3_OUT_ROOT}/{region_nospaces}"
    fs.makedirs(s3_dir, exist_ok=True)

    # CSV — one row with all spinup var values (easy to pd.concat later)
    result = pd.DataFrame([{
        "reg_id":          reg_id,
        "region":          region_name,
        "region_nospaces": region_nospaces,
        "lat":             rep["lat"],
        "lon":             rep["lon"],
        "cropland_frac":   rep["cropland_frac"],
        "rep_distance":    rep["rep_distance"],
        "runoff_proxy":    RUNOFF_PROXY,
        "poro":            float(rep.get("soil_porosity", np.nan)),
        **{v: rep[v] for v in REP_VARS},
        **{f"wt_mean_{v}": wt_mean[v] for v in REP_VARS},
        **{f"wt_std_{v}":  wt_std[v]  for v in REP_VARS},
    }])
    with fs.open(f"{s3_dir}/spinup_vars.csv", "w") as f:
        result.to_csv(f, index=False)

    # Figures
    buf_cell = make_cell_map(gdf, rep, bounds, region_name, reg_id, RUNOFF_PROXY)
    with fs.open(f"{s3_dir}/cell_map.png", "wb") as f:
        f.write(buf_cell.read())

    buf_dist = make_distance_map(gdf, rep, bounds, region_name)
    with fs.open(f"{s3_dir}/distance_map.png", "wb") as f:
        f.write(buf_dist.read())

    # Zarr — full gridded data for all ERA5 cells in the region
    ds_grid = (
        gdf.drop(columns="geometry")
        .reset_index(drop=True)
        .to_xarray()
        .rename({"index": "points"})
    )
    store = fs.get_mapper(f"{s3_dir}/gridded_data.zarr")
    ds_grid.to_zarr(store, mode="w")

    print(f"  → saved to {s3_dir}/")


# ── COMPILE ALL REGIONS ───────────────────────────────────────────────────────

def compile_results(fs, s3_out_root, label="repcells1"):
    """Concatenate all per-region spinup_vars.csv into a single GeoDataFrame parquet.

    The label is embedded in the parquet filename so different selection strategies
    (e.g. repcells1, repcells2) can coexist without overwriting each other.

    If existing CSVs are missing the 'poro' column (written before the save_region
    update), GLDAS porosity is sampled once and back-filled.
    """
    bucket_prefix = s3_out_root.removeprefix("s3://")
    csv_paths = sorted(fs.glob(f"{bucket_prefix}/*/spinup_vars.csv"))
    if not csv_paths:
        print("  No spinup_vars.csv files found — skipping compile.")
        return None

    dfs = []
    for p in csv_paths:
        with fs.open(p) as f:
            dfs.append(pd.read_csv(f))
    df = pd.concat(dfs, ignore_index=True).sort_values("reg_id").reset_index(drop=True)

    # Back-fill porosity for CSVs written before the 'poro' column was added
    if "poro" not in df.columns or df["poro"].isna().any():
        print("  Fetching GLDAS porosity for missing poro values ...")
        gldas_da = xr.open_dataset(
            "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/gldas/SOILS/GLDASp5_porosity_025d.nc4",
            engine="h5netcdf", chunks="auto",
        )["GLDAS_porosity"].squeeze("time")
        poro_vals = gldas_da.sel(
            lon=xr.DataArray(df["lon"].values, dims="points"),
            lat=xr.DataArray(df["lat"].values, dims="points"),
            method="nearest",
        ).compute().values
        if "poro" not in df.columns:
            df["poro"] = poro_vals
        else:
            df["poro"] = df["poro"].where(df["poro"].notna(), poro_vals)

    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["lon"], df["lat"]), crs="EPSG:4326")
    out_path = f"{s3_out_root}/all_regions_spinup_vars_{label}.parquet"
    with fs.open(out_path, "wb") as f:
        gdf.to_parquet(f)
    print(f"  Compiled {len(gdf)} regions → {out_path}")
    return gdf


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    fs = s3fs.S3FileSystem(anon=False)

    print("Loading GCAM shapefile ...")
    gdf_regions = (
        gpd.read_file(S3_GCAM_SHP)[["reg_id", "reg_nm", "geometry"]]
        .to_crs("EPSG:4326")
        .dissolve(by="reg_id")
        .reset_index()
    )
    gdf_regions["region_nospaces"] = gdf_regions["reg_nm"].str.replace(" ", "_", regex=False)

    print("Opening ERA5 icechunk store ...")
    ds_era5, resolved_vars = open_era5(ERA5_VARS)

    print("Opening SoilGrids virtual tiff stack ...")
    sg_ds, _ = open_soilgrids()

    print("Opening ISRIC ...")
    isric_raster, pixel_to_newsuid, depth_tables = open_isric()

    print("Opening MODIS NPP ...")
    modis_url = (
        f"http://files.ntsg.umt.edu/data/NTSG_Products/MOD17/GeoTIFF/MOD17A3/"
        f"GeoTIFF_30arcsec/MOD17A3_Science_NPP_{MODIS_YEAR}.tif"
    )
    modis_da = xr.open_dataset(modis_url, engine="rasterio", chunks="auto")["band_data"].isel(band=0)

    print("Opening GLDAS porosity ...")
    gldas_da = xr.open_dataset(
        "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/gldas/SOILS/GLDASp5_porosity_025d.nc4",
        engine="h5netcdf", chunks="auto",
    )["GLDAS_porosity"].squeeze("time")

    print("Opening cropland zarr ...")
    da_crop = xr.open_zarr(S3_CROPLAND_ZARR)["cropland_frac"]

    # ── Region loop ────────────────────────────────────────────────────────────
    n_regions = len(gdf_regions)
    skipped   = []

    for i, row in gdf_regions.iterrows():
        reg_id           = int(row["reg_id"])
        region_name      = row["reg_nm"]
        region_nospaces  = row["region_nospaces"]
        region_geom      = row["geometry"]
        bounds           = region_geom.bounds

        print(f"\n[{i+1}/{n_regions}] {region_name} (id={reg_id})")

        if SKIP_EXISTING:
            s3_dir        = f"{S3_OUT_ROOT}/{region_nospaces}"
            expected_files = [
                f"{s3_dir}/spinup_vars.csv",
                f"{s3_dir}/cell_map.png",
                f"{s3_dir}/distance_map.png",
                f"{s3_dir}/gridded_data.zarr/.zgroup",
            ]
            if all(fs.exists(p) for p in expected_files):
                print("  → already complete, skipping")
                continue

        try:
            # ERA5
            ltm, lats, lons = era5_ltm_for_region(ds_era5, resolved_vars, bounds, YEAR_RANGE)
            gdf = build_cell_gdf(ltm, lats, lons, region_geom)
            if len(gdf) == 0:
                raise ValueError("No ERA5 cells within region polygon")

            # Cropland
            gdf["cropland_frac"] = sample_cropland(da_crop, gdf)
            if gdf["cropland_frac"].sum() <= 0:
                raise ValueError("No cropland weight in region")

            # MODIS, GLDAS
            gdf["MODIS_NPP_gC_m2_yr"] = sample_modis(modis_da, gdf)
            gdf["soil_porosity"]       = sample_gldas(gldas_da, gdf)

            # SoilGrids
            sg_df = sample_soilgrids(sg_ds, gdf)
            for col in sg_df.columns:
                gdf[col] = sg_df[col].values

            # ISRIC
            isric_df = sample_isric(isric_raster, pixel_to_newsuid, depth_tables, gdf)
            for col in isric_df.columns:
                gdf[col] = isric_df[col].values

            # Derived spinup vars
            gdf = derive_spinup_vars(gdf, RUNOFF_PROXY)

            # Representative cell
            gdf, rep, wt_mean, wt_std = find_representative_cell(gdf, REP_VARS, VAR_CONSTRAINTS)

            print(
                f"  rep cell: lat={rep['lat']:.4f}  lon={rep['lon']:.4f}  "
                f"cropland={rep['cropland_frac']:.3f}  dist={rep['rep_distance']:.3f}"
            )

            # Save
            save_region(fs, reg_id, region_name, region_nospaces, gdf, rep, wt_mean, wt_std, bounds)

        except Exception as exc:
            print(f"  SKIP: {exc}")
            traceback.print_exc()
            skipped.append((reg_id, region_name, str(exc)))

    # ── Summary ────────────────────────────────────────────────────────────────
    print(f"\n{'═'*60}")
    print(f"Done. Processed {n_regions - len(skipped)}/{n_regions} regions.")
    if skipped:
        print("Skipped:")
        for reg_id, name, reason in skipped:
            print(f"  id={reg_id}  {name}: {reason}")
    print(f"Results at: {S3_OUT_ROOT}/")

    # ── Compile all regions into a single parquet ───────────────────────────
    print("\nCompiling all-region parquet ...")
    compile_results(fs, S3_OUT_ROOT)


if __name__ == "__main__":
    main()

# %%
