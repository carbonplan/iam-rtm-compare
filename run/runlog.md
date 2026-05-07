
# Track batch files and progres

## Key questions
Does the richard's equation explain why seasonality gives you a base saturation pulse, but annual mean doesn't?
- lrOAT_noAmnt (v1.0.2 seasonality)
- lrOAT_ann (v1.0.2 annual mean)
- lrOAT_noAmnt_v101 (v1.0.1 seasonality)
- lrOAT_ann_v101 (v1.0.1 annual mean)

Does the buildup of slow wx minerals impact porosity differently from, say, gbas?
- lrOAT_ann_psdfull_baek23nodp (full feedstock)
- lrOAT_ann_psdfull (single feedstock)
- (Could also look at psdfull vs psdbulk for a multi-sp feedstock)

Does the base saturation pulse change the "lag times" between seasonal and annual runs?
- cao_ann (annual clim cao, roughly after Kanzaki et al. 2025)
- cao_szn (seasonal clim cao, roughly after Kanzaki et al. 2025)

Do we need constant application to maintain weathering rate?
- lrOAT_ann_noporo (annual applications for 100 years with same runs as 50on-50off)
- lrOAT_50on-50off (first 50 years are same as lrOAT_ann, then application stops)


## [ CLIM RESOLUTION ]
*GENERAL NOTES*  
<span style="color:red;">!! RE-RUN WITH UPDATED RUNOFF !! Earlier simulations didn't convert m to mm, and ERA5 runoff is crazy low anyway unless we use the Budyko correction from P/PET. </span>     

- We set an `EXP_GROUP = clim-resolution` such that output is on S3 in `../scepter_output/clim-resolution/..`   
- We seem to be having issues with some time resolutions failing. Not sure if there's an issue with the climate files or something else?   
- CDR outcomes appear to be pretty identical while the initial soil-bound acidity is neutralized --- then you can see divergence by time-resolution.   
### <span style="color:#F2668B;">clim-resolution_cc_sites_multi-clim_v0.csv</span>  
**STATUS**: (finnished running)      
**EXP_SET**: 'clim-resolution_cc_sites'  
    - v1.0.2; varying across timescale of climate forcing; 10 years  
    - PSD off; poro off; 
    - 3 sites (no central valley); 3 apprates; no fertrates  
### <span style="color:#F2668B;">clim-resolution_gbas_sites_multi-clim_v0.csv</span>  
**STATUS**: (finnished running)      
**EXP_SET**: 'clim-resolution_gbas_sites'  
    - as `clim-resolution_cc_sites` but for gbas      


## [ CEC simplerun ]
*GENERAL NOTES*
- We set an `EXP_GROUP = cec_simpleTitrate` such that output is on S3 in `../scepter_output/<EXP_GROUP>/..`   
- These are set up to save a series of json files which are then read in as the rundict (with no default dict) -- one json file per simulation in a batch. 
- We're looking at the base sat response to pH, the degree of alkalinity losses, and the cation-dependent timescales to breakthrough, and effect of alpha param.      

### <span style="color:#B23221;">caOnly</span>  
**STATUS**: (postproc complete)      
**EXP_SET**: 'caOnly'  
    - v1.0.2; seasonality off; 5 yrs    
    - inert soil profile with fixed CEC params;   
    - add water with dissolved Ca at varying NO3 levels to simulate pH variability   

### <span style="color:#B23221;">caOnly+cecvary</span>  
**STATUS**: (postproc complete)      
**EXP_SET**: 'caOnly+cecvary'  
    - same as `caOnly` but also varying alpha values.

### <span style="color:#B23221;">naOnly+cecvary</span>  
**STATUS**: (postproc complete)      
**EXP_SET**: 'naOnly+cecvary'  
    - same as `caOnly+cecvary` but using a fixed na+ concentration rather than ca2+ (at the same charge equivalent).   

## [ GCAM ]
*GENERAL NOTES*    
Create GCAM-derived SCEPTER setups by running `overlay-gcam+ag-map.ipynb`, then in era5 `create_clim_inputs_one-cell-per-region.ipynb`, then `calculate-rock-fluxes-from-GCAM.ipynb`, and finally `create-dust_temp-inputs.ipynb`.   

- R = running   
- Re = rerunning failed runs    
- Rc = run complete, needs postproc       
- Y = complete    

*Spinups*
| spinup | note                                          |
|--------|-----------------------------------------------|
| v1     | PmE runoff; pw pH 15-30cm;                    |


|                  | exe_central1p5 | exe_lowCDR1p5  | exe_noWeathering1p5 | exe_techBehaviorChange1p5 |
|------------------|----------------|----------------|---------------------|---------------------------|
| psdfull_nocec    | Y              | Y              | Y                   | Y                         |
| psdfull          | Y              | Y              | Y                   | Y                         |
| nopsd            | Y              |                |                     |                           |
| nopsd_10xrough   | Y              |                |                     |                           |
| psdfull_cornbelt | R              |                |                     |                           |
| max10tha         | Y              |                |                     |                           |
| 2xRock           | R              |                |                     |                           |
| nopsd_cao        | Re             |                |                     |                           |
| nopsd_cao_nocec  | (not run?)     |                |                     |                           |
| nopsd_smooth     | Y              |                |                     |                           |
| nopsd_3xcec      | Y              |                |                     |                           |
| low*             | Y              |                |                     |                           |

*smooth, 200 micron diam, 264 spinup with amnt, fickian mix

Table below is `gcam_opr2_*` (repeating under more controlled conditions)    

|                  | exe_central1p5 | exe_lowCDR1p5  | exe_noWeathering1p5 | exe_techBehaviorChange1p5 |
|------------------|----------------|----------------|---------------------|---------------------------|
| psdfull_nocec    | Y              |                | Y                   |                           |
| psdfull_lowcec   | Y              |                | Y                   |                           |
| psdfull_hicec    | Y              |                | Y                   |                           |
| psdfull_10xrough | Y              |                | Y                   |                           |
| psdfull_max5tha  | Y              |                | Y                   |                           |
| psdfull_midcec_10xrough  | Y              |                | Y                  |                           |
| nopsd            | Y              |                | Y                   |                           |
| nopsd_10xrough   | Y              |                | Re                  |                           |
| psdfull_5xrough  | Y              |                | Y                   |                           |
| psdfull_midcec_5xrough  | Y              |                | Y                   |                           |
| psdfull_midcec_5xrough_max5tha  | Y              |                | Y                   |                           |

Table below is `gcam_opr3_*` (repeating under more controlled conditions)    

|                  | exe_central1p5 | exe_lowCDR1p5  | exe_noWeathering1p5 | exe_techBehaviorChange1p5 |
|------------------|----------------|----------------|---------------------|---------------------------|
| psdfull_lowcec   | Y              |                | Y                   |                           |
| psdfull_hicec    | Y              |                | Y                   |                           |
| psdfull_smooth   | Y              |                | Y                   |                           |
| psdfull_hicec_smooth  | Y              |                | Y                   |                           |
| psdfull_loWx     | Y              |                | Y                   |                           |
| psdfull_25diam   | Y              |                | Y                   |                           |


## [ LAGS ]
*GENERAL NOTES*   
<span style="color:red;">!! RE-RUN WITH UPDATED RUNOFF !! Earlier simulations didn't convert m to mm, and ERA5 runoff is crazy low anyway unless we use the Budyko correction from P/PET. </span>     

- We set an `EXP_GROUP = lags` such that output is on S3 in `../scepter_output/lags/..`   
### <span style="color:#F2668B;">CARB VS SIL</span>  
**Note**: I think the ctrl run has no fertilizer, so fertilizer runs have an extra apparent emission due to the effect of adding fertilizer relative to no-fert counterfactual.   
### <span style="color:#F2668B;">lags_gbas_sites_ann_ltm_v0.csv</span>  
**STATUS**: (postproc complete)      
**EXP_SET**: 'lags_gbas_sites'  
    - v1.0.2; seasonality off; 100 years  
    - PSD off; poro off; 
    - 3 sites (no central valley); 3 apprates; 3 fertrates  

### <span style="color:#F2668B;">lags_cc_sites_ann_ltm_v0.csv</span>  
**STATUS**: (postproc complete)      
**EXP_SET**: 'lags_cc_sites'  
    - as `lags_gbas_sites` but for calcite feedstock   

### <span style="color:#F2668B;">cecOnOff_\*_\*_ann_ltm_v0.csv</span>  
**STATUS**: (run complete)      
**EXP_SET**: 'cecOnOff_\*_\*'  
    - Four total short batches for combos of cec On/Off and gbas vs cc   
    - 100 yr, no seasonality, no psd, no poro    
    - northeast and cornbelt sites    



### <span style="color:#05F6F2;">SEASONAL VS ANNUAL</span>  
### <span style="color:#05F6F2;">cao_ann_monthly_ltm_v0.csv</span>  
**STATUS**: (postproc complete)  
**EXP_SET**: 'cao_ann'  
    - annual climate  
    - cao added uniformly to mixdep for 90% of year  
    - varies across dustrate, secondary mineral choice   


### <span style="color:#05F6F2;">cao_szn_monthly_ltm_v0.csv</span>  
**STATUS**: (postproc complete)  
**EXP_SET**: 'cao_szn'  
    - seasonal climate  
    - cao added uniformly to mixdep for 90% of year  


### <span style="color:#eeff41;">COLUMN DEPTH</span>  
### <span style="color:#eeff41;">lags_depthvary_5app_ann_ltm_5yrs_v0.csv</span>  
**STATUS**: (running)  
**EXP_SET**: 'lags_depthvary_5app'  
    - annual climate  
    - cao added annually for 5 years (2 t/ha/yr), then 100 years of nothing  
    - varies across soil depth (based on spinup) and CEC (from 5-25)
    - no fert   

### <span style="color:#eeff41;">lags_depthvary_1app_ann_ltm_1yrs_v0.csv</span>  
**STATUS**: (running)  
**EXP_SET**: 'lags_depthvary_1app'  
    - as `lags_depthvary_5app` but only apply cao in the first year and then never again  


### <span style="color:#66ff00;">SINGLE APP HYSTERESIS</span>  
### <span style="color:#66ff00;">lags_singleApp_cc_sites_ann_ltm_v0.csv</span>  
**STATUS**: (running)  
**EXP_SET**: 'lags_singleApp_cc_sites'  
    - check for same breakthrough - hysteresis curves we see in Kate's runs    
    - annual clim; calcite added in year 1, then nothing for 14 years  
    - varies across rock app amount and fert rate   

### <span style="color:#66ff00;">lags_singleApp_singlerun_gbas_sites_ann_ltm_v0.csv</span>  
**STATUS**: (running)  
**EXP_SET**: 'lags_singleApp_singlerun_gbas_sites'  
    - check for same breakthrough - hysteresis curves we see in Kate's runs    
    - annual clim; gbas added in year 1, then nothing for 14 years  
    - varies across rock app amount and fert rate   


## [ SURFACE AREA ]
*GENERAL NOTES*
- Exploring: (1) different outcomes from using grain size vs rock amount as the knob for varying surface area; (2) Differences in removal outcomes associated with uncertainties in the roughness factor and PSD parameterization (e.g., SA evol over time) (3) Differences in how surface area is applied over time (single app, multiple smaller apps)    
- can focus mostly on gbas here since cc is so insensitive      

### varyRough
pattern: `surface_area_fullPSD_varyRough_rad<XXX>um_ann_ltm_<FS>_v0.csv`
example: `surface_area_fullPSD_varyRough_rad100um_ann_ltm_baek23_v0.csv`
> EXP_GROUP = `"surface_area"`

> EXP_SET = `f"fullPSD_varyRough_rad{int(dustrad_init)}um"`
> EXP_SET = `f"fullPSD-noCEC_varyRough_rad{int(dustrad_init)}um"`

- R = running   
- Re = rerunning failed runs    
- Rc = run complete, needs postproc       
- Y = complete

| (fullPSD)          | gbas | cc  | bridge | baek23 |
|--------------------|------|-----|--------|--------|
| dustrad_init = 25  | Y    | Y   | Y      | Y      |
| dustrad_init = 100 | Y    | Y   | Y      | Y      |
| dustrad_init = 300 | Y    | Y   | Y      | R      |
| dustrad_init = 600 | Y    | Y   | Y      | R      |

| (fullPSD-noCEC)    | gbas | cc  | bridge | baek23 |
|--------------------|------|-----|--------|--------|
| dustrad_init = 25  | Y    | Rc  | Rc     | Re     |
| dustrad_init = 100 | Y    | Rc  | Rc     | Re     |
| dustrad_init = 300 | Rc   | Rc  | Rc     | Re     |
| dustrad_init = 600 | Rc   | Rc  | Rc     | Re     |

### <span style="color:#EDAB54;">surface_area_fixedSA_noPSD_sites_ann_ltm_v0.csv</span>  
**STATUS**: (run complete)  
**EXP_GROUP**: 'surface_area'
**EXP_SET**: 'fixedSA_noPSD_sites'  
    - runs using the same initial surface area but different grain size - rock amount combinations.   
    - No PSD_full; No PSD_bulk; w/ roughness;   
    - vary across particle radius for a given target sa (collectively used to compute a single rock amount)   

### <span style="color:#EDAB54;">surface_area_fixedSA_fullPSD_sites_ann_ltm_v0.csv</span>  
**STATUS**: (tried running with bulkPSD but got "p80 not found" error... haven't figured out PSD with scepter_richards I guess..)  
**EXP_GROUP**: 'surface_area'
**EXP_SET**: 'fixedSA_fullPSD_sites'  
    - runs using the same initial surface area but different grain size - rock amount combinations.   
    - PSD_full (PSD_bulk doesn't work in v1.0.2); w/ roughness;   
    - vary across particle radius for a given target sa (collectively used to compute a single rock amount)  



## [ LONGRUN ]
*GENERAL NOTES*
- OAT stands for one-at-a-time (rather than all combinations of all parameters we perturb, we just change things one-at-a-time from a control simulation)


### <span style="color:#FFF59D;">lrOAT_monthly_ltm_v0.csv</span>  
**STATUS**: running(3) + postproc | frontsDat     
**EXP_SET**: 'lrOAT'  
    - v1.0.2 ; seasonality on  
    - PSD off  
    - amnt added via `OM_rain.in`  
    - multiyear = False  
    - 4 sites sent from Yoshi (central valley; cornbelt; SE; NE)   


### <span style="color:#FFF59D;">lrOAT_noAmnt_monthly_ltm_v0.csv</span>  
**STATUS**: running(2) + postproc | frontsDat   
**EXP_SET**: 'lrOAT_noAmnt'  
    - PSD off  
    - `spinup_tag = "_no-amnt"` uses spinups with no amnt in `OM_rain.in`  
    - otherwise same as `lrOAT`  


### <span style="color:#FFF59D;">lrOAT_noAmnt_v101_monthly_ltm_v0.csv</span>  
**STATUS**: postproc(58)     
**EXP_SET**: 'lrOAT_noAmnt_v101'  
    - PSD off  
    - v101 so we don't have Richard's equation   
    - `spinup_tag = "_no-amnt"` uses spinups with no amnt in `OM_rain.in`  
    - otherwise same as `lrOAT`  


### <span style="color:#FFF59D;">lrOAT_noAmnt_cbas_monthly_ltm_v0.csv</span>  
**STATUS**: postproc | frontData   
**EXP_SET**: 'lrOAT_noAmnt_cbas'  
    - PSD off  
    - cbas instead of gbas (crystalline instead of glassy)  
    - `spinup_tag = "_no-amnt"` uses spinups with no amnt in `OM_rain.in`  
    - otherwise same as `lrOAT`  


### <span style="color:#FFF59D;">lrOAT_psdbulk_monthly_ltm_v0.csv</span>  
**STATUS**: (Ready for full postproc) | running(17) + postproc    
**EXP_SET**: 'lrOAT_psdbulk'  
    - PSDbulk on  
    - otherwise same as `lrOAT_noAmnt`  


### <span style="color:#FFF59D;">lrOAT_psdfull_monthly_ltm_v0.csv</span>  
**STATUS**: running(14) + postproc  
**EXP_SET**: 'lrOAT_psdfull'  
    - PSDfull on  
    - otherwise same as `lrOAT_noAmnt`  


### <span style="color:#FFF59D;">lrOAT_saRule1_monthly_ltm_v0.csv</span>  
**STATUS**: running(1) + postproc | frontsDat   
**EXP_SET**: 'lrOAT_saRule1'  
    - SA rule 1 is on; 2 is off; (others used SA rule 2)  
    - otherwise same as `lrOAT_noAmnt`  


### <span style="color:#FFF59D;">lrOAT_ann_monthly_ltm_v0.csv</span>  
**STATUS**: postproc   
**EXP_SET**: 'lrOAT_ann'  
    - seasonality turned off    
    - otherwise same as `lrOAT_noAmnt` 


### <span style="color:#FFF59D;">lrOAT_ann_v101_monthly_ltm_v0.csv</span>  
**STATUS**: running   
**EXP_SET**: 'lrOAT_ann_v101'  
    - seasonality turned off  
    - v101 so we don't have richard's eq'n     
    - otherwise same as `lrOAT_noAmnt` 


### <span style="color:#FFF59D;">lrOAT_ann_cbas_monthly_ltm_v0.csv</span>  
**STATUS**: postproc   
**EXP_SET**: 'lrOAT_ann_cbas'  
    - seasonality turned off   
    - cbas instead of gbas (crystalline instead of glassy)   
    - otherwise same as `lrOAT_noAmnt` 


### <span style="color:#FFF59D;">lrOAT_ann_psdfull_monthly_ltm_v0.csv</span>  
**STATUS**: processed   
**EXP_SET**: 'lrOAT_ann_psdfull'  
    - seasonality turned off    
    - psdfull turned on  
    - otherwise same as `lrOAT_noAmnt`  


### <span style="color:#FFF59D;">lrOAT_ann_psdfull_baek23nodp_monthly_ltm_v0.csv</span>  
**STATUS**: processed(21)   
**EXP_SET**: 'lrOAT_ann_psdfull_baek23nodp'  
    - seasonality turned off    
    - psdfull turned on  
    - baek23nodp feedstock (to look at buildup of slow wx mins)
    - otherwise same as `lrOAT_noAmnt` 


### <span style="color:#FFF59D;">lrOAT_ann_psdbulk_monthly_ltm_v0.csv</span>  
**STATUS**: [FAILED]   
**EXP_SET**: 'lrOAT_ann_psdbulk'  
    - seasonality turned off  
    - psdbulk turned on    
    - otherwise same as `lrOAT_noAmnt` 


### <span style="color:#FFF59D;">lrOAT_fickmix_monthly_ltm_v0.csv</span>  
**STATUS**: postproc  
**EXP_SET**: 'lrOAT_fickmix'  
    - `imix`=1 fickian mixing (rather than 3, tilling)  
    - seasonality turned off    
    - otherwise same as `lrOAT_noAmnt` 


### <span style="color:#FFF59D;">lrOAT_hommix_monthly_ltm_v0.csv</span>  
**STATUS**: postproc  
**EXP_SET**: 'lrOAT_hommix'  
    - `imix`=2 homogeneous mixing (rather than 3, tilling)  
    - seasonality turned off    
    - otherwise same as `lrOAT_noAmnt` 


### <span style="color:#FFF59D;">lrOAT_ann_noporo_monthly_ltm_v0.csv</span>  
**STATUS**: running   
**EXP_SET**: 'lrOAT_ann_noporo'   
    - same as `lrOAT_ann` but no porosity variability    


### <span style="color:#FFF59D;">lrOAT_50on-50off_monthly_ltm_v0.csv</span>  
**STATUS**: postproc(3)   
**EXP_SET**: 'lrOAT_50on-50off'      
    - same as `lrOAT_ann_noporo` but we use multiyear mode to only apply rock for the first 50 years, then nothing the next 50  


### <span style="color:#FFF59D;">lrOAT_100on-50off_monthly_ltm_v0.csv</span>  
**STATUS**: postproc   
**EXP_SET**: 'lrOAT_100on-50off'      
    - noporo, no PSD, no seasonality (not sure if there's an equivalent other thing)     


### <span style="color:#FFF59D;">pr1_ann_monthly_ltm_v0.csv</span>  
**STATUS**: postproc  
**EXP_SET**: 'pr1_ann'  
    - seasonality turned off  
    - PSD off     
    - uses spinup with full mineralogy (rather than solely inert)
    - one site, single cliamte case


### <span style="color:#FFF59D;">pr1_ann2_monthly_ltm_v0.csv</span>  
**STATUS**: postproc  
**EXP_SET**: 'pr1_ann2'  
    - as `pr1_ann` but never changing porosity (since this upsets the bulksoil calculations)  


### <span style="color:#FFF59D;">pr1_ann2_cbas_monthly_ltm_v0.csv</span>  
**STATUS**: postproc  
**EXP_SET**: 'pr1_ann2_cbas'  
    - as `pr1_ann2` with cbas basalt feedstock


### <span style="color:#FFF59D;">pr1_ann2_cc_monthly_ltm_v0.csv</span>  
**STATUS**: postproc (with some missing)  
**EXP_SET**: 'pr1_ann2_cc'  
    - as `pr1_ann2` with calcite feedstock  


### <span style="color:#FFF59D;">pr1_ann2_baek23nodp_monthly_ltm_v0.csv</span>  
**STATUS**: postproc  
**EXP_SET**: 'pr1_ann2_baek23nodp'  
    - as `pr1_ann2` with baek23nodp basalt feedstock  


### <span style="color:#FFF59D;">pr9_ann2_monthly_ltm_v0.csv</span>  
**STATUS**: postproc  
**EXP_SET**: 'pr9_ann2'  
    - as `pr1_ann` but never changing porosity (since this upsets the bulksoil calculations) and using pr9 spinup    




---


## [ LONGRUN - deprecated ]
*GENERAL NOTES*
- OAT stands for one-at-a-time (rather than all combinations of all parameters we perturb, we just change things one-at-a-time from a control simulation)


### <span style="color:#FFF59D;">pr1_ann2_bridge_monthly_ltm_v0.csv</span>  
**STATUS**: running | RUNS FAILED  
**EXP_SET**: 'pr1_ann2_bridge'  
    - as `pr1_ann2` with blue ridge basalt feedstock  


*longrunOAT_v2_monthly_ltm_annual_v0.csv*
STATUS: first-round running
EXP_SET: 'longrunOAT_v2'
    - v1.0.2 ; seasonality on
    - multiyear = False
    - 4 sites sent from Yoshi (central valley; cornbelt; SE; NE)
    - (note `longrunOAT` had errors, consider it deprecated)

*longrunOAT_noAmnt_monthly_ltm_annual_v0.csv*
STATUS: first-round running
EXP_SET: 'longrunOAT_noAmnt'
    - v1.0.2 ; seasonality on 
    - Same as `longrunOAT_monthly_ltm_annual_v0.csv` but "amnt" is removed from `OM_rain.in`

*longrunOAT_annClim_monthly_ltm_annual_v0.csv*
STATUS: first-round complete
EXP_SET: 'longrunOAT_annClim'
    - v1.0.2 ; seasonality off
    - otherwise same as `longrunOAT_monthly_ltm_annual_v0.csv`

*longrunOAT_oldSCEPT_monthly_ltm_annual_v2.csv*
STATUS: first-round complete
EXP_SET: 'longrunOAT_oldSCEPT'
    - v1.0.1 ; seasonality off
    - otherwise same as `longrunOAT_monthly_ltm_annual_v0.csv`
