# SOLAR FEWLS TOOL 

Last Updated: 05/02/2023
Primary Contact: Jacob Stid (stidjaco@msu.edu)

Release Notes: 
* February 18, 2025: Included crop_values_CA_05_18.csv in Data/Downloaded and all scenario outputs in Outputs.
* October 21, 2023: Initial release.


General Notes: 
This is the general information README document for code, data, and logic behind the Solar Food, Energy, Water, and Economic Lifecycle Scenario (FEWLS) Tool. 
The FEWLS tool was developed and used in the research article, **From Fields to Photovoltaics: Effects of Agrisolar Co-Location on Food, Energy, Water, and 
Economic Security**. In general, this code takes in a ground-mounted solar PV shape file (with some auxiliary information described below) and generates 
lifespan predictions for food, energy, water, and economic effects due to offsetting agricultural land with solar PV energy generation 
(*agrisolar co-location*). The code and data repository can be found at the github/stidjaco/FEWLS_tool page, the Zenodo Repository 
(https://doi.org/10.5281/zenodo.10023281), or can be shared upon request. 


Running FEWLS: 
To run the default model with default variables, datasets, and attributes, simply run the R script, *FEWLS_run.R*. This script calls all necessary scripts, 
functions, packages, and operations and automatically outputs result dataframes and figures to *FEWLS_tool/Outputs*. To run FEWLS for a new set of solar 
installations, ensure that all attributes in the below section *Solar PV Spatiotemporally Characterized Dataset* are included. Other variables can be altered
in *FEWLS_inputs.R* and/or in *Downloaded* or *Derived* datasets. FEWLS has not yet been stressed tested in other regions, and may require code alteration to 
be operable elsewhere. 


Folder level annotation:
* "#" indicates top level directory
* "##" indicates second level directory
* "###" indicates third level directory
* "####" indicates fourth level directory
* Deeper directories are indicated within "####" using "-" to delineate depth. 


Acronyms: 
* FEWLS: Food, Energy, Water, and Economic Lifecycle Scenario
* PV: photovoltaic
* CDL: Cropland Data Layer
* NEM: Net Energy Metering
* NSRDB: National Solar Radiation Database
* GEE: Google Earth Engine
* CCV: California's Central Valley
* USDA: United States Department of Agriculture
* NASS: National Agriculture Statistics Service
* LCMAP: Land Change Monitoring, Assessment, and Projection
* NREL: National Renewable Energy Laboratory
* ATB: Annual Technology Baseline
* CAPEX: Capital Expenditure
* O&M: Operation and Maintenance
* FRIS: Farm and Ranch Irrigation Survey (2013)
* IWMS: Irrigation and Water Management Survey (2018)
* DTW: Depth to water
* PG&E: Pacific Gas & Electric 
* SMUD: Sacramento Municipal Utility District
* SCE: Southern California Edison


Solar PV Spatiotemporally Characterized Dataset:
This shape file is the driving data behind the FEWLS model. It is stored in *FEWLS_tool/CV_Agrisolar_df/CV_Agrisolar_df.shp*, and has been 
pre-processed to include arrays present in Kruitwagen et al. (2021) that were omitted in our early work, Stid et al. (2022), and the essential attributes below
using the *FEWLS_tool/GEE_scripts*. It contains 935 solar PV arrays with the following essential attributes for the FEWLS tool: 
* Index: Unique dataset array index for reference throughout FEWLS tool - [Integer]
* Class: Mount type (fixed or south oriented, single tracking or e/w oriented) - [String]
* dir_a: Total area taken up by panels and space between for single system - [square meters]
* pnl_a: Total area of panels within single system (proximal PV panels of the same mount type installed in the same year) - [square meters]
* sys_C_#: The USDA Cropland Data Layer mode value for the array "total area" coverage for "#" year(s) prior to the year of installation. This is assumed to be the foregone croptype for the array and is acquired from *LC_RP_MultiYear.txt* & *PV_Stid_Omit_prep.txt* - [CDL Crop Values]
* Crp_n_#: The crop name which corresponds to the crop value described by the Cropland Data Layer, this is acquired by relating sys_C_# with NASS_Classifications.csv "Values" and "Crop" - [string]
* Capacity: Capacity given temporal panel efficiency, panel area, packing factor, and latitude (as described by Martin-Chivelet [2015]) - [MW]
* Year: Estimated solar PV installation year based on method from Stid et al. (2022) and *YOD_LT_get.txt* & *YOD_Manual_Interpretation.txt*. - [integer]
* Gen_20##: pvlib modeled electricity generation for a given system during a given year - [MWh]
* Econ_20##: modeled NEM based returns for electrictiy produced. Incporporates differences in NEM rate strucutres (depending on "Year") and returns for offset (respective rate) and surplus credit (depends on NEM guidelines). This varible must exist but should be "0" for all arrays and year not included in NEM (<1MW). In FEWLS_model.R, you can also select the getResource_energySingleRate(), to use elec_price variable in FEWLS_inputs.R if modeled rates aer not available. See the manuscript and supplemental information for more information. - [$USD nominal] 
* Source: Source for array geospatial coordinates, either "Stid" or "Kruitwagen" - [String]
* irrig: A binary indication of whether or not the direct area of an array was irrigated the year prior to installation using the LanID irrigated area dataset (Xie et al. [2021]), and the lanid_irrigation_export.txt GEE script - [binary integer: 1 = irrigated, = 0 non-irrigated]
* geometry: WGS84 geometry of array "total direct area" encompassed by panels (PV_ID_panels.shp) of the same mount technology and same installation year proximal to one another. Full characterization and grouping methods defined in Stid et al. (2022). - [geometry]
As long as an auxiliary dataset contains the same formatting, column nomenclature, and essential attributes, the tool can be used. 







# FEWLS_tool
// Description: Parent directory for all code and data. FEWLS_tool will automatically call the working directory for this file for automated processes. Detailed descriptions start at line 10 of each file.
* Files: 
  * FEWLS_run.R: R script to compile and execute FEWLS model. All other (necessary) scripts are called here. 
  * FEWLS_inputs.R: R script containing FEWLS model input variables to be used in the model. This includes all required file locations (without which, FEWLS cannot run). 
  * FEWLS_model.R: R script containing all model functions and operations. This file itself does not execute the model, only compiles model operations. 
  * FEWLS_ResourceModelResults.R: R script to execute generation of *FEWLS_tool/Outputs* result dataframes and figures for agrisolar food, energy, and water resource effects. 
  * FEWLS_EconModelResults.R: R script to execute generation of *FEWLS_tool/Outputs* result dataframes and figures for agrisolar food, energy, and water resource effects. 
  * FEWLS_NonUserInputs.R: R script to save within-model variables for per-array model runs. This file should not be altered for model runs.
  * FEWLS_mmPerYrIrrig.R: R script to generate *FEWLS_tool/Data/Derived/mm_per_yearIrrig.csv*, a dataframe containing irrigation depth for crops from the FRIS (2013) and IWMS (2018). 
  * FEWLS_cropRevenue.R: R script to generate *FEWLS_tool/Data/Derived/crop_revenue_kg.csv*, a dataframe containing $/kg prices recieved for CDL crop types derived from USDA NASS Monthly Agricultural Prices Report (see crop_values_CA_06_18.csv). 
  * FEWLS_cropCost.R: R script to generate *FEWLS_tool/Data/Derived/operational_cost_m2.csv*, a dataframe containing operational crop costs from UC Davis Crop Ext. 
  * FEWLS_irrigEnergyReq.R: R script to generate *FEWLS_tool/Data/Derived/irrigEnergyReq_adjusted.csv*, a dataframe of county level irrigation energy requirements (GWh/m3).
  * FEWLS_kcalperm2.R: R script to generate *FEWLS_tool/Data/Derived/kcal_per_m2.csv*, a dataframe of crop specific caloric density and spatial caloric density from USDA NASS Agricultural Yield Surveys and FoodData nutrition values for the 2022 report. Acquired via [USDA QuickStats service](https://quickstats.nass.usda.gov/) and [FoodData](https://fdc.nal.usda.gov/)
  * FEWLS_FarmElecBudgExport.R: R script to generate three dataframes, *FEWLS_tool/Data/Derived/SolarFEWE_generation_input.csv*, *FEWLS_tool/Data/Derived/SolarFEWE_LoadDf_annual.csv*, and *FEWLS_tool/Data/Derived/SolarFEWE_econ_inputs.csv*, and one figure *FEWLS_tool/Outputs/Figures/Surplus_Deficits##LS_##MW.pdf*, supplemental dataframes and figure used in the Agrisolar manuscript and described below in *Data/Derived*. 
  * FEWLS_supplemental.R: R script to generate supplemental figures and all other processes intended for the manuscript, not necessarily essential for the model. If something is missing or unclear from the present methods or files, check here. 
  * ElectricityRateScrape.RMD: RMarkdown file generating hourly electric rate matrices used in net billing and net metering analysis for the six different utilities assessed here.  







## CV_Agrisolar_df
// Description: Folder containing the input solar PV installation shape file
* Files: 
 * CV_Agrisolar_df.shp: The driving data behind the FEWLS model. It contains 935 solar PV arrays with the spatiotemporal attributes described above in *Solar PV Spatiotemporally Characterized Dataset*.







## FEWLS_tool/GEE_scripts
// Description: Sub-directory for all Google Earth Engine scripts used in the preprocessing for FEWLS_tool. 
// This can also be accessed via the shared GEE repository: [Photovoltaics GEE Repo](https://code.earthengine.google.com/?accept_repo=users/stidjaco/PhotoVoltaics)
* Files: 
  * Active_AgArea.txt: JavaScript GEE code to estimate annual active agricultural land area from CDL in CCV. Exports dataframe of year and active ag area.
  * AdjacentLCC_PostInstall.txt: JavaScript GEE code to acquire proximal landcover change prior to and post solar installation. 
  * FallowIdle_area.txt: JavaScript GEE code to estimate annual fallowed agricultural land area from CDL in CCV. Exports dataframe of year and fallow ag area.
  * gridMET_get.txt: JavaScript GEE code to acquire avg county level annual precipitation (mm) from the Idaho EPSCOR GridMET product. Used in irrigation estimations. 
  * lanid_irrigation_export.txt: JavaScript GEE code to acquire information on whether the prior agricultural land use for solar PV installations was irrigated or not. 
  * LC_RP_MultiYear.txt: JavaScript GEE code to acquire CDL and LCMAP land cover types for 'n' years prior to solar PV installation.
  * PV_Stid_Omit_prep.txt: JavaScript GEE code to perform spatiotemporal characterization of arrays in Kruitwagen et al. (2021) that were omitted from Stid et al. (2022). 
  * YOD_LT_get.txt: JavaScript GEE source code to perform spatiotemporal characterization of all arrays. Employs LandTrendr. Functions are exported. 
  * YOD_Manual_Interpretation.txt: JavaScript GEE code to manually interpret omitted installations years from LandTrendr analysis using available satellite/aerial imagery. 







## FEWLS_tool/pvlib_Analysis
// Description: Sub-directory for data and code used to acquire modeled electricity generation and NEM returns for commercial-scale solar PV. Results in output dataframes post-processed in FEWLS_supplemental.R, and included in the base requried in the Solar PV Spatiotemporally Characterized Dataset (in_solar_df.shp).
* Files: 
  * AgrisolarGenNEM.py: Python script that exports matrices for each year of either electricity generation or NEM returns (nominal $). The script is currently coded to return NEM returns, which inherently calculate annual generation. Outputs files to FEWLS_pvlibResults folder


### FEWLS_tool/pvlib_Analysis/FEWLS_pvlibResults
// Description: Sub-directory for all solar electricity generation and NEM return data. 
* Files: 
  * economic_value_20##.csv: NEM return matrix for solar PV arrays installed in the year 20##. In nominal USD, based on reported utility rates derived in ElectricityRateScrape.RMD
  * generation_value_20##.csv: Electricity generation matrix for solar PV arrays installed in the year 20## in kWh. 


### FEWLS_tool/pvlib_Analysis/Weather_files
// Description: Sub-directory for all solar electricity generation and NEM return data. 
* Folders: 
  * 20##: Folder contianing weather files for the year 20## from NSRDB at ## temporal resolution. Coded for location and called in AgrisolarGenNEM.py based on proximity. 


### FEWLS_tool/pvlib_Analysis/utility_rates
// Description: Sub-directory for all solar electricity generation and NEM return data. 
* Files: 
  * 20##: Folder contianing utility-specific fixed (customer), demand, and energy-charge electricity rates for the year 20## from a variety of sources. See *FEWLS_tool/Data/Derived/UtilRates* for more info. 







## FEWLS_tool/Data
// Description: Sub-directory for all data used in FEWLS_tool. 


### FEWLS_tool/Data/Downloaded
// Description: Sub-directory for all downloaded data used in FEWLS_tool. 
* Files: 
  * CPI.csv: Consumer Price Index dataframe for inflation adjustment from [US Bureau of Labor Statistics](https://www.bls.gov/cpi/)
  * USDA_yield_08to18_state.csv: Dataframe containing crop specific annual yields (yield unit / unit area) from USDA NASS Agricultural Yield Surveys. Acquired via [USDA QuickStats service](https://quickstats.nass.usda.gov/).
  * crop_values_CA_05_18.csv: Dataframe containing crop specific annual price received ($ / yield unit) from USDA NASS Monthly Agricultural Prices Report. Acquired via [USDA QuickStats service](https://quickstats.nass.usda.gov/).
  * efficiency_yearly_LBerkeley.csv: Dataframe containing annual efficiency and monocrystalline shares in the US solar market from the Lawrence Berkely Tracking the Sun report. 
  * IrrigWU_CONUScounty85to15.csv: Dataframe containing state and county level annual irrigation volumes for surface water, groundwater, and total irrigation water withdrawls. Acquired via [USDA QuickStats service](https://quickstats.nass.usda.gov/).
  * USDA_FARS_1318.csv: California specific Farm and Ranch Irrigation and Irrigation and Water Management Survey irrigation water depths applied (acrefeet/acre). Acquired via [USDA QuickStats service](https://quickstats.nass.usda.gov/).
  * USDA_FARS_1318_CONUS.csv: United States Farm and Ranch Irrigation and Irrigation and Water Management Survey irrigation water depths applied (acrefeet/acre). This is used to fill crop gaps in the California Dataset, with preference given to California. Acquired via [USDA QuickStats service](https://quickstats.nass.usda.gov/).
  * CONUS_aglandarea_operations.csv: Dataframe containing 5-year average farm sizes (acres) at the county level. Acquired via [USDA QuickStats service](https://quickstats.nass.usda.gov/).
  * NREL_ATB_historical_all.csv: Dataframe containing projected CAPEX O&M costs from [NREL's Electricity Annual Technology Baseline (ATB) report](https://data.openei.org/submissions/5865).
  * NREL_ATB_Historical_Trends_data.csv: Dataframe containing historical ($/kW-DC) median, 80th, and 20th percentile installation costs from [NREL's Electricity Annual Technology Baseline (ATB) report](https://data.openei.org/submissions/5865).
  * CONUS_aglandarea_operations.csv: Dataframe containing county-level agricultural area and number of operations from the 2007, 2012, and 2017 Census of Agriculture. Acquired via [USDA QuickStats service](https://quickstats.nass.usda.gov/).
  * County_IrrigEnergyInputsMcCarthy.csv: Dataframe containing result outputs from McCarthy et al., (2022), to be input and reprocessed with a better CCV depth to water product in FEWLS_irrigEnergyReq.R. This work is in preparation for publishing. 

#### FEWLS_tool/Data/Downloaded/CA_UtilService_shp
// Description: Sub-directory for the California Electric Utility Service Areas Shape file. 
* Files: 
  * California_Electric_Utility_Service_Areas.shp: California Electric Utility Service Areas Shape files downloaded from [CA.gov](https://gis.data.ca.gov/datasets/7d9a1dc059a74bd7a4fd89ff2a78a3bf_0/explore).


#### FEWLS_tool/Data/Downloaded/CV_alluvial_bnd
// Description: Sub-directory for for CCV alluvial boundary from Faunt et al. (2009). 
* Files: 
  * Alluvial_Bnd.shp: CCV alluvial boundary from Faunt et al. (2009). 


#### FEWLS_tool/Data/Downloaded/Depth_to_Water
// Description: Sub-directory for old or preliminary Depth to Water raster products intended for the irrigation energy requirement analysis. 
* Folders: 
  * SGMA: 
    * Files: i08_GroundwaterDepthSeasonal_Contours.shp: Downloaded spring 2018 DTW contours for CCV used to acquire DTW raster in FEWLS_tool/Data/Derived/DTW_raster from [SGMA Data Viewer](https://data.cnra.ca.gov/dataset/periodic-groundwater-level-measurements https://sgma.water.ca.gov/webgis/?appid=SGMADataViewer{\#}currentconditions)
  * Zell_Sanford_2020:
    * Files: conus_MF6_SS_Unconfined_250_dtw.tif: DTW raster downloaded from Zell & Sanford (2020). Old product no longer used in analysis for lack of recently relevant DTW in CCV due to recent drought and overdraft. 


#### FEWLS_tool/Data/Downloaded/FoodData
// Description: Sub-directory for FoodData repository download from [FoodData](https://fdc.nal.usda.gov/). This folder contains all downloaded FoodData database files, only a few of which are used in FEWLS. 
* Files (used in analysis): The R script, *FEWLS_kcalperm2.R* scrapes and prepares the final derived dataframe used in the analysis, *FEWLS_tool/Data/Derived/kcal_per_m2.csv*.
  * food.csv: Dataframe containing food names, units, ID's, and nutrient ranks (kcal is ID = 2047 | 2048 | 1008).
  * food_nutrient.csv: Dataframe containing food ID's, and nutrient values. All kcal values are per 100 grams of food


#### FEWLS_tool/Data/Downloaded/UCDavis_CropExt
// Description: Sub-directory for the California UC Davis Crop Extension Cost and Return Studies from [Current Cost & Return Studies](https://coststudies.ucdavis.edu/en/current/). 
* Files: 
  * 79 total .pdf reports containing crop specific operational costs. The R script, *FEWLS_cropCost.R* scrapes and prepares the final derived dataframe used in the analysis, *FEWLS_tool/Data/Derived/operational_cost_m2.csv*.


#### FEWLS_tool/Data/Downloaded/UtilityRate_Reports
// Description: Sub-directory for downloaded historical utility electricity rates for three of six utilities assessed. Other employed U.S. Utility Rate Database (links and detailed README's in FEWLS_tool/Data/Derived/UtilRates)
* Folders: 
  * Pacific gas and electric company: 2008 to 2018 annual csv's containing PG&E electricity rates delineated by large ag (LgAg), small ag (SmAg), and effective report dates, within file titles. 
  * Sacramento Municipal Utility District: Contains .pdf files that are annual (2010, 2011, 2013, 2015, 2017) General Managers reports for SMUD. 
  * Southern California Edison: Contains sub-folders for each year, containing .pdf's delineating SCE schedules for a given year.


### FEWLS_tool/Data/Derived
// Description: Sub-directory for all derived data used in FEWLS_tool. 
* Files: 
  * mm_per_yearIrrig.csv: Output of *FEWLS_mmPerYrIrrig.R". Default file will be present unless *cleanFEWLS()* is run. State level irrigation depths (mm/yr) for 2013 and 2018 survey years
  * crop_revenue.csv: Output of *FEWLS_cropRevenue.R*. Default file will be present unless *cleanFEWLS()* is run. State level price recieved ($/kg) for CDL crops from NASS Agricultural Prices Reports.
  * irrigEnergyReq_adjusted.csv: Output of *FEWLS_irrigEnergyReq.R*. Default file will be present unless *cleanFEWLS()* is run. County level irrigation energy requirements (GWh/m3) based on DTW raster. 
  * operational_cost_m2: Output of *FEWLS_cropCost.R*. Default file will be present unless *cleanFEWLS()* is run. Crop level operational costs ($/m2) based on UC Davis Crop Ext Cost Studies. 
  * kcal_per_m2: Output of *FEWLS_kcalperm2.R*. Default file will be present unless *cleanFEWLS()* is run. Crop level spatial caloric density (kcal/m2) based on FoodData database and NASS Agricultural Yield Surveys. 
  * NASS_Classifications.csv: NASS CDL crop types, with assumed Farm and Ranch Irrigation Survey and Irrigation and Water Management Survey crop type comparisons. 
  * county_annu_gridMET.csv: GridMET derived annual average county-level precipitation depths (mm) derived from *gridMET_get.txt*.
  * SolarFEWE_generation_input.csv: Dataframe containing the system index, panel area, mount class, installation year, pvlib modeled annual electricity generated by the solar array (MWh/yr)
  * SolarFEWE_econ_inputs.csv: Dataframe containing the system index, panel area, mount class, installation year, and the estimated monthly load for the farm associated with each array (MWh/month) ** this is distributed by the load profiles we had discussed. 
  * SolarFEWE_LoadDF_annual.csv: Dataframe containing three rows for each installation delineating annual load, deficit, and surplus generation of each array (MWh).  
  * Surplus_Deficits##LS_##MW.pdf: Dual sub-plot figure with (top) local-per-array average proportion of annual load met [%] and (bottom) regional surplus and deficit energy [GWh] for irrigated [blue] and non-irrigated [brown] installations.
  * Derating_Loss_Total_08to18.csv: Dataframe containing total temporally relevant losses from NREL's 2020 Cost Benchmark. Losses include preinverter derate (DC losses * wiring, inverter mismatches, shading, snow loading), inverter efficiency losses (AC losses), and soiling losses. Used in economic analysis. 
  * efficiency_LBNL_08to18_monopoly.csv: Dataframe containing annually reported PV module efficiencies, and linearly extrapoloated monocrystalline and polychristline efficiencies (based on monocrystalline share from LBNL Tracking the Sun Report). Used in economic analysis. 

#### FEWLS_tool/Data/Derived/AdjacentLCC_PostInstall
// Description: Sub-directory for post GEE (AdjacentLCC_PostInstall.txt) dataframes. 
* Files: 
  * AdjacentLCC_PostInstall_#.csv: Dataframes containing installation index, year installed, and CDL landcover type for year pre* and post* installation.


#### FEWLS_tool/Data/Derived/DTW_Raster
// Description: Sub-directory for Depth to Water raster used in the irrigation energy requirement analysis. 
* Files: 
  * Spring2018_DTW.tif: DTW raster derived from Periodic Groundwater Level Measurement Contours (FEWLS_tool/Data/Downloaded/Depth_to_Water/SGMA/).


#### FEWLS_tool/Data/Derived/GEE_ActiveAg_FallowLand_df
// Description: Sub-directory for dataframes exported from Active_AgArea.txt & FallowIdle_area.txt GEE files for values used in article. 
* Files: 
  * cv_active_ag_area.csv: Active agricultural area (km2) in CCV from 2008 to 2018
  * cv_active_ag_area_thru2021.csv: Active agricultural area (km2) in CCV from 2019 to 2021
  * cv_fallow_idle.csv: Fallowed/idle agricultural area (km2) in CCV from 2008 to 2018
  * cv_fallow_idle_thru2021.csv: Fallowed/idle agricultural area (km2) in CCV from 2019 to 2021


#### FEWLS_tool/Data/Derived/Kruit_Stid_Omit
// Description: Sub-directory for combining solar PV arrays present in Kruitwagen et al. (2021) that were omitted in our early detection Stid et al. (2022). 
* Files (none of these files are used in the FEWLS model but are included for FAIR data purposes in the repository): 
  * Crop_Rotation (folder): Folder containing GEE output crop rotation .csv, *kruit_stid_yrPriorcrops.csv*, processed from *PV_Stid_Omit_prep.txt*. Also contains input shape file to GEE script. 
  * GEE_Solar_Analysis (folder): Folder containing final GEE output processed in *PV_Stid_Omit_prep.txt* and merged with dataset form Stid et al. (2022). Also contains input shape file to GEE script. 
  * pvlib_prep_and_results (folder): pvlib electricity generation modeled output for Kruitwagen et al. (2021) arrays omitted from Stid et al. (2022).  
  * Stidetal_OmitKruit_PV_ID.shp: Resulting shapefile of intersection between Kruitwagen et al. (2021) & Stid et al. (2022) to acquire omitted geospatial coordinates.


#### FEWLS_tool/Data/Derived/PV_ID_Complete
// Description: Sub-directory for the spatiotemporal solar PV dataset used in the FEWLS analysis. This is the combination of Stid et al. (2022) and Kruitwagen et al. (2021) solar PV arrays in CCV. 
* Files: 
** CV_CoLocated_Solar_df.shp: Shape file containing geospatial coordinates of direct area's for 935 solar PV arrays installed in CCV between 2008 and 2018. This script contains data from all GEE pre-processing scripts necessary for model run. 


#### FEWLS_tool/Data/Derived/UtilRates
// Description: Sub-directory for the final historical electric utility rate matrices used in the solar PV production economic analysis. Sub-folders contain detailed README's surrounding the generation of the dataframes. 
* Folders:  
  * MercedIrrigDistrict: Folder containing derived .csv files of historical annual electricity rate's (demand, energy, fixed, & other) for the Merced Irrigation District Uitility Provider.
  * ModestoIrrigDistrict: Folder containing derived .csv files of historical annual electricity rate's (demand, energy, fixed, & other) for the Modesto Irrigation District Utility Provider.
  * PGE: Folder containing derived .csv files of historical annual electricity rate's (demand, energy, fixed, & other) for the PG&E Utility Provider. 
  * SCE: Folder containing derived .csv files of historical annual electricity rate's (demand, energy, fixed, & other) for the SCE Utility Provider. 
  * SMUD: Folder containing derived .csv files of historical annual electricity rate's (demand, energy, fixed, & other) for the SMUD Utility Provider. 
  * TurlockIrrigDistrict: Folder containing derived .csv files of historical annual electricity rate's (demand, energy, fixed, & other) for the Turlock Irrigation District Utility Provider. 







## FEWLS_tool/Outputs
// Description: Sub-directory for all model results including figure and dataframe generation. This directory is generated if it doesn't already exist.


### FEWLS_tool/Outputs/Dataframes
// Description: Sub-directory for all model results including figure and dataframe generation. 
// This is empty until the model is run. The github repository contains the default model outputs used in the published article. Naming logic is described below. 
* Files:  
  * Dataset_Information_##yrLS_#MW.csv: For a certain lifespan and capacity cutoff, dataframe containing capacity, total area, number, and number previously irrigated of commercial and utility scale installations in the input dataset. 
  * FEW_Resource_Impacts_##yrLS_#MW.csv: For a certain lifespan and capacity cutoff, dataframe containing cumulative end value (full dataset) food, energy, and water resource effects of agrisolar co-location. Units are descirbed in the dataframe. 
  * Economic_Impacts_##yrLS_#MW.csv: For a certain lifespan and capacity cutoff, dataframe containing cumulative end value (full dataset) economic effects of agrisolar co-location. Units are descirbed in the dataframe. 
  * ResourceModel_CommArrayResults_##yrLS_#MW.csv: For a certain lifespan and capacity cutoff, dataframe containing annually predicted resource results of agrisolar co-location for commercial scale arrays, used in generating final dataframes and figures. 
  * ResourceModel_UtilArrayResults_##yrLS_#MW.csv: For a certain lifespan and capacity cutoff, dataframe containing annually predicted resource results of agrisolar co-location for utility scale arrays, used in generating final dataframes and figures. 
  * EconModel_CommArrayResults_##yrLS_#MW.csv: For a certain lifespan and capacity cutoff, dataframe containing annually predicted economic results of agrisolar co-location for commercial scale arrays, used in generating final dataframes and figures. 
  * EconModel_UtilArrayResults_##yrLS_#MW.csv: For a certain lifespan and capacity cutoff, dataframe containing annually predicted economic results of agrisolar co-location for utility scale arrays, used in generating final dataframes and figures. 


### FEWLS_tool/Outputs/Figures
// Description: Sub-directory for all model results including figure and dataframe generation. 
// This is empty until the model is run. The github repository contains the default model outputs used in the published article. Naming logic is described below. 
* Files: 
  * FEW_results_##yrLS_#MW.pdf: For a certain lifespan and capacity cutoff, figure containing predicted cumulative food, energy, and water resource effects of agrisolar co-location.
  * ECON_results_##yrLS_#MW.pdf: For a certain lifespan and capacity cutoff, figure containing predicted cumulative economic effects of agrisolar co-location.
  * tech_size_dist_##yrLS_#MW.pdf: If the function *getPlot_Array_Tech_Dist(in_solar_df)* is run, for a certain lifespan and capacity cutoff, figure containing commercial and utility fixed and single-axis mounted array distributions through time. 
  * Surplus_Deficits_##yrLS_#MW.pdf: If the function *getFarmElec_Budget(in_solar_df)* is run, for a certain lifespan and capacity cutoff, predictions for electricity generation contributing to annual load vs. surplus electricity generation 
