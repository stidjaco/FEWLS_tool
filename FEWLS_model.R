######################################################################################
#  ____        _              _____ _______        ___     ____    _____           _ #
# / ___|  ___ | | __ _ _ __  |  ___| ____\ \      / / |   / ___|  |_   _|__   ___ | |#
# \___ \ / _ \| |/ _` | '__| | |_  |  _|  \ \ /\ / /| |   \___ \    | |/ _ \ / _ \| |#
#  ___) | (_) | | (_| | |    |  _| | |___  \ V  V / | |___ ___) |   | | (_) | (_) | |#
# |____/ \___/|_|\__,_|_|    |_|   |_____|  \_/\_/  |_____|____/    |_|\___/ \___/|_|#
#                                                                                    #
######################################################################################

## De"script"tion
"
Main model script compiling all model functions and sub-functions to project resource and
economic effects of agrisolar land use change and energy production. This script is called 
in the FEWLS_run.R script, and funcitons are employed in the FEWLS_ResourceModelResults.R 
script, FEWLS_EconModelResults.R script, FEWLS_FarmElecBudgExport.R script, and the 
FEWLS_supplemental.R scritpt. It also calls into memory important dataframes employed in 
the model. Note that the discounted cashflow model (DCF) is only applied by the applyDCF 
function within the getEconomic_Commerical() and getEconomic_Utility() functions, not 
within individual rate functions. 

Functions compiled in this script for the resource and economic models: 

## SETUP FUNCTIONS
- create_directory(): Checks for, and creates directories for Outputs, Figures, and Dataframes.
- getContinous_installPeriodVector(): A vestigial function from when we did not perform operations on a per array basis. It has some value if larger or more complex operations are to be parallelized rather than run on each array. It is retained for that reason, but no longer employed with value in the model. Returns dataframe with empty rows the extent of array lifespans. Requires df as input with Year (integer) attribute. 
- getTemporal_Projection(): Generates a vector of an input variable dependent on the lifespan variable set in FEWLS_inputs.R. Requires an input vector (numeric, ie energy production in 2018), input name (string), vector year (integer), and whether or not to accumulate post vector year or to repeat (boolean, accumulate is TRUE)
- getCrop_rotation(): Returns an expanded dataframe with new rows for each crop in recent (5-years) crop sequence. Requires dataframe input with a direct area (numeric) and CDL crop names (string) for previous 5 years. Note that when applying crop rotation function, any process not concerning direct area, the contribution of the variable needs to be manually split in external code since direct area is the only pre-split variable. We do this throughout existing functions. 
- getTemporal_Capacity(): Returns dataframe of capacity through time. This is redundant for the same reasons as getContinous_installPeriodVector(), but is retained as the functionality is still needed. 
- getTemporal_Area(): Returns dataframe of cumulative area converted (m2) for an input dataframe of one array. 
- getCropTemporal_Area(): Returns dataframe of cumulative area converted (m2) for an input dataframe of one array for each crop type in the crop rotation. 
- cleanFEWLS(): Removes five within-model-derived dataframes (Data/Derived/) that are generated based on input data. These are essential to the model, but can be regenerated with scripts within the model (that are automatically called if these dataframes are absent). Those dataframes are: irrigation depth (mm_per_yearIrrig.csv), crop revenue (crop_revenue_kg.csv), spatial caloric density (kcal_per_m2.csv), and irrigation energy requirements (irrigEnergyReq_adjusted.csv). These files are included in the default FEWLS folder.

## FEW RESOURCE FUNCTIONS 
- getResource_foodkcal(): Returns dataframe of cumulative forgone calories (kcal - base, best, and worst case scenarios) over the selected lifespan of solar land use conversion. Requires input dataframe of single array (in_df).
- getResource_irrigWater(): Returns dataframe of cumulative forgone irrigation water use (m3 - base, best, and worst case scenarios) over the selected lifespan of solar land use conversion. Requires input dataframe of single array (in_df), a boolean for calculating total cumulative irrigation energy as opposed to irrigation water volume (irrig_energy), a boolean for calculating irrigation energy per m2 for the economic irrig-electricity analysis (irrig_energy_per_m2), and a boolean for whether or not to include O&M water use in total offset water volume (include_oandm).
- getResource_OandMWaterUsed(): Returns dataframe of cumulative O&M water use (m3 - base, best, and worst case scenarios) over the selected lifespan of solar lifespan. Requires input dataframe of single array (in_df), and a boolean for including O&M water use for arrays that offset non-irrigated land (onsider_nonIrrig).
- getResource_energyProduced(): Returns dataframe of cumulative electricity generated (GWh - base, best, and worst case scenarios) over the selected solar lifespan. Requires input dataframe of single array (in_df).
- getResource_irrigEnergySaved(): Returns dataframe of cumulative forgone energy use by removing irrigation (GWh - base, best, and worst case scenarios) over the selected lifespan of solar land use conversion. Requires input dataframe of single array (in_df).

## ECONOMIC FUNCTIONS 
- getEconomic_food(): Returns dataframe of cumulative forgone food revenue (base, best, and worst case scenarios) over the selected lifespan of solar land use conversion. Requires input dataframe of single array (in_df).
- getEconomic_energy(): Returns dataframe of cumulative economic value of the electricity produced under NEM (base, best, and worst case scenarios) over the selected solar lifespan. Requires input dataframe of single array (in_df).
- getEconomic_energy_SingleRate(): Returns dataframe of cumulative economic value of the electricity produced under NEM (base, best, and worst case scenarios) over the selected solar lifespan. This function uses a single dataset average rate ($/GWh), while getEconomic_energy() uses locally-relevant utility-specific rates for each array through time. This function is used if rate data is not readily available for each array through time, or for utility-scale sensitivity analysis (changing capacity_threshold). Requires input dataframe of single array (in_df).
- getEconomic_water(): Returns dataframe of cumulative economic value of the change in water use (base, best, and worst case scenarios) over the selected solar lifespan. Requires input dataframe of single array (in_df), and a boolean for whether or not the economic value of electricity (Econ_20##) exists, or should be estimated using the single rate function (econ_exist).
- getEconomic_OandM(): Returns dataframe of cumulative economic cost of O&M (base, best, and worst case scenarios) over the selected solar lifespan. Requires input dataframe of single array (in_df).
- getEconomic_install(): Returns dataframe of cumulative economic cost of installation for commercial scale arrays (base, best, and worst case scenarios) over the selected solar lifespan. Requires input dataframe of single array (in_df).
- getEconomic_landlease(): Returns dataframe of cumulative economic value of land lease (base, best, and worst case scenarios) over the selected lifespan of solar land use change. Requires input dataframe of single array (in_df).
- getEconomic_operationalCost(): Returns dataframe of cumulative economic value of offset agricultural operational costs (base, best, and worst case scenarios) over the selected lifespan of solar land use change. Requires input dataframe of single array (in_df).
- applyDCF(): Returns dataframe in which the discounted cash flow model is applied to each vector. Requires dataframe to input, and discount rate (base, best, and worst case scenarios) to apply. 

## GET MODEL RESULTS FUNCTIONS: 
- getResource_CommUtil(): Returns dataframe of cumulative resource functions (base, best, and worst case scenarios) over the selected lifespan of solar land use change. Requires input dataframe of single array (in_df).
- getEconomic_Commercial(): Returns dataframe of cumulative economic value of all commercial-scale cash flow (base, best, and worst case scenarios) over the selected lifespan of solar land use change. Requires input dataframe of single array (in_df). DCF and Discount rate are applied here. 
- getEconomic_Utility(): Returns dataframe of cumulative economic value of all utility-scale cash flow (base, best, and worst case scenarios) over the selected lifespan of solar land use change. Requires input dataframe of single array (in_df). DCF and Discount rate are applied here. 

## OTHER RESULT FUNCTIONS
- getResource_annFarmOperationReq(): Returns dataframe of electricity contributed to annual load, surplus generation, and deficit generation, for each array. Dataframe contains three rows for each array. Requires input dataframe of all arrays (in_df).
- getResource_irrigWaterFarmFallow(): Returns same dataframe as getResource_irrigWater(), but using average farm area instead of solar direct area to quantify volume of irrigation water-use offset (m * m2). Used in FEWLS_supplemental.R for proximal fallowed land analysis.
- getData_info(): Returns dataframe of input dataset info including total capacity, direct area, number of installations, and number of previously irrigated arrays, delineated by capacity threshold. Requires input dataframe of all arrays (in_df).
- getPlot_Array_Tech_Dist(): Returns figure of number of arrays, mount technology, and average size (capacity), split by capacity threshold over time. Requires input dataframe of all arrays (in_df).
- getFarmElec_Budget(): Returns figure of solar electricity contribution to annual load through time, split by previously irrigated and non-irrigated array land use. Requires input dataframe of all arrays (in_df).
"

#_______________________________________________________________________________________________________________________________________ SETUP START

# Get start time
start_time <- Sys.time()

# Packages for Processing 
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(ggplot2)
library(ggpubr)
library(miscTools)
library(lubridate)
library(tigris)
# library(plyr) -- plyr needs to be installed for model results, but is manually called due to conflict with group_by and summarize in dplyr

# Set Working Directory
setwd(wd)

# Call input file for dataframe locations, conversion factors, and variables. 
source("FEWLS_inputs.R")

# Function to create new output folders if they don't already exist
create_directory = function(in_dir, sub_dir){
  setwd(in_dir)
  if (file.exists(sub_dir)){
    # Do nothing
  } else {
    dir.create(file.path(in_dir, sub_dir))
  }
  setwd(wd)
}

# Create Folders
create_directory(wd, "Outputs")
create_directory(file.path(wd, "Outputs"), "Figures")
create_directory(file.path(wd, "Outputs"), "Dataframes")

# Set global options
options(warn = -1)
options(show.error.messages = TRUE)
options(tigris_use_cache = TRUE)
options(stringsAsFactors = TRUE)

# Compile necessary input files and df's
Nass_Classifications <- read.csv(NASS_classes_loc) 
crop_yield <- read.csv(food_yield_loc)
food_dollarValue <- read.csv(food_dollarValue_loc)
eff_df <- read.csv(energy_efficiency_df_loc)
irrig_perCrop_state <- read.csv(water_irrig_perCrop_state_loc)
irrig_perCrop_conus <- read.csv(water_irrig_perCrop_conus_loc)
USGS_wu_conus <- read.csv(USGS_wu_loc)
county_annu_precip <- read.csv(conus_county_annu_precip_loc)
county_aglandarea <- read.csv(agland_operation_conus_loc)
NREL_ATB_future <- read.csv(NREL_ATB_future_loc)
install_cost_df <- read.csv(NREL_ATB_historical_loc)

# Inflation adjustments for inputs -- CPI and PPI From https://data.bls.gov/PDQWeb/wp -- All commodities and All Urban Customers
inflation_ratePPI <- read.csv(inflation_ratesPPI_loc)
inflation_ratePPI$infl_rate_mult <- inflation_ratePPI$Annual[which(inflation_ratePPI$Year==USD_infl_adj_year)] / inflation_ratePPI$Annual
inflation_rateCPI <- read.csv(inflation_ratesCPI_loc)
inflation_rateCPI$infl_rate_mult <- inflation_rateCPI$Annual[which(inflation_rateCPI$Year==USD_infl_adj_year)] / inflation_rateCPI$Annual

# Build source files for building input dataframes -- Only need to run on startup of model, and if cleanFEWLS is run. Model also checks if these files exist and runs these functions if they dont
#source("FEWLS_mmPerYrIrrig.R")
#source("FEWLS_cropRevenue.R")
#source("FEWLS_kcalperm2.R")
#source("FEWLS_irrigEnergyReq.R")
#source("FEWLS_cropCost.R")

#_______________________________________________________________________________________________________________________________________ SETUP END

## ----------------------------------------------------------- ##
##                                                             ##
##                 Build resource predictions                  ##
##                                                             ##
## ----------------------------------------------------------- ##

#_______________________________________________________________________________________________________________________________________ RESOURCE PREDICTION START

"getContinous_installPeriodVector is a vestigial function from when we did not perform operations on a per array basis. It has some value if larger or more complex operations are to be 
paralellized rather than run on each array."
# Build a funciton to make vectors continuous throughout install period if they are not already continuous -- NOTE: its essential to remove these rows after processing, for clarity
getContinous_installPeriodVector = function(in_df){
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  # Ensure in_df years are integers
  in_df$Year <- in_df$Year %>% as.character() %>% as.integer()
  # Get existing years in input df
  years_df <- unique(in_df$Year) %>% as.integer() %>% sort(decreasing = FALSE)
  years_model <- c(start_install_year:end_install_year)
  # If vectors differ, add arbitrary rows to dataframe with missing years
  if(setequal(years_df, years_model)==FALSE){
    missing_years <- setdiff(years_model, years_df)
    for(year in missing_years){
      in_df <- add_row(in_df, Year=year)
    }
    # Order by year for cumsum operations
    in_df <- in_df[order(in_df$Year), ]
  } 
  return(in_df)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Build a projection function based on model length -- Note: Input vector must contain data for all start year to end year years
getTemporal_Projection = function(input_vector, input_name, vector_year, accumulate_post_endYear){
  # Save final install year value to add 
  endInstallYear_resource = tail(input_vector, 1)
  input_df <- data.frame(Year = vector_year, resource = input_vector)
  # Save a version for post lifespan removeal
  input_InstallPeriod_save <- input_df
  # Add rows for projection -- Note the "-1" used hereafter are to account for including install year as first year of operation (rather than thereafter)
  input_df[nrow(input_df)+(system_lifespan-1),] <- NA
  input_df$Year <- c(vector_year:(vector_year+system_lifespan-1))
  # Project end year resource value through endyear + system lifespan, leave observed values (+1 includes install year, second +1 selects the following row to change)
  for(i in c(2:nrow(input_df))){
    input_df[i,]$resource <- ifelse(accumulate_post_endYear==TRUE, input_df[i-1,]$resource + endInstallYear_resource, input_df[i-1,]$resource)}
  output_df <- input_df
  names(output_df)[2] <- input_name
  return(output_df)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"Note that when applying crop rotation function, any process not concerning direct area, the contribution of the variable needs to be manually split in external 
code since direct area is the only pre-split variable. We do this thorughout existing functions."
# Crop rotation funtion 
getCrop_rotation <- function(in_df){
  # ------------- #
  # Crop Rotation # -- Requires input of 5 years of previous crop types in the form of 'Crp_n_#' 
  # ------------- #
  # Ensure all crop names are lower case 
  in_df$Crp_n_1 <- in_df$Crp_n_1 %>% tolower()
  in_df$Crp_n_2 <- in_df$Crp_n_2 %>% tolower()
  in_df$Crp_n_3 <- in_df$Crp_n_3 %>% tolower()
  in_df$Crp_n_4 <- in_df$Crp_n_4 %>% tolower()
  in_df$Crp_n_5 <- in_df$Crp_n_5 %>% tolower()
  # Check that CDL accuracy is not resulting in a switch between orchard crops and non orchard crops. 
  # Given that orchards are important to the analysis, particularly to water use, and likely not an intermediate crop (shift to and then away) because of their long time-to-yield (mature), remove them from crop history
  # For orchard to be shifted to and away, they would be present in year prior 2, 3, or 4 of the crop history, but not in year prior 1 and 5. 
  # If they are in present in year prior 1 or 5, that might indicate a shift to OR away (prior to solar installation), followed by solar installation, that could result from water-scarcity concerns 
  # For example, a farm tried orchards, failed, and tried solar instead. OR a farm switched away from orchard, or cleared the land, and then installed solar. We want to keep these, because they are possible. 
  nass = Nass_Classifications
  nass$Crop <- nass$Crop %>% tolower()
  nass$Irrig_Crop <- nass$Irrig_Crop %>% tolower()
  FARS_orchards <- nass$Crop[which(nass$Irrig_Crop == "orchards")]
  # Check that intermdediate crop history is not orchards if the bounding crop history is also not orchards (numerical CDL crop code)
  in_df$sys_C_2 <- ifelse(in_df$Crp_n_2 %in% FARS_orchards & !in_df$Crp_n_1 %in% FARS_orchards & !in_df$Crp_n_5 %in% FARS_orchards, NA, in_df$sys_C_2)
  in_df$sys_C_3 <- ifelse(in_df$Crp_n_3 %in% FARS_orchards & !in_df$Crp_n_1 %in% FARS_orchards & !in_df$Crp_n_5 %in% FARS_orchards, NA, in_df$sys_C_3)
  in_df$sys_C_4 <- ifelse(in_df$Crp_n_4 %in% FARS_orchards & !in_df$Crp_n_1 %in% FARS_orchards & !in_df$Crp_n_5 %in% FARS_orchards, NA, in_df$sys_C_4)
  # Do the same for crop name
  in_df$Crp_n_2 <- ifelse(in_df$Crp_n_2 %in% FARS_orchards & !in_df$Crp_n_1 %in% FARS_orchards & !in_df$Crp_n_5 %in% FARS_orchards, NA, in_df$Crp_n_2)
  in_df$Crp_n_3 <- ifelse(in_df$Crp_n_3 %in% FARS_orchards & !in_df$Crp_n_1 %in% FARS_orchards & !in_df$Crp_n_5 %in% FARS_orchards, NA, in_df$Crp_n_3)
  in_df$Crp_n_4 <- ifelse(in_df$Crp_n_4 %in% FARS_orchards & !in_df$Crp_n_1 %in% FARS_orchards & !in_df$Crp_n_5 %in% FARS_orchards, NA, in_df$Crp_n_4)
  # Run crop rotation
  in_df$crop_rotation <- NA
  in_df$unique_crops <- NA
  # Get unique crop types
  for(i in c(1:nrow(in_df))){
    # Create list of CDL crop types for each array five years prior to installation
    in_df[i, ]$crop_rotation <- list(c(in_df[i, ]$sys_C_1, in_df[i, ]$sys_C_2, in_df[i, ]$sys_C_3, in_df[i, ]$sys_C_4, in_df[i, ]$sys_C_5))
    # Remove NA from non crop types
    in_df[i, ]$crop_rotation <- lapply(in_df[i, ]$crop_rotation, function(x) x[!is.na(x) & !x%in%cdl_non_ag])
    # Get unique crops within five year crop list
    in_df[i, ]$unique_crops <- list(c(unique(unlist(in_df[i, ]$crop_rotation, use.names = FALSE))))
  }
  # Unnest by uniuqe crops
  in_df <- in_df %>% unnest(unique_crops)
  # Double check to ensure no erroneous crop values (may be result of GEE export rounding issue)
  in_df <- in_df[which(in_df$unique_crops %in% nass$Value), ]
  # Crop proportion
  in_df$crop_prop <- 0
  for(i in c(1:nrow(in_df))){
    in_df[i, ]$crop_prop <- sum(in_df[i, ]$unique_crops == unlist(in_df[i, ]$crop_rotation)) / lengths(in_df[i, ]$crop_rotation)
  }
  # New total area by crop
  in_df$dir_a <- in_df$dir_a * in_df$crop_prop
  
  return(in_df)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get temporal capacity
getTemporal_Capacity = function(in_df){
  # Remove geometry 
  in_df$geometry = NULL
  vector_year = in_df$Year
  # Get capacity grouped by year
  capacity_df <- in_df %>% 
    group_by(Year) %>% 
    summarize(Capacity = sum(Capacity)/1000) %>%
    as.data.frame()
  # Get temporal Capacity projection
  capacity_df <- getTemporal_Projection(input_vector=capacity_df$Capacity, vector_year=vector_year, input_name='Capacity', accumulate_post_endYear=FALSE)
  return(capacity_df)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get temporal area
getTemporal_Area = function(in_df){
  # Remove geometry, get crop rotation, and set variables
  in_df$geometry = NULL
  vector_year = in_df$Year
  in_df <- getCrop_rotation(in_df)
  # Match and group
  nass_temp <- Nass_Classifications
  in_df$Year <- in_df$Year %>% as.character() %>% as.integer()
  in_df$Crop <- nass_temp$Crop[match(in_df$unique_crops, nass_temp$Value)] %>% tolower()
  # Group by crop and year and sum area
  in_df <- in_df[which(!is.na(in_df$Crop)), ] 
  # Get capacity grouped by year
  area_df <- in_df %>% 
    group_by(Year) %>% 
    summarize(area = sum(dir_a)) %>%
    as.data.frame()
  area_df$area <- area_df$area %>% cumsum()
  # Get temporal Capacity projection
  area_df <- getTemporal_Projection(input_vector=area_df$area, vector_year=vector_year, input_name='area', accumulate_post_endYear=TRUE)
  return(area_df)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get crop specific temporal area
getCropTemporal_Area = function(in_df){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  # --------------------------------------------------------------------------- #
  # Get temporal area for each crop and calculate total operational cost offset # 
  # --------------------------------------------------------------------------- #
  # Remove geometry 
  in_df$geometry = NULL
  vector_year = in_df$Year
  # Get crop roatation
  in_df <- getCrop_rotation(in_df)
  # Match and group
  nass_temp <- Nass_Classifications
  in_df$Year <- in_df$Year %>% as.character() %>% as.integer()
  in_df$Crop <- nass_temp$Crop[match(in_df$unique_crops, nass_temp$Value)] %>% tolower()
  # Group by crop and year and sum area
  group_df <- in_df[which(!is.na(in_df$Crop)), ] %>%
    group_by(Crop, Year) %>% 
    summarise(dir_a = sum(dir_a, na.rm=TRUE))
  # For every crop, calculate total operational offset
  # Create enmpty dataframe to append to, find unique crops, and for every crop in df, run predictions
  df = group_df[which(group_df$Crop=="Non Existent"), ] %>% as.data.frame()
  unique_crops = unique(group_df$Crop) 
  for(crop in unique_crops){
    # Subset Food for each crop
    crop_temp <- group_df[which(group_df$Crop==crop), ] %>% as.data.frame()
    crop_temp <- getContinous_installPeriodVector(crop_temp)
    # Find total area used for solar in each year and project (This initial cumsum is because land converted to solar in one year continues to be converted the next)
    crop_temp$dir_a = cumsum(coalesce(crop_temp$dir_a, 0)) + crop_temp$dir_a*0 # -- replace missing values with 0 and cumsum
    # Fill in missing years (added with getContinous_installPeriodVector function) with appropriate assumptions (this has been cross referenced with the non crop specific method and matches)
    first_crop_year <- as.integer(crop_temp$Year[which(crop_temp$dir_a==first(na.omit(crop_temp$dir_a)))])
    for(i in 1:nrow(crop_temp)){
      crop_temp[i, ]$dir_a <- ifelse(crop_temp[i, ]$Year<first_crop_year, 0, 
                                     ifelse(is.na(crop_temp[i, ]$dir_a), crop_temp[i-1, ]$dir_a, crop_temp[i, ]$dir_a))
    }
    # Find total area used for solar in each year and project (This initial cumsum is because land converted to solar in one year continues to be converted the next)
    crop_temp$Crop = NULL
    crop_temp <- getTemporal_Projection(input_vector=crop_temp$dir_a, vector_year=vector_year, input_name='dir_a', accumulate_post_endYear=TRUE)
    crop_temp$Crop = crop
    df <- rbind(df, crop_temp)
  }
  return(df)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Food kcal projections
getResource_foodkcal = function(in_df){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  
  # Set Variables for yield through time
  food_yield_base = c( (1-yield_deficit_base) + yield_time_base * c(1:(model_length)))
  food_yield_best = c( (1-yield_deficit_best) + yield_time_best * c(1:(model_length)))
  food_yield_worst= c( (1-yield_deficit_worst) + yield_time_worst * c(1:(model_length)))
  
  # -------------------------- #
  # Generate kcal/m2 dataframe # -- Given regional inputs, and the completion of the FEWLS_kcalperm2.R corrections, check if df exists and generate if not 
  # -------------------------- #
  if(file.exists("Data/Derived/kcal_per_m2.csv")){
    # If file already generated in previous model run, call in file
    kcal_sqm <- read.csv("Data/Derived/kcal_per_m2.csv", stringsAsFactors=FALSE)
  } else {
    # If file has not been generated, fun FEWLS_kcalperm2.R to generate, then call in
    source("FEWLS_kcalperm2.R")
    kcal_sqm <- read.csv("Data/Derived/kcal_per_m2.csv", stringsAsFactors=FALSE)
  }
  
  # ------------------------------------ #
  # Crop Rotation grouped by year & crop # -- Results in df with more than initial number of rows because multiple crops per array
  # ------------------------------------ #
  # Remove geometry, get year, get crop rotation
  in_df$geometry = NULL
  vector_year = in_df$Year
  in_df <- getCrop_rotation(in_df)
  # Match and group
  nass_temp <- Nass_Classifications
  in_df$Year <- in_df$Year %>% as.character() %>% as.integer()
  in_df$Crop <- nass_temp$Crop[match(in_df$unique_crops, nass_temp$Value)] %>% tolower()
  # Group by crop and year and sum area
  group_df <- in_df[which(!is.na(in_df$Crop)), ] %>%
    group_by(Crop, Year) %>% 
    summarise(dir_a = sum(dir_a))
  
  # --------------- #
  # Total Crop Kcal # -- Making adjustments for feed/silage crops which go to beef and dairy, and relative seed oil crops kcal for human consumption has been done in FEWLS_kcalperm2.R
  # --------------- #
  # Get matching kcal to 
  group_df$kcal_sqm <- kcal_sqm$kcal_per_m2[match(group_df$Crop, kcal_sqm$crop)]
  nass_temp$Irrig_Crop <- tolower(nass_temp$Irrig_Crop)
  nass_temp$Crop <- tolower(nass_temp$Crop)
  group_df$FARS_crop <- nass_temp$Irrig_Crop[match(group_df$Crop, nass_temp$Crop)] %>% tolower()
  group_df$kcal <- group_df$kcal_sqm * group_df$dir_a

  # Group by year
  Food <- group_df %>% 
    group_by(Year) %>%
    summarize(kcal = sum(kcal))
  # Find total area used for solar in each year and project (This initial cumsum is because land converted to solar in one year continues to be converted the next)
  Food$kcal = Food$kcal %>% cumsum()
  Food <- getTemporal_Projection(input_vector=Food$kcal, input_name='kcal', vector_year=vector_year, accumulate_post_endYear=TRUE)
    
  # ----------------------------- #
  # Best and Worst Case Scenarios # 
  # ----------------------------- #
  Food$kcal_min <- Food$kcal * food_yield_best
  Food$kcal_max <- Food$kcal * food_yield_worst
  Food$kcal <- Food$kcal * food_yield_base
  return(Food)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Water Use Projections and irrigation energy requirements
getResource_irrigWater = function(in_df, irrig_energy = FALSE, irrig_energy_per_m2 = FALSE, include_oandm = FALSE){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 

  # ------------------------ #
  # Get Water use by Crop df # -- Generate if doesnt exist
  # ------------------------ #
  if(file.exists("Data/Derived/mm_per_yearIrrig.csv")){
    # If file already generated in previous model run, call in file
    irrig <- read.csv("Data/Derived/mm_per_yearIrrig.csv", stringsAsFactors=FALSE)
  } else {
    # If file has not been generated, fun FEWLS_kcalperm2.R to generate, then call in
    source("FEWLS_mmPerYrIrrig.R")
    irrig <- read.csv("Data/Derived/mm_per_yearIrrig.csv", stringsAsFactors=FALSE)
  }
  
  # --------------------------------------------------------------------------------- #
  # Prep temporal gridMET precip and USGS water use dataframe to create annual scalar # 
  # --------------------------------------------------------------------------------- #
  # USGS water use by irrigation datframe -- Water Use Data for States: All years + County + All Counties + Irrigation, Total & Irrigation, Crop
  USGS_wu <- USGS_wu_conus[which(USGS_wu_conus$State_Code==statefp), ]
  USGS_wu <- data.frame(Year = USGS_wu$Year, COUNTYFP = USGS_wu$Area_Code, wu = (USGS_wu$IR.WGWFr+USGS_wu$IR.WSWFr))
  # Make Parameter that adjusts Water use from Mgal/day into Acrefeet/year into m3/year
  USGS_wu$wu <- USGS_wu$wu*mm3_year_Mgal_day
  # Select USGS Years
  USGS_yrs = c(1985, 1990, 1995, 2000, 2005, 2010, 2015)
  county_precip_USGS <- county_annu_precip[which(county_annu_precip$STATEFP==statefp & county_annu_precip$Year %in% USGS_yrs), ]
  # Merge dataframes by county and year, set rows to NA if water use = 0 (no data)
  county_precip_USGS$COUNTYFP <- county_precip_USGS$COUNTYFP %>% as.character() %>% as.numeric()
  USGS_wu$COUNTYFP <- USGS_wu$COUNTYFP %>% as.character() %>% as.numeric()
  county_ppt_wu <- merge(USGS_wu, county_precip_USGS, by=c('Year', 'COUNTYFP'))
  # Convert to proper datatypes (numeric)
  county_ppt_wu$wu <- as.numeric(county_ppt_wu$wu)
  # Group by county and calculate linear model for each county. Remove counties with less than 3 reported wu years (no slope or R2)
  county_ppt_wu <- county_ppt_wu %>%
    group_by(COUNTYFP) %>%
    #summarise(PPTslope = coefficients(lm(wu ~ ppt))[2], intercept = coefficients(lm(wu ~ ppt))[1], R_squared = summary(lm(wu ~ ppt))[8])
    summarise(PPTslope = coefficients(lm(wu ~ ppt + Year))[2], YEARslope = coefficients(lm(wu ~ ppt + Year))[3], 
              intercept = coefficients(lm(wu ~ ppt + Year))[1], R_squared = summary(lm(wu ~ ppt + Year))[8], num = n())
  county_ppt_wu$R_squared <- as.numeric(county_ppt_wu$R_squared)
  county_ppt_wu <- county_ppt_wu[which(!is.na(county_ppt_wu$PPTslope)), ]
  # Select map years and USGS years
  all_yrs = c(1985, 1990, 1995, 2000, 2005, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020)
  ppt_df <- county_annu_precip[which(county_annu_precip$STATEFP==statefp & county_annu_precip$Year %in% all_yrs), ]
  ppt_df$COUNTYFP <- ppt_df$COUNTYFP %>% as.character() %>% as.numeric()
  ppt_df$Year <- as.integer(ppt_df$Year)
  ppt_df$ppt <- as.numeric(ppt_df$ppt)
  # Save ppt_df and df for gridMET precip plot
  ppt_df_year <- ppt_df %>%
    group_by(Year) %>%
    summarize(ppt = mean(ppt, na.rm = TRUE))
  ppt_studyperiod <- mean(ppt_df_year$ppt[which(ppt_df_year$Year >= 2008 & ppt_df_year$Year <= 2018)], na.rm = TRUE)
  # Simplify df name, and ensure countyfp are ordered
  df <- county_ppt_wu[order(county_ppt_wu$COUNTYFP), ]
  
  # ---------------------------------------------------------------------- #
  # Get irrig scalar by WU ~ Precip linear model relative to 2013 and 2018 # 
  # ---------------------------------------------------------------------- #
  # Apply regression to df precip of all years relative to 2013 and 2018 (predicting water use in each year)
  relative_func <- function(Year){
    ppt_year <- ppt_df[which(ppt_df$Year==Year), ]
    ppt_year <- ppt_year[order(ppt_year$COUNTYFP), ]
    regression <- df$intercept + df$PPTslope * ppt_year$ppt[ppt_year$COUNTYFP==df$COUNTYFP] # + df$YEARslope * Year
    return(regression)
  }
  # Water Use Relative to 2013
  rl13 <- relative_func(2013)
  df$rl08to13 <- relative_func(2008) / rl13
  df$rl09to13 <- relative_func(2009) / rl13
  df$rl10to13 <- relative_func(2010) / rl13
  df$rl11to13 <- relative_func(2011) / rl13
  df$rl12to13 <- relative_func(2012) / rl13
  df$rl13to13 <- relative_func(2013) / rl13
  df$rl14to13 <- relative_func(2014) / rl13
  df$rl15to13 <- relative_func(2015) / rl13
  df$rl16to13 <- relative_func(2016) / rl13
  df$rl17to13 <- relative_func(2017) / rl13
  df$rl18to13 <- relative_func(2018) / rl13
  # Water Use Relative to 2018
  rl18 <- relative_func(2018)
  df$rl08to18 <- relative_func(2008) / rl18
  df$rl09to18 <- relative_func(2009) / rl18
  df$rl10to18 <- relative_func(2010) / rl18
  df$rl11to18 <- relative_func(2011) / rl18
  df$rl12to18 <- relative_func(2012) / rl18
  df$rl13to18 <- relative_func(2013) / rl18
  df$rl14to18 <- relative_func(2014) / rl18
  df$rl15to18 <- relative_func(2015) / rl18
  df$rl16to18 <- relative_func(2016) / rl18
  df$rl17to18 <- relative_func(2017) / rl18
  df$rl18to18 <- relative_func(2018) / rl18
  rm(rl13, rl18, relative_func, ppt_df, county_ppt_wu)
  # Define scalar by average for each year
  df_13 <- data.frame(
    COUNTYFP = df$COUNTYFP,
    R2 = df$R_squared,
    scalar_2008 = df$rl08to13, 
    scalar_2009 = df$rl09to13, 
    scalar_2010 = df$rl10to13, 
    scalar_2011 = df$rl11to13, 
    scalar_2012 = df$rl12to13, 
    scalar_2013 = df$rl13to13,
    scalar_2014 = df$rl14to13, 
    scalar_2015 = df$rl15to13,
    scalar_2016 = df$rl16to13, 
    scalar_2017 = df$rl17to13, 
    scalar_2018 = df$rl18to13)
  df_18 <- data.frame(
    COUNTYFP = df$COUNTYFP,
    R2 = df$R_squared,
    scalar_2008 = df$rl08to18, 
    scalar_2009 = df$rl09to18, 
    scalar_2010 = df$rl10to18, 
    scalar_2011 = df$rl11to18, 
    scalar_2012 = df$rl12to18, 
    scalar_2013 = df$rl13to18,
    scalar_2014 = df$rl14to18, 
    scalar_2015 = df$rl15to18,
    scalar_2016 = df$rl16to18, 
    scalar_2017 = df$rl17to18, 
    scalar_2018 = df$rl18to18)
  # Set logical scaling limit -- maximum change in reported irrigation between surveys 
  scal_upp_limit = range(irrig$Value[which(irrig$Year==2013)] / irrig$Value[which(irrig$Year==2018)], na.rm = TRUE)[2] %>% as.numeric()
  scal_low_limit = range(irrig$Value[which(irrig$Year==2018)] / irrig$Value[which(irrig$Year==2013)], na.rm = TRUE)[1] %>% as.numeric()
  
  # ------------------------------------------------ #
  # Generate county irrigation energy reqs dataframe # -- Only for CV at the moment 
  # ------------------------------------------------ #
  if(file.exists("Data/Derived/irrigEnergyReq_adjusted.csv")){
    # If file already generated in previous model run, call in file
    irrig_energy_vars <- read.csv("Data/Derived/irrigEnergyReq_adjusted.csv", stringsAsFactors=FALSE)
  } else {
    # If file has not been generated, fun FEWLS_kcalperm2.R to generate, then call in
    source("FEWLS_irrigEnergyReq.R")
    irrig_energy_vars <- read.csv("Data/Derived/irrigEnergyReq_adjusted.csv", stringsAsFactors=FALSE)
  }
  
  # -------------------------------------------- #
  # Get per array irrigation energy requirements #  
  # -------------------------------------------- # 
  # Pull in energy use Dataframe which has been adjusted with improved depth to water from Ben Mcarthy's thesis
  Energy_Use <- irrig_energy_vars[which(irrig_energy_vars$STATE_FIPS==statefp), ]
  # Get county specific GWh / m3
  Energy_Use <- data.frame(COUNTYFP = Energy_Use$FIPS, GWh_m3 = (Energy_Use$Total_Energy/joules_GWh) / (Energy_Use$GW_Total + Energy_Use$SW_Total))
  # Get County data to match to county FIPS
  counties <- counties(state, TRUE) %>% st_as_sf() %>% st_set_crs(st_crs(in_df))
  intersection <- st_intersection(st_centroid(in_df), counties)
  in_df$County <- intersection$NAME[match(in_df$Index, intersection$Index)]
  in_df$COUNTYFP <- counties$COUNTYFP[match(in_df$County, counties$NAME)]
  in_df$COUNTYFP <- gsub("(?<![0-9])0+", "", in_df$COUNTYFP, perl = TRUE) # Extracts non-zero led countyfp numerically from a string
  counties$geometry = NULL
  # Match county names to energy use df
  Energy_Use$COUNTYFP <- str_remove(as.character(Energy_Use$COUNTYFP), paste("^", statefp, "+", sep="")) 
  counties$COUNTYFP <- as.character(counties$COUNTYFP)
  Energy_Use$County <- counties$NAME[match(Energy_Use$COUNTYFP, counties$COUNTYFP)] 
  rm(counties)
  # Match GWh per m3 by counties in crop_df, and multiply by irrig (m3)
  in_df$irrig_energy <- Energy_Use$GWh_m3[match(in_df$County, Energy_Use$County)] 

  # Remove geometry 
  in_df$geometry = NULL
  
  # ------------------------------------ #
  # Crop Rotation grouped by year & crop # 
  # ------------------------------------ #
  vector_year = in_df$Year
  in_df <- getCrop_rotation(in_df)
  # Match and group
  nass_temp <- Nass_Classifications
  in_df$Year <- in_df$Year %>% as.character() %>% as.integer()
  in_df$Crop <- nass_temp$Crop[match(in_df$unique_crops, nass_temp$Value)] %>% tolower()
  # Get non missing croptypes
  crop_water <- in_df[which(!is.na(in_df$Crop)), ]
  
  # ----------------------------------------------------------------------------------------------------------- #
  # Assign Irrig depths based on crop and year, scale using scalar, and take average of 2013 & 2018 predictions # 
  # ----------------------------------------------------------------------------------------------------------- #
  # Call in NASS to get irrig crop types
  nass_temp$Crop <- nass_temp$Crop %>% tolower()
  crop_water$Crop <- nass_temp$Irrig_Crop[match(crop_water$Crop, nass_temp$Crop)] %>% tolower
  
  # Get appropriate scalar values and multiply by area (end result is m3)
  crop_water$irrig_rl13 = 0
  crop_water$irrig_rl18 = 0
  crop_water$Year <- crop_water$Year %>% as.character() %>% as.integer()
  for(i in c(1:nrow(crop_water))){
    # Get necessary variables
    crop <- crop_water[i, ]$Crop
    Year <- crop_water[i, ]$Year
    countyfp <- crop_water[i, ]$COUNTYFP
    # Get scalars
    scalar_2013 <- df_13[[paste("scalar_", Year, sep="")]][which(df_13$COUNTYFP==countyfp)]
    scalar_2018 <- df_18[[paste("scalar_", Year, sep="")]][which(df_18$COUNTYFP==countyfp)]
    # Check if scalars are greater than limits
    scalar_2013 <- ifelse(scalar_2013 >= scal_upp_limit, scal_upp_limit, ifelse(scalar_2013 <= scal_low_limit, scal_low_limit, scalar_2013))
    scalar_2018 <- ifelse(scalar_2013 >= scal_upp_limit, scal_upp_limit, ifelse(scalar_2018 <= scal_low_limit, scal_low_limit, scalar_2018))
    # Get irrig values
    irrig_rl13 <- irrig$Value[which(irrig$Crop==crop & irrig$Year==2013)] * scalar_2013 / mm_m
    irrig_rl18 <- irrig$Value[which(irrig$Crop==crop & irrig$Year==2018)] * scalar_2018 / mm_m
    # Assign scalar and irrig value to crop
    crop_water[i, ]$irrig_rl13 <- crop_water[i, ]$dir_a * irrig_rl13 # m3
    crop_water[i, ]$irrig_rl18 <- crop_water[i, ]$dir_a * irrig_rl18 # m3
  }
  # Get the average prediction
  crop_water$irrig_fg <- (crop_water$irrig_rl13 + crop_water$irrig_rl18) / 2
  
  # ----------------------------- #
  # Best and Worst Case Scenarios # 
  # ----------------------------- #
  crop_water$irrig_rl13 = 0
  crop_water$irrig_rl18 = 0
  crop_water$irrig_fg_dry = 0
  wet_dry_yrs <- ppt_df_year[which(ppt_df_year$Year >= start_install_year & ppt_df_year$Year <= end_install_year), ]
  wet_year = wet_dry_yrs$Year[which(wet_dry_yrs$ppt==max(wet_dry_yrs$ppt))] %>% as.numeric()
  dry_year = wet_dry_yrs$Year[which(wet_dry_yrs$ppt==min(wet_dry_yrs$ppt))] %>% as.numeric()
  for(i in c(1:nrow(crop_water))){
    # Get necessary variables
    crop <- crop_water[i, ]$Crop
    countyfp <- crop_water[i, ]$COUNTYFP
    # Get scalars
    scalar_2013 <- df_13[[paste("scalar_", dry_year, sep="")]][which(df_13$COUNTYFP==countyfp)]
    scalar_2018 <- df_18[[paste("scalar_", dry_year, sep="")]][which(df_18$COUNTYFP==countyfp)]
    # Check if scalars are greater than limits
    scalar_2013 <- ifelse(scalar_2013 >= scal_upp_limit, scal_upp_limit, ifelse(scalar_2013 <= scal_low_limit, scal_low_limit, scalar_2013))
    scalar_2018 <- ifelse(scalar_2013 >= scal_upp_limit, scal_upp_limit, ifelse(scalar_2018 <= scal_low_limit, scal_low_limit, scalar_2018))
    # Get irrig values
    irrig_rl13 <- irrig$Value[which(irrig$Crop==crop & irrig$Year==2013)] * scalar_2013 / mm_m
    irrig_rl18 <- irrig$Value[which(irrig$Crop==crop & irrig$Year==2018)] * scalar_2018 / mm_m
    # Assign scalar and irrig value to crop
    crop_water[i, ]$irrig_rl13 <- crop_water[i, ]$dir_a * irrig_rl13 # m3
    crop_water[i, ]$irrig_rl18 <- crop_water[i, ]$dir_a * irrig_rl18 # m3
    # Rather than average prediction, assign largest potential offset (best case) for dry year
    crop_water[i, ]$irrig_fg_dry <- max(crop_water[i, ]$irrig_rl13, crop_water[i, ]$irrig_rl18, na.rm = TRUE)
  }
  # Get the average
  #crop_water$irrig_fg_dry <- (crop_water$irrig_rl13 + crop_water$irrig_rl18) / 2
  crop_water$irrig_rl13 = 0
  crop_water$irrig_rl18 = 0
  crop_water$irrig_fg_wet = 0
  # Scenario for wet projection (2017)
  for(i in c(1:nrow(crop_water))){
    # Get necessary variables
    crop <- crop_water[i, ]$Crop
    countyfp <- crop_water[i, ]$COUNTYFP
    # Get scalars
    scalar_2013 <- df_13[[paste("scalar_", wet_year, sep="")]][which(df_13$COUNTYFP==countyfp)]
    scalar_2018 <- df_18[[paste("scalar_", wet_year, sep="")]][which(df_18$COUNTYFP==countyfp)]
    # Check if scalars are greater than limits
    scalar_2013 <- ifelse(scalar_2013 >= scal_upp_limit, scal_upp_limit, ifelse(scalar_2013 <= scal_low_limit, scal_low_limit, scalar_2013))
    scalar_2018 <- ifelse(scalar_2013 >= scal_upp_limit, scal_upp_limit, ifelse(scalar_2018 <= scal_low_limit, scal_low_limit, scalar_2018))
    # Get irrig values
    irrig_rl13 <- irrig$Value[which(irrig$Crop==crop & irrig$Year==2013)] * scalar_2013 / mm_m
    irrig_rl18 <- irrig$Value[which(irrig$Crop==crop & irrig$Year==2018)] * scalar_2018 / mm_m
    # Assign scalar and irrig value to crop
    crop_water[i, ]$irrig_rl13 <- crop_water[i, ]$dir_a * irrig_rl13 # m3
    crop_water[i, ]$irrig_rl18 <- crop_water[i, ]$dir_a * irrig_rl18 # m3
    # Rather than average prediction, assign smallest potential offset (best case) for wet year
    crop_water[i, ]$irrig_fg_wet <- min(crop_water[i, ]$irrig_rl13, crop_water[i, ]$irrig_rl18, na.rm = TRUE)
  }
  # Get the average
  #crop_water$irrig_fg_wet <- (crop_water$irrig_rl13 + crop_water$irrig_rl18) / 2
  crop_water$irrig_rl13 = NULL
  crop_water$irrig_rl18 = NULL

  # ------------------------------------------------- #
  # Decide Irrigation water use saved or energy saved # 
  # ------------------------------------------------- #
  if(irrig_energy==TRUE){
    # ----------------- # -- 
    # Get O&M Water Use # -- There is also an independent way of looking at this in a below function, but to include a total change in water use, 
    # ----------------- # -- O&M needs to be account for here as well (only in regard to energy--if we are only interested in irrig water offset, discount this)
    # Operation and Maintance water usage
    if(include_oandm==TRUE){
      # Get county specific irrig energy requirements and group -- Regardless of O&M inclusion -- Crop Prop still used because capacity is split between every crop in rotation
      crop_water$irrig_fg <- crop_water$irrig_fg - (crop_water$Capacity * crop_water$crop_prop * baseOandM_wu)
      crop_water$irrig_fg_dry <- crop_water$irrig_fg_dry - (crop_water$Capacity * crop_water$crop_prop  * lowOandM_wu)
      crop_water$irrig_fg_wet <- crop_water$irrig_fg_wet - (crop_water$Capacity * crop_water$crop_prop  * highOandM_wu)
    }
    
    # Get county specific irrig energy requirements and group -- Regardless of O&M inclusion
    crop_water$irrig_fg <- crop_water$irrig_fg * crop_water$irrig_energy
    crop_water$irrig_fg_dry <- crop_water$irrig_fg_dry * crop_water$irrig_energy
    crop_water$irrig_fg_wet <- crop_water$irrig_fg_wet * crop_water$irrig_energy
  }
  
  # Save dataframe for export to other functions
  crop_water_save = crop_water %>% 
    group_by(Index) %>%
    summarize(irrig_fg = sum(irrig_fg), irrig_fg_dry = sum(irrig_fg_dry), irrig_fg_wet = sum(irrig_fg_wet), area = sum(dir_a))

  # Group by year and crop save in different df to save crop specific info
  Water <- crop_water %>% 
    group_by(Year) %>%
    summarize(irrig_fg = sum(irrig_fg), irrig_fg_dry = sum(irrig_fg_dry), irrig_fg_wet = sum(irrig_fg_wet))
  # Get cumulative sum of water use, get temproal projection, and save to df
  Water$irrig_fg = Water$irrig_fg %>% cumsum()
  Water$irrig_fg_dry = Water$irrig_fg_dry %>% cumsum()
  Water$irrig_fg_wet = Water$irrig_fg_wet %>% cumsum()
  Water_base <- getTemporal_Projection(input_vector=Water$irrig_fg, input_name='irrig_fg', vector_year=vector_year, accumulate_post_endYear=TRUE)
  Water_dry <- getTemporal_Projection(input_vector=Water$irrig_fg_dry, input_name='irrig_fg_dry', vector_year=vector_year, accumulate_post_endYear=TRUE)
  Water_wet <- getTemporal_Projection(input_vector=Water$irrig_fg_wet, input_name='irrig_fg_wet', vector_year=vector_year, accumulate_post_endYear=TRUE)
  Water_base$irrig_fg_dry <- Water_dry$irrig_fg_dry
  Water_base$irrig_fg_wet <- Water_wet$irrig_fg_wet
  Water <- Water_base
  rm(Water_dry, Water_wet, Water_base, Water_temp)
  
  # If get irrig_energy_per_m2 for every array is true (for electricity economic analysis)
  if(irrig_energy_per_m2==TRUE){
    Water = crop_water_save %>% as.data.frame()
  } else {
    # If non-irrigated and non energy dependent, return empty df # This computational saving method had to be removed since the basis for annual on-farm load estiamtion requires irrigation requirements even for unirrigated land
    if(first(in_df$irrig)==0){Water = data.frame(Year = c(in_df$Year:(in_df$Year+system_lifespan-1)), irrig_fg = rep(0, system_lifespan), irrig_fg_dry = rep(0, system_lifespan), irrig_fg_wet = rep(0, system_lifespan))}
  }

 # Return df
 return(Water)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Operation and Maintance Water Use Prediction
getResource_OandMWaterUsed = function(in_df, consider_nonIrrig = FALSE){
  # Get temporal capacity #
  OandM_wu <- getTemporal_Capacity(in_df)
  # Operation and Maintance water usage
  OandM_wu$oandm_base <- c(baseOandM_wu * OandM_wu$Capacity*1000) %>% cumsum()
  OandM_wu$oandm_min <- c(lowOandM_wu * OandM_wu$Capacity*1000)  %>% cumsum()
  OandM_wu$oandm_max <- c(highOandM_wu * OandM_wu$Capacity*1000)  %>% cumsum()
  OandM_wu$Capacity = NULL
  # If we are not considering O&M water use for non-irrigated arrays, results in df of zeros 
  if(consider_nonIrrig == FALSE){
    if(in_df$irrig==0){
      OandM_wu$oandm_base <- rep(0, nrow(OandM_wu))
      OandM_wu$oandm_min <- rep(0, nrow(OandM_wu))
      OandM_wu$oandm_max <- rep(0, nrow(OandM_wu))
    }
  }
  # Return Dataframe
  return(OandM_wu)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Solar Energy Production Projections
getResource_energyProduced = function(in_df){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  
  # For degradation, need to know mount technology - No longer used because Jordan et al., 2022 reports no stat difference between mounts
  #mount = in_df$Class
  
  # ------------------------------------ #
  # Crop Rotation grouped by year & crop # 
  # ------------------------------------ #
  # Remove geometry 
  in_df$geometry = NULL
  # Select only irrigated crops
  in_df <- getCrop_rotation(in_df)
  # New total area by crop -- because cropRotation fucntion only splits direct area by crop_prop
  in_df[, grepl("Gen", colnames(in_df))] <- in_df[, grepl("Gen", colnames(in_df))] * in_df$crop_prop
  # Match and group
  nass_temp <- Nass_Classifications
  in_df$Year <- in_df$Year %>% as.character() %>% as.integer()
  in_df$Crop <- nass_temp$Crop[match(in_df$unique_crops, nass_temp$Value)] %>% tolower()
  in_df$crop_prop = NULL
  in_df$unique_crops = NULL
  in_df$crop_rotation = NULL
  in_df <- in_df %>% dplyr::select(matches("Year|Crop|Gen"))
  # Group by crop and year and sum dir_a
  crop_energy <- in_df[which(!is.na(in_df$Crop)), ] %>%
    group_by(Crop, Year) %>% 
    summarise_all(sum, na.rm=TRUE)
  
  # ----------------------------------------------- #
  # Get linear model of mono vs multi for scenarios #  
  # ----------------------------------------------- #
  # Linear model of efficiency with mono share
  lm <- lm(eff_df$Efficiency ~ eff_df$Mono_share + eff_df$Year)
  r2 <- summary(lm(eff_df$Efficiency~eff_df$Mono_share  + eff_df$Year))$adj.r.squared
  predict <- predict(lm)
  # Set mono share equal to 0 and 1
  mono_share <- 0
  mono_none <- lm$coefficients[1] + lm$coefficients[2]*(mono_share) + lm$coefficients[3]*(c(start_install_year:end_install_year))
  mono_share <- 1
  mono_all <- lm$coefficients[1] + lm$coefficients[2]*(mono_share) + lm$coefficients[3]*(c(start_install_year:end_install_year))
  # Group into df of mono and multi efficiencies
  eff_df_moml <- data.frame(Year = c(2008:2018), efficiency_mono = mono_all, efficiency_multi = mono_none, 
                            efficiency_reported = eff_df$Efficiency[which(eff_df$Year>=start_install_year & eff_df$Year<=end_install_year)])
  
  # -------------------------------------------------------------- #
  # For every crop, perform electricity projection, and save to df # -- This can be done because in the above crop rotation, all generations are split by crop_prop
  # -------------------------------------------------------------- #
  # Create enmpty dataframe to append to, find unique crops, and for every crop in df, run predictions
  df <- data.frame(Year=integer(), Crop=character(), energy_base=numeric(), energy_min=numeric(), energy_max=numeric())
  namesdf <- names(df)
  unique_crops = unique(crop_energy$Crop) 
  for(crop in unique_crops){ 
    # Subset Water for each crop
    energy_temp <- crop_energy[which(crop_energy$Crop==crop), ] %>% as.data.frame()
    #energy_temp <- getContinous_installPeriodVector(energy_temp)
    # Mono entirety function
    mono_function <- function(YOD){
      # Save install year and grab all generations
      Year <- in_df$Year
      in_df <- dplyr::select(in_df, contains("Gen"))
      # New efficiency / Old efficiency = generation multiplier
      eff_converstion <- eff_df_moml$efficiency_mono[eff_df_moml$Year==YOD] / eff_df$Efficiency[eff_df$Year==YOD]
      for(j in c(1:install_period)){
        in_df <- in_df * eff_converstion}
      # Return Year
      in_df$Year <- Year
      return(in_df)
    }
    # Multi entirety function
    multi_function <- function(YOD){
      # Save install year and grab all generations
      Year <- in_df$Year
      in_df <- dplyr::select(in_df, contains("Gen"))
      in_df$Year <- Year
      # New efficiency / Old efficiency = generation multiplier
      eff_converstion <- eff_df_moml$efficiency_multi[eff_df_moml$Year==YOD] / eff_df$Efficiency[eff_df$Year==YOD]
      for(j in c(1:install_period)){
        in_df <- in_df * eff_converstion}
      # Return Year
      in_df$Year <- Year
      return(in_df)
    }
    # Mono vs multi limits
    in_df = energy_temp
    yr_list <- c(end_install_year:start_install_year)
    mono_run <- t(lapply(yr_list, FUN = mono_function))
    multi_run <- t(lapply(yr_list, FUN = multi_function))
    Energy_mono <- do.call("rbind", mono_run)
    Energy_multi <- do.call("rbind", multi_run)
    Energy_base <- in_df %>% dplyr::select(matches("Gen|Year"))

    # ------------------------------------------------------- #
    # Project electricity generation forward with degradation # -- This is a special projection, so we don't use the getTemporal_Projection function
    # ------------------------------------------------------- #
    # Account for fully multi or fully monocrystalline
    energy_project <- function(df){
      # Convert to GWh
      year = df$Year %>% first() %>% as.numeric()
      df[, grepl("Gen", names(df))] <- df[, grepl("Gen", names(df))] / 1000
      # Create generation dataframe with gen and year of inst contribution
      Energy <- df %>%  
        group_by(Year) %>% 
        summarize_all(sum, na.rm=TRUE)
      Energy$Year = NULL
      Energy <- Energy %>% t() %>% as.data.frame()
      # Get first row of generation & First year-row to project from
      gen_yrs_start <- which(rownames(Energy)==paste("Gen_", year, sep="")) %>% as.numeric()
      project_start_yearRow = nrow(Energy)+1
      # Add rows to 2042 -- 25 years post 2018 installation, and year column
      Energy[c((nrow(Energy)+1):(gen_yrs_start+system_lifespan-1)),] <- 0 # Saves all modeled generation and generates empty df for system lifespan
      # Project 2018 generation (GWh) through 2042, leave modeled values
      # New degradation rates based on median rates delineated by mount tech from Jordan et al., 2022, slightly more conservative and updated (https://doi.org/10.1002/pip.3566)
      for(i in project_start_yearRow:nrow(Energy)){
        Energy[i,] <- Energy[gen_yrs_start, ] - ( (i - gen_yrs_start) * efficiency_degradation_rate * Energy[gen_yrs_start, ]) # -0.75%/yr
      }
      # Limit to first year of generation and Return df
      Energy <- Energy[c(gen_yrs_start:nrow(Energy)), ] %>% as.data.frame()
      return(Energy)
    }
    # Get projected energy
    Energy_mono <- energy_project(Energy_mono) %>% cumsum()
    Energy_base <- energy_project(Energy_base) %>% cumsum()
    Energy_multi <- energy_project(Energy_multi) %>% cumsum()
    # Combine into single df
    Energy <- data.frame(Year = c(tail(yr_list, 1):(tail(yr_list, 1)+system_lifespan-1)), base = Energy_base, max = Energy_mono, min = Energy_multi)
    names(Energy) = c("Year", "energy_base", "energy_max", "energy_min")
    Energy$Crop <- crop
    df <- rbind(df, Energy) 
    }

  # Get data from df 
  df$Crop = NULL
  Energy_final <- df %>% group_by(Year) %>% summarise_all(sum, na.rm=TRUE)
  return(as.data.frame(Energy_final))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Irrigation energy water use savings (budget of )
getResource_irrigEnergySaved = function(in_df){ # Includes O&M water use in total water budget for energy and economic implications
  # -------------------------------------------------- #
  # Run irrigation water function with energy variable #  
  # -------------------------------------------------- #
  irrig_energy_df <- getResource_irrigWater(in_df, irrig_energy = TRUE, irrig_energy_per_m2 = FALSE, include_oandm = TRUE) # Currently, O&M water use is still included in utility-scale estimations, because from an electricity offset point of view (even if not from the farmers point of view), this is still required electricity and water
  return(irrig_energy_df)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Crop specific forgone food prices
getEconomic_food = function(in_df){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  
  # Set Variables for yield through time
  food_yield_base = c( (1-yield_deficit_base) + yield_time_base * c(1:(model_length)))
  food_yield_best = c( (1-yield_deficit_best) + yield_time_best * c(1:(model_length)))
  food_yield_worst= c( (1-yield_deficit_worst) + yield_time_worst * c(1:(model_length)))
  
  # --------------------------- #
  # Get USDA Survey Crop Prices # 
  # --------------------------- #
  if(file.exists("Data/Derived/crop_revenue_kg.csv")){
    # If file already generated in previous model run, call in file
    cropdollar <- read.csv("Data/Derived/crop_revenue_kg.csv", stringsAsFactors=FALSE)
  } else {
    # If file has not been generated, fun FEWLS_cropRevenue.R to generate, then call in
    source("FEWLS_cropRevenue.R")
    cropdollar <- read.csv("Data/Derived/crop_revenue_kg.csv", stringsAsFactors=FALSE)
  }
  
  # ------------------------------------------ #
  # Generate kcal/m2 dataframe with kg/m2 data # -- Redundant as this is done above, but checks again to ensure file exists and updated
  # ------------------------------------------ #
  if(file.exists("Data/Derived/kcal_per_m2.csv")){
    # If file already generated in previous model run, call in file
    kcal_sqm <- read.csv("Data/Derived/kcal_per_m2.csv", stringsAsFactors=FALSE)
  } else {
    # If file has not been generated, fun FEWLS_kcalperm2.R to generate, then call in
    source("FEWLS_kcalperm2.R")
    kcal_sqm <- read.csv("Data/Derived/kcal_per_m2.csv", stringsAsFactors=FALSE)
  }
  
  # Group and average by food prices post install year. This is a simplification but technically provides the same total price * area for the period of ag-census data
  df <- data.frame(dollar_kg = numeric(), Crop = character())
  for(crop in unique(cropdollar$Crop)){
    # Get crop
    cropdollar_crop <- cropdollar[which(cropdollar$Crop==crop), ]
    # Get crop prices after installation year
    cropdollar_temp <- cropdollar_crop[which(cropdollar_crop$Year>=start_install_year), ] %>% group_by(Crop) %>% summarize(dollar_kg = mean(dollar_kg, na.rm=TRUE))
    # If price data is missing for years after install, use historical average
    cropdollar_missPrice <- cropdollar_crop %>% group_by(Crop) %>% summarize(dollar_kg = mean(dollar_kg, na.rm=TRUE))
    
    # Decide which to use, if price is reported post install
    dollar_kg <- ifelse(nrow(cropdollar_temp)==0, cropdollar_missPrice$dollar_kg, cropdollar_temp$dollar_kg)
    # Return to new df 
    df_temp <- data.frame(dollar_kg = dollar_kg, Crop = crop)
    df <- rbind(df, df_temp)
  }
  
  # Return cropdollar
  cropdollar = df
  
  # -------------------------------------- #
  # Crop Rotation and merging of and kg/m2 # 
  # -------------------------------------- #
  # Remove geometry and run crop rotation
  in_df$geometry = NULL
  vector_year = in_df$Year
  in_df <- getCrop_rotation(in_df)
  # Match and group
  nass_temp <- Nass_Classifications
  in_df$Year <- in_df$Year %>% as.character() %>% as.integer()
  in_df$Crop <- nass_temp$Crop[match(in_df$unique_crops, nass_temp$Value)] %>% tolower()
  nass_temp$Irrig_Crop <- tolower(nass_temp$Irrig_Crop)
  nass_temp$Crop <- tolower(nass_temp$Crop)
  in_df$FARS_crop <- nass_temp$Irrig_Crop[match(in_df$Crop, nass_temp$Crop)] %>% tolower()
  # Match yield of crops to tot_a df -- not temporal so is consistent
  in_df$kg_per_m2 = 0
  in_df$dollar_kg = 0
  for(i in c(1:nrow(kcal_sqm))){
    # Get forgone crop type and year from survey data
    crop <- kcal_sqm[i, ]$crop
    kg_per_m2 <- kcal_sqm[i, ]$kg_per_m2 %>% as.numeric()
    # Get matching crop type from dataset
    in_df[which(grepl(crop, in_df$Crop)), ]$kg_per_m2 <- kg_per_m2
  }
  
  # -------------------------------------------- #
  # Get annual $/kg and apply to array over time # -- This is done for each array and each arrays area is projected forward and the price for subsequent years is grabbed for temporal precision
  # -------------------------------------------- #
  # For every array, get forgone crop, project forgone crop area forward, grab relative price, and regroup. Crop and year specific prices
  in_df$dollar <- 0
  in_df$dollar_m2 <- 0
  for(i in c(1:nrow(in_df))){
    # Get continous vector for array
    in_df_temp <- getContinous_installPeriodVector(in_df[i, ])
    # Fill in missing years (added with getContinous_installPeriodVector function) with appropriate assumptions for dir_a and yield (kg_per_m2)
    first_crop_year <- in_df[i, ]$Year %>% as.numeric()
    # Redundant, from when we processed the dataset as a whole rather than on a per array basis (getContinousInstallVector solution)
    for(ii in 1:nrow(in_df_temp)){
      in_df_temp[ii, ]$dir_a <- ifelse(in_df_temp[ii, ]$Year<first_crop_year, 0, ifelse(is.na(in_df_temp[ii, ]$dir_a), in_df_temp[ii-1, ]$dir_a, in_df_temp[ii, ]$dir_a))
      in_df_temp[ii, ]$kg_per_m2 <- ifelse(in_df_temp[ii, ]$Year<first_crop_year, 0, ifelse(is.na(in_df_temp[ii, ]$kg_per_m2), in_df_temp[ii-1, ]$kg_per_m2, in_df_temp[ii, ]$kg_per_m2))
    }
    # Get cumulative area for the array through install period
    #in_df_temp$dir_a <- cumsum(coalesce(in_df_temp$dir_a, 0)) + in_df_temp$dir_a*0 -- Prices not reflected "cumulatively" initial error in concept 
    in_df_temp$Crop <- c(rep(in_df[i, ]$Crop, nrow(in_df_temp)))
    
    # Subset crop dollar to reduce computation -- WITH CDL CROP
    crop_string <- scan(text = as.character(in_df[i, ]$Crop), what = " ")
    string_ext = as.character("")
    for(crop in crop_string){
      string_ext = paste(string_ext, "|", crop, sep = "")  }
    string_ext <- substring(string_ext, 2) %>% tolower() # removes | at start of string
    cropdollar_temp <- cropdollar[which(grepl(string_ext, cropdollar$Crop)), ]
    if(nrow(cropdollar_temp)>0){
      # Match dollar value of crops to tot_a df, multiply to get ($/m2), and multiply by tot_a in m2 to get $
      for(ii in c(1:nrow(cropdollar_temp))){
        # Get forgone crop type and year from survey data
        crop <- cropdollar_temp[ii, ]$Crop
        #year <- cropdollar_temp[ii, ]$Year
        dollar <- cropdollar_temp[ii, ]$dollar_kg
        # Get matching crop type from dataset, and fill in cost per kg for crop every year from start year to end year
        in_df_temp[which(grepl(crop, in_df_temp$Crop)), ]$dollar_kg <- dollar
        #in_df_temp[which(grepl(crop, in_df_temp$Crop) & in_df_temp$Year==year), ]$dollar_kg <- dollar
      }
    }
    # Subset crop dollar to reduce computation -- WITH FARS CROP
    crop_string <- scan(text = as.character(in_df[i, ]$FARS_crop), what = " ", quiet = TRUE)
    string_ext = as.character("")
    for(crop in crop_string){
      string_ext = paste(string_ext, "|", crop, sep = "")  }
    string_ext <- substring(string_ext, 2) %>% tolower() # removes | at start of string
    cropdollar_tempFARS <- cropdollar[which(grepl(string_ext, cropdollar$Crop)), ]
    if(nrow(cropdollar_tempFARS)>0){
      # Match dollar value of crops to tot_a df, multiply to get ($/m2), and multiply by tot_a in m2 to get $
      for(ii in c(1:nrow(cropdollar_tempFARS))){
        # Get forgone crop type and year from survey data
        crop <- cropdollar_tempFARS[ii, ]$Crop
        #year <- cropdollar_tempFARS[ii, ]$Year
        dollar <- cropdollar_tempFARS[ii, ]$dollar_kg
        # Get matching crop type from dataset, and fill in cost per kg for crop every year from start year to end year
        in_df_temp[which(grepl(crop, in_df_temp$FARS_crop)), ]$dollar_kg <- dollar
        #in_df_temp[which(grepl(crop, in_df_temp$FARS_crop) & in_df_temp$Year==year), ]$dollar_kg <- dollar
      }
    } else {
      # Manually assume equivalents for still missing crops after for loop
    }
    # Calculate total dollar amount for this array and assign to appropriate row in original df
    in_df[i, ]$dollar_m2 <- mean(in_df_temp$dollar_kg, na.rm=TRUE) * mean(in_df_temp$kg_per_m2[which(in_df_temp$kg_per_m2>0)], na.rm=TRUE) # this simply saves the value for manual additions of NA price received crops later
    in_df[i, ]$dollar <- sum(in_df_temp$dir_a * in_df[i, ]$dollar_m2, na.rm = TRUE)
  }
  
  # Manually assume equivalents for still missing crops -- These may change depending on if new crops arise -- multiplying the direct area times the number of years through end year of analysis gives cumulative impact of each array to end year
  # Remove from all because of crop price grouping -- & cropdollar$Year==in_df$Year
  in_df$dollar <- ifelse(in_df$Crop=="grassland/pasture", cropdollar$dollar_kg[which(cropdollar$Crop=="hay")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="clover/wildflowers", cropdollar$dollar_kg[which(cropdollar$Crop=="hay")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="vetch", cropdollar$dollar_kg[which(cropdollar$Crop=="hay")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="triticale", cropdollar$dollar_kg[which(cropdollar$Crop=="wheat")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="small grains", cropdollar$dollar_kg[which(cropdollar$Crop=="wheat")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="pomegranates", cropdollar$dollar_kg[which(cropdollar$Crop=="oranges")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="citrus", cropdollar$dollar_kg[which(cropdollar$Crop=="oranges")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="peas", cropdollar$dollar_kg[which(cropdollar$Crop=="broccoli")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="soybeans", cropdollar$dollar_kg[which(cropdollar$Crop=="beans")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="watermelons", cropdollar$dollar_kg[which(cropdollar$Crop=="melons")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="cantaloupes", cropdollar$dollar_kg[which(cropdollar$Crop=="melons")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="switchgrass", cropdollar$dollar_kg[which(cropdollar$Crop=="hay")] * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="berry totals", mean(cropdollar$dollar_kg[which(cropdollar$Crop %in% c("blueberries", "boysenberries", "strawberries", "raspberries"))], na.rm=TRUE) * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="misc vegs & fruits", mean(cropdollar$dollar_kg[which(cropdollar$Crop %in% c("broccoli","lettuce","cabbage","celery","cucumbers","cherries","nectarines","peaches"))], na.rm=TRUE) * in_df$dir_a, in_df$dollar)
  in_df$dollar <- ifelse(in_df$Crop=="other tree crops", mean(cropdollar$dollar_kg[which(cropdollar$Crop %in% c("almonds","cherries","nectarines","peaches"))], na.rm=TRUE) * in_df$dir_a, in_df$dollar)
  
  # Add double crop sums
  dbl_crop_func <- function(crop_1, crop_2){
    # If crop is double crop, dollar is sum of dollar_m2 * dir_a (again, multiplying the direct area times the number of years through end year of analysis gives cumulative impact of each array to end year)
    in_df$dollar <- ifelse(in_df$FARS_crop==paste(crop_1, crop_2, sep = "_"), (cropdollar$dollar_kg[which(cropdollar$Crop==crop_1)] + 
                                                                               cropdollar$dollar_kg[which(cropdollar$Crop==crop_2)]) * in_df$dir_a, in_df$dollar)
    return(in_df) }
  # Run double crop function -- These are all possible double crop pairings
  in_df <- dbl_crop_func("barley", "corn")
  in_df <- dbl_crop_func("barley", "sorghum")
  in_df <- dbl_crop_func("lettuce", "barley")
  in_df <- dbl_crop_func("lettuce", "cotton")
  in_df <- dbl_crop_func("lettuce", "wheat")
  in_df <- dbl_crop_func("oats", "corn")
  in_df <- dbl_crop_func("wheat", "corn")
  in_df <- dbl_crop_func("wheat", "cotton")
  in_df <- dbl_crop_func("wheat", "sorghum")
  
  # Add double crop sums for crops not included in crop surveys and in the derived dataframe "cropdollar"
  dbl_crop_funcAUX <- function(crop_1, crop_2, crop_3, crop_4){
    # If crop is double crop, dollar is sum of dollar_m2 * dir_a (again, multiplying the direct area times the number of years through end year of analysis gives cumulative impact of each array to end year)
    in_df$dollar <- ifelse(in_df$FARS_crop==paste(crop_1, crop_2, sep = "_"), (cropdollar$dollar_kg[which(cropdollar$Crop==crop_3)] + 
                                                                                 cropdollar$dollar_kg[which(cropdollar$Crop==crop_4)]) * in_df$dir_a, in_df$dollar)
    return(in_df) }
  in_df <- dbl_crop_funcAUX("barley", "soybeans", "barley", "beans") # for instance, no soybeans, so double crop use beans value
  in_df <- dbl_crop_funcAUX("corn", "soybeans", "corn", "beans")
  in_df <- dbl_crop_funcAUX("soybeans", "cotton", "beans", "cotton")
  in_df <- dbl_crop_funcAUX("soybeans", "oats", "beans", "oats")
  in_df <- dbl_crop_funcAUX("wheat", "soybeans", "wheat", "beans")
  in_df <- dbl_crop_funcAUX("small grains", "corn", "wheat", "corn")
  in_df <- dbl_crop_funcAUX("lettuce", "cantaloupes", "lettuce", "melons")
  in_df <- dbl_crop_funcAUX("lettuce", "cantaloupes", "lettuce", "melons")
  in_df <- dbl_crop_funcAUX("lettuce", "cantaloupes", "lettuce", "melons")
  
  # This is only for us to see if we have missed crops -- Failsafe that worked
  for(i in c(1:nrow(in_df))){
    in_df[i, ]$dollar <- ifelse(min(in_df$dollar)==0, NA, in_df[i, ]$dollar)
  }
  
  # -------------------------------------------- #
  # Group and calculate total offset and project # 
  # -------------------------------------------- #
  # Group by year and create new column for $/m2
  Food_D <- in_df %>% 
    group_by(Year) %>% 
    summarise(dollar = sum(dollar, na.rm = FALSE), dir_a = sum(dir_a, na.rm=TRUE)) # NA.RM = False lets us see if we have missed crops
  
  # Description is from regional assessment. Ignore, but kept for reference
  # Find total area used for solar in each year and project (This initial cumsum is because land converted to solar in one year continues to be converted the next)
  # HOWEVER -- Above, by taking each crop of each array and projecting/sum forward through the end year, we have already made the Food_D$dollar cumulative through 2018, this must be accounted for
  #Food_D$dollar = Food_D$dollar %>% cumsum() # -- Already cumulative from above code, however, still accumlate post end year
  Food_D <- getTemporal_Projection(input_vector=Food_D$dollar, input_name='dollar', vector_year=vector_year, accumulate_post_endYear=TRUE)
  
  # Assume that food prices will scale directly with energy prices
  energy_price_growth_rate_temp <- c(energy_price_growth_rate, rep(mean(energy_price_growth_rate, na.rm=TRUE), 100)) %>% cumsum() # extends average growth rate by 100 years
  energy_price_growth_rate_temp <- energy_price_growth_rate_temp[1:system_lifespan]
  for(i in c((install_period+1):model_length)){
    Food_D[i, ]$dollar <- Food_D[i, ]$dollar*(1+energy_price_growth_rate_temp[i-install_period])
  }
  
  # Create food vectors, for right now just vary yield by 20%
  Food_D$food_min <- Food_D$dollar*food_yield_best
  Food_D$food_max <- Food_D$dollar*food_yield_worst
  Food_D$food_base <- Food_D$dollar*food_yield_base
  Food_D$dollar = NULL
  return(Food_D)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Current version work -- only run for "commercial" arrays, or where modeled electricity returns exist (variables in input dataset of "Econ_20##")
getEconomic_energy = function(in_df){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  
  # ------------------------------ #
  # Get total electricity produced # -- Also acquires mono/mulit silicate scenarios and adjust for inflation 
  # ------------------------------ #
  # Get energy produced over lifetime of array
  produced <- getResource_energyProduced(in_df)
  # Get multiplier for scenarios
  scenario_mult <- produced
  scenario_mult$energy_max <- scenario_mult$energy_max / scenario_mult$energy_base
  scenario_mult$energy_min <- scenario_mult$energy_min / scenario_mult$energy_base
  scenario_mult$energy_base <- rep(1, nrow(scenario_mult))
  # Get Econ Variables Only
  indx <- in_df$Index
  yr_inst <- in_df$Year
  in_df$geometry = NULL
  crop <- in_df$Crp_n_1 %>% tolower()
  energy_econ <- in_df %>% dplyr::select(matches("Econ")) %>% t() %>% as.data.frame()
  colnames(energy_econ) <- c("energy_base")
  # Get year list from col names
  year_list <- rownames(energy_econ) %>% parse_number()
  energy_econ$Year <- year_list # list of modeled years -- continous
  # Adjust for inflation -- prices modeled as nominal -- CPI because in this case, farmers are consumers 
  for(Year in c(unique(energy_econ$Year))){
    infl_rate_mult <- inflation_rateCPI$infl_rate_mult[which(inflation_rateCPI$Year==Year)]
    energy_econ$energy_base <- ifelse(energy_econ$Year==Year, energy_econ$energy_base * infl_rate_mult, energy_econ$energy_base)
  }
  
  # ---------------------------------------------------------------------------------------------------------- #
  # Get projected prices based on final model year price by EIA exspected rate of change in electricity prices #
  # ---------------------------------------------------------------------------------------------------------- #
  # Get electric rate (per $/GWh)
  energy_econ_modelTailRate <- ( energy_econ$energy_base %>% tail(1) %>% as.numeric() ) / ( in_df$Gen_2018*kW_MW %>% as.numeric() ) # $/kWh, modeled generation in MWh so have to convert to kWh
  energy_econ_modelTailRate <- energy_econ_modelTailRate * kW_GW # $/GWh -- the output of getResource_energyProduced() is in GWh

  # Temporal changes in price of electricity
  energy_price_growth_rate_temp <- c(energy_price_growth_rate, rep(mean(energy_price_growth_rate, na.rm=TRUE), 100)) %>% cumsum() # extends average growth rate by 100 years for scenarios
  energy_price_growth_rate_temp <- energy_price_growth_rate_temp[1:system_lifespan] # Energy growth rate from start year of growth projection through lifespan
  # Apply compound interest to adjust expected rate changes in electricity
  #projected_price = data.frame(price = numeric(), Year = integer())
  #for(n in c(1:system_lifespan)){
  #  projected_price_temp <- c(energy_econ_modelTailRate*(1+energy_price_growth_rate_temp[n])^n)
  #  year_temp = tail(year_list, 1)+n-1
  #  df_temp = data.frame(price = projected_price_temp, Year = year_temp)
  #  projected_price <- rbind(projected_price, df_temp)
  #}
  # NON COMPOUND INTEREST - adjust expected rate changes in electricity -- rate of annual change is the current input
  projected_price_temp <- c(energy_econ_modelTailRate*(1+energy_price_growth_rate_temp))
  projected_price <- data.frame(price = projected_price_temp, Year = c(tail(year_list, 1):(tail(year_list, 1)+system_lifespan-1)))
  
  # ---------------------------------------- #
  # Project Econ by Rate of Final Model year # -- And by EIA espected rate of change in electricity prices
  # ---------------------------------------- #
  # Get cumulative sum for array, and prepare dataframe for projection
  energy_econ$energy_base <- energy_econ$energy_base %>% cumsum()
  econ_yrs_start <- which(rownames(energy_econ) %>% parse_number()==yr_inst) %>% as.numeric() # Get first row of generation & First year-row to project from
  energy_econ[c((nrow(energy_econ)+1):(econ_yrs_start+system_lifespan-1)),] <- 0 # Saves all modeled generation and generates empty df for system lifespan
  energy_econ <- energy_econ[c(econ_yrs_start:nrow(energy_econ)), ] %>% as.data.frame() # Limit to first year of generation and Return df
  energy_econ$Year <- c(yr_inst:(yr_inst+system_lifespan-1))
  
  # Steal GWh production from getResource_EnergyProduced() function called above, and use the temporally adjusted rate ($/GWh) to get dollar amount in a give year
  projected_price <- projected_price[which(projected_price$Year<=max(energy_econ$Year)), ]
  energy_econ[which(energy_econ$Year > tail(year_list, 1)), ]$energy_base <- produced$energy_base[which(produced$Year > tail(year_list, 1))] * projected_price$price[which(projected_price$Year > tail(year_list, 1))] # for year after model is finsihed
  
  # ---------------------------------------- #
  #  Apply electricity production scenarios  # 
  # ---------------------------------------- # 
  energy_econ$energy_max <- energy_econ$energy_base * scenario_mult$energy_max
  energy_econ$energy_min <- energy_econ$energy_base * scenario_mult$energy_min
  energy_econ$energy_base <- energy_econ$energy_base * scenario_mult$energy_base
  rownames(energy_econ) <- c(1:nrow(energy_econ))
  # Return df
  return(energy_econ)
}

# Economic implications of change in energy budget
getEconomic_energy_SingleRate = function(in_df){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  
  # ------------------------------ #
  # Get total electricity produced # 
  # ------------------------------ #
  # Only for commerical arrays 
  produced <- getResource_energyProduced(in_df)
  
  # ---------------------------------- #
  # Calculate value of produced energy # 
  # ---------------------------------- #
  # Assume temporal changes in price of electricity
  energy_price_growth_rate_temp <- c(energy_price_growth_rate, rep(mean(energy_price_growth_rate, na.rm=TRUE), 100)) %>% cumsum() # extends average growth rate by 100 years
  energy_price_growth_rate_temp <- energy_price_growth_rate_temp[1:system_lifespan]
  energy_price <- elec_price * kW_GW # get into GWh
  energy_price <- c(rep(energy_price, 10), energy_price*(1+energy_price_growth_rate_temp))
  energy_price <- energy_price[1:model_length]
  # Create final dataframe
  produced <- data.frame(Year = produced$Year, 
                         energy_base = (produced$energy_base*energy_price), 
                         energy_max = (produced$energy_max*energy_price), 
                         energy_min = (produced$energy_min*energy_price))
  # Return df
  return(produced)
}  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Economic implications of change in water use -- ECON exists is true if economic modeling has been performed (in_solar_df contains "Econ_20##" variables). Otherwises, uses set average.
getEconomic_water = function(in_df, econ_exist = TRUE){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  
  # ------------------------------------ #
  # Get total change in water use budget # -- And energy produced, and cost of energy
  # ------------------------------------ #
  irrig_saved <- getResource_irrigEnergySaved(in_df) # Run water budget effect functions
  energy_produced <- getResource_energyProduced(in_df) # Get energy produced over lifetime of array
  if(econ_exist==TRUE){
    energy_cost <- getEconomic_energy(in_df)
  } else {
    energy_cost <- getEconomic_energy_SingleRate(in_df) # Get value of energy produced, use single rate if econ model does not exist for array
  } 
  
  # ----------------------------------------- #
  # Calculate value of conserved irrig energy # 
  # ----------------------------------------- #
  # Get annual electricity rates
  energy_price <- energy_cost$energy_base / energy_produced$energy_base 
  
  # Adjust water right by inflation
  waterRightDollar_m3_temp = waterRightDollar_m3 * inflation_rateCPI$infl_rate_mult[which(inflation_rateCPI$Year==waterRightDollar_year)]
  
  # Calculate value of conserved irrig energy WITH water right contract per m3
  water_base = (irrig_saved$irrig_fg*energy_price + irrig_saved$irrig_fg*waterRightDollar_m3_temp) 
  water_max = (irrig_saved$irrig_fg_dry*energy_price + irrig_saved$irrig_fg_dry*waterRightDollar_m3_temp)
  water_min = (irrig_saved$irrig_fg_wet*energy_price + irrig_saved$irrig_fg_wet*waterRightDollar_m3_temp)
  
  # Combine into df
  irrig_saved <- data.frame(Year = irrig_saved$Year, water_base = water_base, water_max = water_max, water_min = water_min)
  return(irrig_saved)
} 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Economic implications of Operating and Maintaining PV
getEconomic_OandM = function(in_df){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  
  # ------------------ #
  # Get historical ATB # 
  # ------------------ #
  # Get O&M Costs from ATB Database
  NREL_ATB_OandM <- NREL_ATB_future[grep("O&M", NREL_ATB_future$core_metric_parameter), ] # Fixed and Variable O&M Costs
  NREL_ATB_OandM <- NREL_ATB_OandM[grep("CommPV", NREL_ATB_OandM$technology), ] # relevant for this model, get only commercial scale O&M
  
  # Get Cost Return Period Closest to System Lifespan model variable
  crpyears <- unique(NREL_ATB_OandM$crpyears) %>% as.character() %>% as.numeric()
  near_crpyear = crpyears[which(abs(system_lifespan-crpyears)==min(abs(system_lifespan-crpyears)))] %>% min() # Get closest census year to installation (minmum will be most concervative)
  NREL_ATB_OandM <- NREL_ATB_OandM[grep(near_crpyear, NREL_ATB_OandM$crpyears), ] # Fixed and Variable O&M Costs
  
  # Group across years and scenarios and correct to MW
  NREL_ATB_OandM <- NREL_ATB_OandM %>% group_by(core_metric_variable, scenario) %>% summarise(value = mean(value, na.rm=TRUE)*kW_GW) # $/kW/yr to $/GW/yr
  OandM_cost_df <- data.frame(Year = unique(NREL_ATB_OandM$core_metric_variable), 
                              Conservative = NREL_ATB_OandM$value[which(NREL_ATB_OandM$scenario=="Conservative")], 
                              Moderate = NREL_ATB_OandM$value[which(NREL_ATB_OandM$scenario=="Moderate")], 
                              Advanced = NREL_ATB_OandM$value[which(NREL_ATB_OandM$scenario=="Advanced")])
  # Get model years
  model_years = c(start_install_year:(end_install_year+system_lifespan-1))
  
  # Historically reported and modeled values -- Feldman et al., 2021; NREL, 2022a; Ramasamy et al., 2021
  hist_base = 18 # $/kWdc/yr
  hist_max  = 40 # $/kWdc/yr
  hist_min  = 0  # $/kWdc/yr
  
  # Extend bounds using start and end year (no assumptions of change)
  ext_yr = 250 # How many years to arbitrarily extend assumptions. Not adding any new data, just filling gap with bound values
  NREL_ATB_yrs = unique(OandM_cost_df$Year) %>% as.numeric()
  start_ext_yr = NREL_ATB_yrs[1] - ext_yr
  end_ext_yr = tail(NREL_ATB_yrs, 1) + ext_yr
  OandM_cost_df <- data.frame(Year = c(start_ext_yr:end_ext_yr), 
                              Conservative = c(rep(hist_max*kW_GW, ext_yr-1), OandM_cost_df$Conservative, rep(tail(OandM_cost_df$Conservative, 1), ext_yr+1)), 
                              Moderate = c(rep(hist_base*kW_GW, ext_yr-1), OandM_cost_df$Moderate, rep(tail(OandM_cost_df$Moderate, 1), ext_yr+1)), 
                              Advanced = c(rep(hist_min*kW_GW, ext_yr-1), OandM_cost_df$Advanced, rep(tail(OandM_cost_df$Advanced, 1), ext_yr+1)))
  #OandM_cost_df <- data.frame(Year = c(start_ext_yr:end_ext_yr), 
  #                            Conservative = c(rep(head(OandM_cost_df$Conservative, 1), ext_yr), OandM_cost_df$Conservative, rep(tail(OandM_cost_df$Conservative, 1), ext_yr)), 
  #                            Moderate = c(rep(head(OandM_cost_df$Moderate, 1), ext_yr), OandM_cost_df$Moderate, rep(tail(OandM_cost_df$Moderate, 1), ext_yr)), 
  #                            Advanced = c(rep(head(OandM_cost_df$Advanced, 1), ext_yr), OandM_cost_df$Advanced, rep(tail(OandM_cost_df$Advanced, 1), ext_yr)))
  
  # Limit to years of interest
  OandM_cost_df <- OandM_cost_df[which(OandM_cost_df$Year %in% model_years),]
  
  # Adjust for Inflation
  OandM_USD_adj <- inflation_ratePPI$infl_rate_mult[which(inflation_ratePPI$Year==NREL_ATB_yr)]
  OandM_cost_df[, c(2:4)] <- OandM_cost_df[, c(2:4)] * OandM_USD_adj
                              
  # ----------------------------------------- #
  # Get total O&M Cost from Temporal Capacity # 
  # ----------------------------------------- #
  # Operations and Maintenance dollar value
  Capacity <- getTemporal_Capacity(in_df)
  OandM_D <- OandM_cost_df$Moderate * Capacity$Capacity %>% cumsum() 
  OandM_Dmax <- OandM_cost_df$Conservative * Capacity$Capacity %>% cumsum() 
  OandM_Dmin <- OandM_cost_df$Advanced * Capacity$Capacity %>% cumsum()
  OandM_df <- data.frame(Year = Capacity$Year, oandm_base = OandM_D, oandm_max = OandM_Dmax, oandm_min = OandM_Dmin)
  return(OandM_df)
} 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Economic implications of Installing PV
getEconomic_install = function(in_df){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  
  # ------------------------------- #
  # Get Installation Costs from ATB # 
  # ------------------------------- #
  inst_USD_adj <- inflation_rateCPI$infl_rate_mult[which(inflation_rateCPI$Year==inst_USD_yr)] # https://atb.nrel.gov/electricity/2022/commercial_pv
  
  # Set variables
  commercial_instCost_base <- install_cost_df$Value[which(install_cost_df$Year %in% c(start_install_year:end_install_year) & install_cost_df$Type=="median")] / 1000 * inst_USD_adj * (1-solar_ITC)
  commercial_instCost_high <- install_cost_df$Value[which(install_cost_df$Year %in% c(start_install_year:end_install_year) & install_cost_df$Type=="80th percentile")] / 1000 * inst_USD_adj * (1-solar_ITC)
  commercial_instCost_low  <- install_cost_df$Value[which(install_cost_df$Year %in% c(start_install_year:end_install_year) & install_cost_df$Type=="20th percentile")] / 1000 * inst_USD_adj * (1-solar_ITC)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Just for comparing
  
  # Comparison to NREL cost benchmark values -- 2020 USD
  NREL_CB <- c(6.67, 6.13, 4.08, 3.26, 3.21, 2.64, 2.50, 2.07, 1.84) %>% rev() * inflation_rateCPI$infl_rate_mult[which(inflation_rateCPI$Year==2020)] # Ramassamay et al., 2021
  TrackSun <- install_cost_df$Value[which(install_cost_df$Year %in% c(2010:2018) & install_cost_df$Type=="median")] / 1000 * inst_USD_adj # Barbose et al., 2021
  mean_dif <- mean(TrackSun - NREL_CB)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Just for comparing
  
  # Get commercial base, 80th percentile, and 20th percentile prices
  base_vect_temp <- commercial_instCost_base
  max_vect_temp <- commercial_instCost_high
  min_vect_temp <- commercial_instCost_low
  # Run capacity calc
  Capacity <- getTemporal_Capacity(in_df)
  
  # Installation cost -- convert from $/kw to $/GW for monosilicate prices, multisilicate  prices, and averages
  instcost_D <- c(c(Capacity$Capacity[1], diff(Capacity$Capacity)) * base_vect_temp*watt_GW) %>% cumsum()
  instcost_D <- c(instcost_D[1:system_lifespan]*(1-solar_ITC), rep(instcost_D[system_lifespan]*(1-solar_ITC), install_period))
  instcost_Dmax <- c(c(Capacity$Capacity[1], diff(Capacity$Capacity)) * max_vect_temp*watt_GW) %>% cumsum()
  instcost_Dmax <- c(instcost_Dmax[1:system_lifespan]*(1-solar_ITC), rep(instcost_Dmax[system_lifespan]*(1-solar_ITC), install_period))
  instcost_Dmin <- c(c(Capacity$Capacity[1], diff(Capacity$Capacity)) * min_vect_temp*watt_GW) %>% cumsum()
  instcost_Dmin <- c(instcost_Dmin[1:system_lifespan]*(1-solar_ITC), rep(instcost_Dmin[system_lifespan]*(1-solar_ITC), install_period)) 
  install_df <- data.frame(Year = Capacity$Year, 
                           install_base = instcost_D, # - interconnection_fee_df$base
                           install_max = instcost_Dmax, # - interconnection_fee_df$worst
                           install_min = instcost_Dmin) #  - interconnection_fee_df$best
  return(install_df)
} 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Economic implications of Land Leases for Utility Scale Arrays
getEconomic_landlease = function(in_df){
  # ------------------------------------------------- #
  # Get temporal area and calculate land lease impact # 
  # ------------------------------------------------- #
  # Get temporal area
  landlease_df <- getTemporal_Area(in_df)
  # Include temporal increases in land lease rates
  landlease_price_growth_rate_temp <- landlease_price_growth_rate %>% cumsum()
  
  # Multiply by land lease rates
  landlease_df$lease_base <- landlease_df$area * landlease_per_m2_avg * (1+landlease_price_growth_rate_temp) 
  landlease_df$lease_max <- landlease_df$area * landlease_per_m2_max * (1+landlease_price_growth_rate_temp)
  landlease_df$lease_min <- landlease_df$area * landlease_per_m2_min * (1+landlease_price_growth_rate_temp)
  landlease_df$area = NULL
  return(landlease_df)
} 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Economic implications of offsetting annual agricultural operation costs per square meter
getEconomic_operationalCost = function(in_df){
  # ------------------------- #
  #    Get Non User Inputs    # 
  # ------------------------- #
  # Non user variables 
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  
  # SETUP Generated Variables -- no input necessary
  install_period = end_install_year - start_install_year
  model_length = install_period + system_lifespan 
  
  # Set Variables for yield through time
  food_yield_base = c( (1-yield_deficit_base) + yield_time_base * c(1:(model_length)))
  food_yield_best = c( (1-yield_deficit_best) + yield_time_best * c(1:(model_length)))
  food_yield_worst= c( (1-yield_deficit_worst) + yield_time_worst * c(1:(model_length)))
  
  # -------------------------------------------- #
  # Generate cost/m2 dataframe with cost/m2 data # -- Only if operaiontal_cost_m2.csv has not been previously defined
  # -------------------------------------------- #
  if(file.exists("Data/Derived/operational_cost_m2.csv")){
    # If file already generated in previous model run, call in file
    cost_sqm <- read.csv("Data/Derived/operational_cost_m2.csv", stringsAsFactors=FALSE)
  } else {
    # If file has not been generated, fun FEWLS_kcalperm2.R to generate, then call in
    source("FEWLS_cropCost.R")
    cost_sqm <- read.csv("Data/Derived/operational_cost_m2.csv", stringsAsFactors=FALSE)
  }
  
  # ------------------------------- #
  # Get crop specific temporal area # 
  # ------------------------------- #
  # Run function
  df <- getCropTemporal_Area(in_df)
  # Get crop operational costs from operational_cost_m2.csv and multiply
  df$cost_m2 = 0
  df$cost_m2 <- cost_sqm$cost_m2[match(df$Crop, cost_sqm$Crop)]
  df$cost <- df$dir_a * df$cost_m2
  # Group by year
  df <- df %>% group_by(Year) %>% summarize(operation_base = sum(cost, na.rm = TRUE))
  
  # -------------------------------- #
  # Get operational offset scenarios # 
  # -------------------------------- #
  # Growth in operational costs based on changed in electricity requirements
  energy_price_growth_rate_temp <- c(energy_price_growth_rate, rep(mean(energy_price_growth_rate, na.rm=TRUE), 100)) %>% cumsum() # extends average growth rate by 100 years
  energy_price_growth_rate_temp <- energy_price_growth_rate_temp[1:system_lifespan]
  for(i in c((install_period+1):model_length)){
    df[i, ]$operation_base <- df[i, ]$operation_base*(1+energy_price_growth_rate_temp[i-install_period])
  }
  
  # Projected increase in operational costs of ag (scenarios based on variation in yield)
  food_operational_cost_base <- (food_yield_base) 
  food_operational_cost_worst <- (food_yield_best) # These are flipped because offset food yield is a negative but offset operational cost is a positive
  food_operational_cost_best <- (food_yield_worst) # These are flipped because offset food yield is a negative but offset operational cost is a positive

  # Calculate operational costs through time with projection changes
  df$operation_max <- (df$operation_base * food_operational_cost_best)  
  df$operation_min <- (df$operation_base * food_operational_cost_worst) 
  df$operation_base <- (df$operation_base * food_operational_cost_base)
  #for(i in c((install_period+2):model_length)){
  #  df[i, ]$operation_min <- df[i, ]$operation_base*(1+food_operational_cost_best_temp[i-install_period-1])
  #  df[i, ]$operation_max <- df[i, ]$operation_base*(1+food_operational_cost_worst_temp[i-install_period-1])
  #  df[i, ]$operation_base <- df[i, ]$operation_base*(1+food_operational_cost_base_temp[i-install_period-1])
  #}
  return(df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Resource implications of agrisolar co-location -- non-scale specific
getResource_CommUtil = function(in_df){
  # ----------------------------------------------------- #
  # Generate all economic variables for commercial arrays # 
  # ----------------------------------------------------- #
  # Get cumulative resource predictions
  kcal <- getResource_foodkcal(in_df) # FOOD in kcal 
  energy <- getResource_energyProduced(in_df) # ENERGY PRODUCED in (MWh)
  irrEngy <- getResource_irrigEnergySaved(in_df) # ENERGY CONSERVED FROM IRRIGATION in (MWh)
  irrig <- getResource_irrigWater(in_df) # WATER  CONSERVED FROM IRRIGATION (m3)
  OandM <- getResource_OandMWaterUsed(in_df, TRUE) # WATER USED FOR O&M -- For all arryas, consider_nonIrrig == TRUE (m3)
  OandM_irrig <- getResource_OandMWaterUsed(in_df) # WATER USED FOR O&M -- For only irrigated arrays, consider_nonIrrig == FALSE (m3)
  
  # Change names of columns to be unique
  names(irrEngy) <- c("Year", "irrEngy_fg", "irrEngy_fg_dry", "irrEngy_fg_wet")
  names(OandM_irrig) <- c("Year", "oandmIr_base", "oandmIr_min", "oandmIr_max")
  
  # Save as single df and remove duplicate year columns
  df <- cbind(kcal, energy, irrEngy, irrig, OandM, OandM_irrig)
  df <- df[, !duplicated(colnames(df))]
  
  # Save valuable variables
  df$Index = in_df$Index
  df$Crop = in_df$Crp_n_1 %>% tolower()
  df$Yr_inst = in_df$Year # Have to have different variable name for installations year since year is also year of model
  df$Capacity = in_df$Capacity
  return(df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to adjust invividual resources to NPV with DCF for continous costs/profits
applyDCF <- function(df, discount_rate_nominal){
  # First, calculate real discount rate -- Fisher equation (1896)
  discount_rate_real = ( (1+discount_rate_nominal) / (1+inflation_rate_future) ) - 1
  #For each column, apply discount rate to each year's cashflow contribution
  for(j in c(1:ncol(df))){
    vector <- c(df[,j][1], diff(df[,j]))
    for(i in c(1:length(vector))){ # initial year is not discounted
      # Calcualte discount
      vector[i] <- vector[i] / ((1+discount_rate_real)^i)}
    # Cumsum to return a cumulative value
    df[,j] <- vector %>% cumsum()
  }
  return(df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Economic implications of commercial scale installations (Installation-Based Function)
getEconomic_Commercial = function(in_df){
  # ----------------------------------------------------- #
  # Generate all economic variables for commercial arrays # 
  # ----------------------------------------------------- #
  # Get cumulative economic resource predictions
  food <- getEconomic_food(in_df) 
  OandM <- getEconomic_OandM(in_df)
  install <- getEconomic_install(in_df)
  operation <- getEconomic_operationalCost(in_df)
  
  # Special case for electricity economics, and water economics (since it requries electricity price). For out dataset, we only calculated economic values for commercial scale installations (<1MW). If running a scenraio with a different threshold, have to use single electric rate value. Also, for running sensitivity analysis
  # First, see if electricity return column(s) exist for array
  econ_exist_temp = in_df %>% dplyr::select(matches("Econ"))
  econ_exist_temp$geometry = NULL
  econ_exist_temp = econ_exist_temp %>% rowSums() %>% as.numeric()
  econ_exist_temp <- ifelse(is.na(econ_exist_temp), 0, econ_exist_temp)
  # Then call in appropriate functions and settings to run based on presence of econ variables and desired output
  if(econ_exist_temp > 0 & use_single_elecRate == FALSE){
    energy <- getEconomic_energy(in_df)
    water <- getEconomic_water(in_df, econ_exist = TRUE)
  } else {
    energy <- getEconomic_energy_SingleRate(in_df)
    water <- getEconomic_water(in_df, econ_exist =  FALSE)
  }
  
  # Bind and   # Remove duplicate year columns
  df <- cbind(food, energy, water, OandM, install, operation)
  df <- df[, !duplicated(colnames(df))]
  
  # Dicount rate function (continous costs/profits) -- remove install because that is cost of initial investment (day 0, so no discount)
  df[, which(!names(df) %in% c("Year", "install_base", "install_max", "install_min"))] <- applyDCF(df = df[, which(!names(df) %in% c("Year", "install_base", "install_max", "install_min"))], discount_rate_nominal = discount_rate_nominal)
  
  # Calculate total budget and return index, crop, and year of installation
  df$Tot_Budg_base <- df$energy_base + df$water_base + df$operation_base - df$food_base - df$oandm_base - df$install_base 
  df$Tot_Budg_max <- df$energy_max + df$water_max + df$operation_max - df$food_min - df$oandm_min - df$install_min 
  df$Tot_Budg_min <- df$energy_min + df$water_min + df$operation_min - df$food_max - df$oandm_max - df$install_max 
  
  # Save valuable variables
  df$Index = in_df$Index
  df$Crop = in_df$Crp_n_1 %>% tolower()
  df$Yr_inst = in_df$Year # Have to have different variable name for installations year since year is also year of model
  df$Capacity = in_df$Capacity
  return(df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Economic implications of utility scale installations
getEconomic_Utility = function(in_df){
  # ----------------------------------------------------- #
  # Generate all economic variables for commercial arrays # 
  # ----------------------------------------------------- #
  # Get cumulative economic resource predictions
  food <- getEconomic_food(in_df)
  water <- getEconomic_water(in_df, econ_exist = FALSE) # econ exist is set to false because economic modeling was not performed for utility scale arrays (inappropriate rates for size of system) -- uses commercial average
  operation <- getEconomic_operationalCost(in_df)
  landlease <- getEconomic_landlease(in_df)
  df <- cbind(food, water, operation, landlease)
  # Remove duplicate year columns
  df <- df[, !duplicated(colnames(df))]
  
  # Dicount rate function
  df[, which(!names(df) %in% c("Year"))] <- applyDCF(df = df[, which(!names(df) %in% c("Year"))], discount_rate_nominal = discount_rate_nominal)
  
  # Calculate total budget
  df$Tot_Budg_base <- df$lease_base + df$water_base + df$operation_base - df$food_base
  df$Tot_Budg_max <- df$lease_max + df$water_max + df$operation_max - df$food_min
  df$Tot_Budg_min <- df$lease_min + df$water_min + df$operation_min - df$food_max
  
  # Save valuable variables
  df$Index = in_df$Index
  df$Crop = in_df$Crp_n_1 %>% tolower()
  df$Yr_inst = in_df$Year # Have to have different variable name for installations year since year is also year of model
  df$Capacity = in_df$Capacity
  return(df)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get annual on-farm operational requirements (MWh load)
getResource_annFarmOperationReq = function(in_df){
  # ------------- #
  # Get farm size # 
  # ------------- #
  # Only for commerical arrays -- check FEWLS_FarmElecBudgExport to check that theres not a erronious capacity_threshold
  in_df <- in_df[which(in_df$Capacity<=capacity_threshold), ]
  
  # Get counties for each array
  counties <- counties(state, TRUE) %>% st_as_sf() %>% st_set_crs(st_crs(in_df))
  intersection <- st_intersection(st_centroid(in_df), counties)
  in_df$County <- intersection$NAME[match(in_df$Index, intersection$Index)]
  in_df$COUNTYFP <- counties$COUNTYFP[match(in_df$County, counties$NAME)]
  in_df$COUNTYFP <- gsub("(?<![0-9])0+", "", in_df$COUNTYFP, perl = TRUE) # Extracts non-zero led countyfp numerically from a string
  rm(counties)
  
  # Get county level average farm size
  state_county_agarea <- county_aglandarea[which(tolower(county_aglandarea$State)==tolower(state)), ]
  ag_farmsize <- state_county_agarea[which(state_county_agarea$Data.Item=="AG LAND - ACRES"), ]
  ag_farmsize$Value <- as.numeric(gsub(",","",ag_farmsize$Value))
  ag_farmsize$area <- ag_farmsize$Value * acres_m2
  ag_farmsize$farms <- state_county_agarea$Value[which(state_county_agarea$Data.Item=="AG LAND - NUMBER OF OPERATIONS")]
  ag_farmsize$farms <- as.numeric(gsub(",","", ag_farmsize$farms))
  ag_farmsize$avg_area <- ag_farmsize$area / ag_farmsize$farm
  ag_farmsize <- ag_farmsize[which(!is.na(ag_farmsize$avg_area)), ]
  
  # Save county average farm size to array based on county and survey year closest to installation -- farm_size here results in m2
  census_yrs = unique(ag_farmsize$Year) %>% as.numeric() 
  in_df$farm_size = median(ag_farmsize$avg_area[which(ag_farmsize$Year==max(census_yrs, na.rm = TRUE))], na.rm = TRUE) # save avearge farm size across state in latest census yearincase any missing data
  for(i in c(1:nrow(in_df))){
    county = in_df[i, ]$County %>% tolower()
    install_yr = in_df[i, ]$Year
    near_census_yr = census_yrs[which(abs(install_yr-census_yrs)==min(abs(install_yr-census_yrs)))] # Get closest census year to installation
    avg_farm_size = ag_farmsize$avg_area[which(as.numeric(ag_farmsize$Year)==near_census_yr & tolower(ag_farmsize$County)==county)]
    avg_farm_size = ifelse(length(avg_farm_size)==0, median(ag_farmsize$avg_area[which(tolower(ag_farmsize$County)==county)], na.rm = TRUE), avg_farm_size) # If nearest census year missing, choose most recent county census year
    in_df[i, ]$farm_size <- ifelse(length(avg_farm_size)==0, in_df[i, ]$farm_size, avg_farm_size) # If still missing, use state wide average
  }
  
  # -------------------------- #
  # Get Farm Irrigation Energy # -- (Ei) Annual on-farm irrigation energy requirments: For every farm regardless of irrigation status, Ei is crop and spatially relevent
  # -------------------------- #
  # Get irrig energy per m2
  irrigEnergy_m2 <- data.frame(Index = integer(), irrig_fg = numeric(), irrig_fg_dry = numeric(), irrig_fg_wet = numeric())
  for(i in c(1:nrow(in_df))){
    irrigEnergy_m2_temp <- getResource_irrigWater(in_df[i, ], irrig_energy = TRUE, irrig_energy_per_m2 = TRUE, include_oandm = FALSE) # get irrig_fg (GWh) and area (m2) for each array, and does NOT include O&M water use additions because we are looking for a irrigation energy density here
    irrigEnergy_m2 <- rbind(irrigEnergy_m2, irrigEnergy_m2_temp) 
    print(paste("System Index: ", i, sep="")) }
  irrigEnergy_m2$GWh_m2 <- irrigEnergy_m2$irrig_fg / irrigEnergy_m2$area * MW_GW # To GWh
  in_df$geometry = NULL
  irrigEnergy_m2 <- merge(in_df, irrigEnergy_m2, by="Index")
  # Get farm level estimate of irrigating the average county level farm
  irrigEnergy_m2$FarmE_irrig_GWh <- irrigEnergy_m2$GWh_m2 * irrigEnergy_m2$farm_size # (Ei) in GWh
  
  # --------------------- #
  # Get Farm Other Energy # -- (Eo) Annual on-farm non-irrigation energy requiremetns: For every farm, Eo is only crop releveant and attempt to remove spatial importance
  # --------------------- #
  # Get crop rotation and split annual farm irrigation energy requirements across crops in rotation
  nonirrigEnergy_m2 <- getCrop_rotation(irrigEnergy_m2)
  nonirrigEnergy_m2$FarmE_irrig_GWh <- nonirrigEnergy_m2$FarmE_irrig_GWh * nonirrigEnergy_m2$crop_prop # because cropRotation function only splits direct area by crop_prop
  nonirrigEnergy_m2$GWh_m2 <- nonirrigEnergy_m2$GWh_m2 * nonirrigEnergy_m2$crop_prop
  # Set starting value for Eo
  nonirrigEnergy_m2$FarmE_other_GWh = 0
  # For every crop, calculate the average Ei
  unique_crops = unique(nonirrigEnergy_m2$unique_crops) 
  for(crop in unique_crops){
    # Subset for every crop
    Eo_temp <- nonirrigEnergy_m2[which(nonirrigEnergy_m2$unique_crops==crop), ] %>% as.data.frame() # Select for crop
    Ei_avg <- median(Eo_temp$GWh_m2, na.rm = TRUE) # Get Irrigation Energy Req spatial density (same as Ei / AreaFarm ~ EarrayIrrig / ArArray)
    # Correct for proportion of farm level energy contributed to irrigation
    IrrigEProp <- (1-irrig_energymix_prop) / irrig_energymix_prop
    # Calculate Eo
    nonirrigEnergy_m2$FarmE_other_GWh[which(nonirrigEnergy_m2$unique_crops==crop)] <- Ei_avg * IrrigEProp * nonirrigEnergy_m2$farm_size
  }
  # Group by array index and sum Eo for each crop within rotation -- then add to irrig dataframe
  nonirrigEnergy_m2 <- nonirrigEnergy_m2 %>% group_by(Index) %>% summarize(FarmE_other_GWh = sum(FarmE_other_GWh, na.rm = TRUE))
  irrigEnergy_m2$FarmE_other_GWh <- nonirrigEnergy_m2$FarmE_other_GWh[match(irrigEnergy_m2$Index, nonirrigEnergy_m2$Index)]
  
  # If irrigated (Et = Ei + Eo). If non irrigated (Et = Eo)
  irrigEnergy_m2$ann_operation_MWh <- ifelse(irrigEnergy_m2$irrig == 1, irrigEnergy_m2$FarmE_irrig_GWh + irrigEnergy_m2$FarmE_other_GWh, irrigEnergy_m2$FarmE_other_GWh)
  #irrigprop = mean(irrigEnergy_m2$FarmE_irrig_GWh / irrigEnergy_m2$ann_operation_MWh) # 65% for irrig arrays
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Values not saved on export, but important to report
  
  # Susbset for irrig vs non irrig
  df <- irrigEnergy_m2
  df <- irrigEnergy_m2[which(irrigEnergy_m2$irrig==1),] # irrigated
  df <- irrigEnergy_m2[which(irrigEnergy_m2$irrig==0),] # non-irrigated
  
  # Get some values for estimated requirements pre-budget assessment
  annOppMWh_mean = median(df$ann_operation_MWh, na.rm=TRUE) # [which(df$Year==2018)] -- USED IN PAPER
  annOppMWh_summary = summary(df$ann_operation_MWh) # USED IN PAPER
  # Save for paper but not export
  prop2018 <- (df$Gen_2018 / df$ann_operation_MWh) #%>% mean()
  prop2017 <- (df$Gen_2017 / df$ann_operation_MWh) #%>% mean()
  prop2016 <- (df$Gen_2016 / df$ann_operation_MWh) #%>% mean()
  prop2015 <- (df$Gen_2015 / df$ann_operation_MWh) #%>% mean()
  prop2014 <- (df$Gen_2014 / df$ann_operation_MWh) #%>% mean()
  prop2013 <- (df$Gen_2013 / df$ann_operation_MWh) #%>% mean()
  prop2012 <- (df$Gen_2012 / df$ann_operation_MWh) #%>% mean()
  prop2011 <- (df$Gen_2011 / df$ann_operation_MWh) #%>% mean()
  prop2010 <- (df$Gen_2010 / df$ann_operation_MWh) #%>% mean()
  prop2009 <- (df$Gen_2009 / df$ann_operation_MWh) #%>% mean()
  prop2008 <- (df$Gen_2008 / df$ann_operation_MWh) #%>% mean()
  prop_all <- c(prop2018, prop2017, prop2016, prop2015, prop2014, prop2013, prop2012, prop2011, prop2010, prop2009, prop2008)
  prop_all <- prop_all[which(prop_all>0 & prop_all<Inf)]
  # Get number of overproducers
  mean_prop_all <- mean(prop_all)
  #mean_prop_2008 <- mean(prop2008[which(prop2008>0 & prop2008<Inf)])
  #mean_prop_2018 <- mean(prop2018[which(prop2018>0 & prop2018<Inf)])
  num_overprod <- length(prop2018[prop2018>1]) # USED IN PAPER
  overprod <- length(prop2018[prop2018>1]) / length(prop2018) 
  num_underprod <- length(prop2018[prop2018<1])
  underprod <- length(prop2018[prop2018<1]) / length(prop2018) # USED IN PAPER
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Values not saved on export, but important to report
  
  # For every array, calculate annual production value compared to predicted operation
  irrigEnergy_m2 <- irrigEnergy_m2 %>% dplyr::select(matches("Year|Gen|Index|irrig|ann_operation_MWh")) 
  irrigEnergy_m2$irrig_fg = NULL
  irrigEnergy_m2$irrig_fg_dry = NULL
  irrigEnergy_m2$irrig_fg_wet = NULL
  irrigEnergy_m2$FarmE_irrig_GWh = NULL
  df <- data.frame(Year = integer(), Index = integer(), irrig = integer(), ann_operation_MWh = numeric(), Gen = numeric())
  for(i in c(1:nrow(irrigEnergy_m2))){
    # Get continous vector for each array
    index <- irrigEnergy_m2[i, ]$Index
    array_temp <- irrigEnergy_m2[i, ]
    # Get generations and transpose into one column
    gen <- array_temp %>% dplyr::select(matches("Gen")) %>% t() %>% as.data.frame()
    gen <- gen[, 1]
    # Add transposed column to array temp
    array_temp <- array_temp %>% dplyr::select(-matches("Gen"))
    array_temp$Gen <- gen[which(gen>0)] %>% mean(na.rm=TRUE)
    # Fill in missing years (added with getContinous_installPeriodVector function) with appropriate values (NA for prior, and continued for post installation)
    #first_crop_year <- irrigEnergy_m2[i, ]$Year %>% as.numeric()
    #for(ii in c(1:nrow(array_temp))){
    #  array_temp[ii, ]$irrig <- ifelse(array_temp[ii, ]$Year<first_crop_year, NA, ifelse(is.na(array_temp[ii, ]$irrig), array_temp[ii-1, ]$irrig, array_temp[ii, ]$irrig))
    #  array_temp[ii, ]$ann_operation_MWh <- ifelse(array_temp[ii, ]$Year<first_crop_year, NA, ifelse(is.na(array_temp[ii, ]$ann_operation_MWh), array_temp[ii-1, ]$ann_operation_MWh, array_temp[ii, ]$ann_operation_MWh))
    #}
    # Remove rows for years before array was installed (irrelevant--wasted memory)
    array_temp <- array_temp[which(!is.na(array_temp$irrig)), ]
    array_temp$Index = index
    # Rbind with df
    df <- rbind(df, array_temp)
  }
  
  # Calculate overproducers vs under producers, split each array production for each year into Annual load vs overproduction
  df$overprod <- df$Gen - df$ann_operation_MWh
  df$perc_dev_annload <- (df$Gen / df$ann_operation_MWh) * 100 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Values not saved on export, but important to report
  
  # Calc some values for total over and under
  sum_overunder_irrig <- sum(df[which(df$irrig==1), ]$overprod) / MW_GW # USED IN SUPP
  sum_overunder_nonirrig <- sum(df[which(df$irrig==0), ]$overprod)  / MW_GW # USED IN SUPP
  # Calc totals for only those arrays that reach load
  sum_over_irrig <- sum(df[which(df$irrig==1 & df$perc_dev_annload>=100), ]$overprod) / MW_GW # USED IN SUPP
  sum_over_nonirrig <- sum(df[which(df$irrig==0 & df$perc_dev_annload>=100), ]$overprod)  / MW_GW # USED IN SUPP
  # Calc totals for only those arrays that DO NOT reach load
  sum_under_irrig <- sum(df[which(df$irrig==1 & df$perc_dev_annload<100), ]$overprod) / MW_GW # USED IN SUPP
  sum_under_nonirrig <- sum(df[which(df$irrig==0 & df$perc_dev_annload<100), ]$overprod)  / MW_GW # USED IN SUPP
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Values not saved on export, but important to report
  
  #df$perc_dev_annload <- (df$Gen - df$ann_operation_MWh) / df$ann_operation_MWh * 100 
  load_df <- data.frame(Year = integer(), Index = integer(), irrig = integer(), ann_operation_MWh = numeric(), Gen = numeric(), overprod = numeric(), perc_dev_annload = numeric(), value = numeric())
  for(i in c(1:nrow(df))){
    # Get row
    temp_df <- df[i, ]
    # Create an index of the rows you want with duplications & Use that index to duplicate
    temp_df <- temp_df[rep(1:nrow(temp_df), 3),]
    # Give new load column
    temp_df$load <- c("Annual Load", "Surplus", "Deficit")
    # Apply appropraite values
    temp_df$value <- ifelse(temp_df$load=="Annual Load" & temp_df$overprod<0, temp_df$Gen, 
                     ifelse(temp_df$load=="Annual Load" & temp_df$overprod>0, temp_df$ann_operation_MWh, 
                     ifelse(temp_df$load=="Surplus" & temp_df$overprod<0, 0, 
                     ifelse(temp_df$load=="Surplus" & temp_df$overprod>0, temp_df$overprod, 
                     ifelse(temp_df$load=="Deficit" & temp_df$overprod<0, temp_df$overprod, 
                     ifelse(temp_df$load=="Deficit" & temp_df$overprod>0, 0, NA))))))
    # Rbind with load df
    load_df <- rbind(load_df, temp_df)
  }
  
  # Export
  load_df$Year <- as.numeric(load_df$Year)
  return(load_df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate irrigation water offset if farm of average size is fallowed 
getResource_irrigWaterFarmFallow <- function(in_df){
  
  #~~~~~~~~~~~~~ Get crop relevant irrigation water depth (m3/m2)
  
  # Irrig water results
  irrigWater_result <- getResource_irrigWater(in_df)
  # Get rate of irrigWater offset (m3/m2) from result
  irrigWater_result[, 2:4] <- irrigWater_result[, 2:4] / in_df$dir_a
  
  #~~~~~~~~~~~~~ Get Farm Size
  
    # Get counties for each array
  counties <- counties(state, TRUE) %>% st_as_sf() %>% st_set_crs(st_crs(in_df))
  intersection <- st_intersection(st_centroid(in_df), counties)
  in_df$County <- intersection$NAME[match(in_df$Index, intersection$Index)]
  in_df$COUNTYFP <- counties$COUNTYFP[match(in_df$County, counties$NAME)]
  in_df$COUNTYFP <- gsub("(?<![0-9])0+", "", in_df$COUNTYFP, perl = TRUE) # Extracts non-zero led countyfp numerically from a string
  rm(counties)
  # Get county level average farm size
  state_county_agarea <- county_aglandarea[which(tolower(county_aglandarea$State)==tolower(state)), ]
  ag_farmsize <- state_county_agarea[which(state_county_agarea$Data.Item=="AG LAND - ACRES"), ]
  ag_farmsize$Value <- as.numeric(gsub(",","",ag_farmsize$Value))
  ag_farmsize$area <- ag_farmsize$Value * acres_m2
  ag_farmsize$farms <- state_county_agarea$Value[which(state_county_agarea$Data.Item=="AG LAND - NUMBER OF OPERATIONS")]
  ag_farmsize$farms <- as.numeric(gsub(",","", ag_farmsize$farms))
  ag_farmsize$avg_area <- ag_farmsize$area / ag_farmsize$farm
  ag_farmsize <- ag_farmsize[which(!is.na(ag_farmsize$avg_area)), ]
  # Save county average farm size to array based on county and survey year closest to installation -- farm_size here results in m2
  census_yrs = unique(ag_farmsize$Year) %>% as.numeric() 
  in_df$farm_size = median(ag_farmsize$avg_area[which(ag_farmsize$Year==max(census_yrs, na.rm = TRUE))], na.rm = TRUE) # save avearge farm size across state in latest census yearincase any missing data
  for(i in c(1:nrow(in_df))){
    county = in_df[i, ]$County %>% tolower()
    install_yr = in_df[i, ]$Year
    near_census_yr = census_yrs[which(abs(install_yr-census_yrs)==min(abs(install_yr-census_yrs)))] # Get closest census year to installation
    avg_farm_size = ag_farmsize$avg_area[which(as.numeric(ag_farmsize$Year)==near_census_yr & tolower(ag_farmsize$County)==county)]
    avg_farm_size = ifelse(length(avg_farm_size)==0, median(ag_farmsize$avg_area[which(tolower(ag_farmsize$County)==county)], na.rm = TRUE), avg_farm_size) # If nearest census year missing, choose most recent county census year
    in_df[i, ]$farm_size <- ifelse(length(avg_farm_size)==0, in_df[i, ]$farm_size, avg_farm_size) # If still missing, use state wide average
  }
  
  #~~~~~~~~~~~~~ Find offset water use for average farm size fallowing
  
  # Apply farm size to irrigation water use offset rather than direct solar area (returns to m3)
  irrigWater_result[, 2:4] <- irrigWater_result[, 2:4] * in_df$farm_size
  irrigWater_result$farm_size <- in_df$farm_size
  # Return dataframe
  return(irrigWater_result)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to clean out generated data and dataframes. Allows for regeneration of these data with new inputs 
cleanFEWLS = function(){
  # Check if irrigation depth dataframe exists
  if(file.exists("Data/Derived/mm_per_yearIrrig.csv")){
    file.remove("Data/Derived/mm_per_yearIrrig.csv") # If file has been generated, delete
  }
  
  # Check if crop revenue dataframe exists
  if(file.exists("Data/Derived/crop_revenue_kg.csv")){
    file.remove("Data/Derived/crop_revenue_kg.csv") # If file has been generated, delete
  }
  
  # Check if operational cost dataframe exists
  if(file.exists("Data/Derived/operational_cost_m2.csv")){
    file.remove("Data/Derived/operational_cost_m2.csv") # If file has been generated, delete
  }
  
  # Check if kcal/m2 dataframe exists
  if(file.exists("Data/Derived/kcal_per_m2.csv")){ 
    file.remove("Data/Derived/kcal_per_m2.csv") # If file has been generated, delete
  }
  
  # Check if irrigEnergyReq_adjusted dataframe exists
  if(file.exists("Data/Derived/irrigEnergyReq_adjusted.csv")){
    file.remove("Data/Derived/irrigEnergyReq_adjusted.csv") # If file has been generated, delete
  }
  print("FEWLS Memory Cleaned, Ready to Re-start.")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save important auxilary output data
getData_info = function(in_df){
  # --------------------------- #
  # Get dataset info and figure # 
  # --------------------------- #
  # Get lat long coords
  coords <- as.data.frame(st_coordinates(st_centroid(in_df$geometry)))
  in_df$Lat <- coords$Y
  in_df$Long <- coords$X
  
  # Get commercial and utility
  Commercial <- in_df[which(in_df$Capacity<capacity_threshold),]
  Utility <- in_df[which(in_df$Capacity>=capacity_threshold),]
  
  # Get mean lat for both scales -- NOT SAVED
  comm_lat <- mean(Commercial$Lat, na.rm=TRUE)
  util_lat <- mean(Utility$Lat, na.rm=TRUE)
  avg_dist <- (comm_lat - util_lat)*60*60 # Into arc seconds
  
  # Create dataframe to export
  df <- data.frame(Info = character(), Commercial = numeric(), Utility = numeric(), Unit = character())
  df <- add_row(df, Info = "Capacity", Commercial = round(sum(Commercial$Capacity), 2), Utility = round(sum(Utility$Capacity), 2), Unit = "MW")
  df <- add_row(df, Info = "Area", Commercial = round(sum(Commercial$dir_a/1000000), 2), Utility = round(sum(Utility$dir_a/1000000), 2), Unit = "km2")
  df <- add_row(df, Info = "Number", Commercial = as.integer(nrow(Commercial)), Utility = as.integer(nrow(Utility)), Unit = "arrays")
  df <- add_row(df, Info = "Previously Irrig", Commercial = as.integer(nrow(Commercial[which(Commercial$irrig==1), ])), Utility = as.integer(nrow(Utility[which(Utility$irrig==1), ])), Unit = "arrays")
  write.csv(df, paste("Outputs/Dataframes/Dataset_Information_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""), row.names=FALSE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Prep array tech size and number plots
getPlot_Array_Tech_Dist = function(in_df){
  # Change single double axis array to single axis for simplicity
  in_df$Class <- ifelse(in_df$Class=="double_axis", "single_axis", in_df$Class)
  # Remove geometry
  in_df$geometry = NULL
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  =  in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  # Make plot function
  getPlot_Array_Tech_Dist_temp = function(plot_df){
    #___________________________ Prep
    # Technology grouping
    groupPV_class_df <- plot_df %>%
      group_by(Year, Class) %>% 
      summarise(num = n(),
                Average_size = mean(Capacity))
    groupPV_class_df$Year = as.numeric(groupPV_class_df$Year)
    # Capacity group
    groupPV_df <- plot_df %>%
      group_by(Year) %>% 
      summarise(num = n(),
                dir_a = sum(Capacity), 
                Average_size = mean(Capacity))
    groupPV_df$Year = as.numeric(groupPV_df$Year)
    # Pull from class split
    single <- groupPV_class_df[which(groupPV_class_df$Class=='single_axis'), ]
    groupPV_df$Avgsize_single <- 0
    groupPV_df$Avgsize_single <- single$Average_size[match(groupPV_df$Year, single$Year)]
    fixed <- groupPV_class_df[which(groupPV_class_df$Class=='fixed_axis'), ] 
    groupPV_df$Avgsize_fixed <- 0
    groupPV_df$Avgsize_fixed <- fixed$Average_size[match(groupPV_df$Year, fixed$Year)] 
    groupPV_class_df$Class[which(groupPV_class_df$Class=="fixed_axis")] <- "Fixed Axis"
    groupPV_class_df$Class[which(groupPV_class_df$Class=="single_axis")] <- "Single Axis"
    # Get graph params
    coeff = (max(groupPV_df$num, na.rm = TRUE) / max(groupPV_df$Avgsize_single, na.rm = TRUE)) %>% plyr::round_any(10, f=floor)
    coeff = ifelse(coeff==0, (max(groupPV_df$num, na.rm = TRUE) / max(groupPV_df$Avgsize_single, na.rm = TRUE)) %>% plyr::round_any(1, f=floor), coeff)
    coeff = ifelse(coeff==0, 0.8, coeff)
    max_y = plyr::round_any(max(groupPV_df$num), 10, f=ceiling)
    graph_name = deparse(substitute(plot_df))
    #___________________________ Plot  
    tech_size_plot <- ggplot(groupPV_df, aes(Year, num)) +
      geom_bar(data=groupPV_class_df, stat ="identity", aes(fill=Class), color = 'black') +
      geom_line(aes(y = Avgsize_fixed*coeff), size=0.5, color="black", linetype=5) + # alpha = 1/3
      geom_line(aes(y = Avgsize_single*coeff), size=0.5, color="black", linetype=5) + # alpha = 1/3
      geom_point(aes(y = Avgsize_fixed*coeff), size=4, color="#999999", shape=19) +
      geom_point(aes(y = Avgsize_fixed*coeff), size=4, color="black", shape=21) +
      geom_point(aes(y = Avgsize_single*coeff), size=4, color="#E69F00", shape=19) +
      geom_point(aes(y = Avgsize_single*coeff), size=4, color="black", shape=21) +
      scale_fill_manual(values=c("#999999", "#E69F00")) +
      scale_colour_viridis_c(option = "inferno")  +
      scale_y_continuous(name = "Installed Arrays", sec.axis = sec_axis(~./coeff, name="Array Capacity (MW)"), limits = c(0, max_y), expand = c(0,0)) +
      scale_x_continuous(breaks=seq(start_install_year, end_install_year, 1), limits = c(start_install_year-0.5, end_install_year+0.5), expand = c(0,0)) + 
      theme_minimal() + 
      theme(legend.position = "none", plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), 
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.border = element_rect(colour="black", fill=NA, size = 1), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      annotate("text", x = start_install_year-0.5, y=max_y*0.95, label = graph_name, hjust=0, fontface=2)
  }
  # Get commercial and utility
  Commercial <- in_df[which(in_df$Capacity<capacity_threshold),]
  Utility <- in_df[which(in_df$Capacity>=capacity_threshold),]
  # Get plots
  commercial_plot <- getPlot_Array_Tech_Dist_temp(Commercial) + theme(axis.text.x = element_blank())
  utility_plot <- getPlot_Array_Tech_Dist_temp(Utility)
  fig <- ggarrange(commercial_plot, utility_plot, nrow = 2, ncol = 1, align = c("v"), heights = c(9,10)) 
  fig <- annotate_figure(fig, bottom = text_grob("Installation Year", size = 12), left = text_grob("Number of Installations", size = 12, rot = 90),
                         right = text_grob("Average Array Capacity - - - (MW)", size = 12, rot = 90))
  ggsave(path = 'Outputs/Figures', filename = paste("tech_size_dist_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".pdf", sep=""), width = 4.75, height = 5, units = "in", useDingbats=FALSE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Auxilary Graph if desired. Dispalys surplus, deficit, and average on farm budgets by calculating on farm energy requirements
getFarmElec_Budget = function(in_df){
  # ----------------------------------------- #
  # Get on-farm reqs, and surplis and deficit # 
  # ----------------------------------------- #
  # A double-check to ensure we are only looking at commercial scale arrays
  in_df <- in_df[which(in_df$Capacity<capacity_threshold), ]
  # Run annual farm operational requirement budget function
  load_df <- getResource_annFarmOperationReq(in_df)
  irrig <- load_df[which(load_df$irrig==1), ]
  nonirrig <- load_df[which(load_df$irrig==0), ]
  
  # Get start and end install year
  start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
  end_install_year  =  in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Plot with proportion of annual load met
  
  # Plot vars
  margin_plot <- unit(c(0.5,0.5,0,0),"cm")
  alpha = 0.3
  lwd = 1.5
  
  # Annual load budget (if net-electricity is desired as opposed to proportion of load met -- helps know on farm level how much each farm over or under produces on average)
  #irriggrouped <- irrig %>% group_by(Year) %>% summarise(median = as.numeric(summary(overprod)[3])/MW_GW, firstq = as.numeric(summary(overprod)[2])/MW_GW, thirdq = as.numeric(summary(overprod)[5])/MW_GW) %>% as.data.frame()
  #nonirriggrouped <- nonirrig %>% group_by(Year) %>% summarise(median = as.numeric(summary(overprod)[3])/MW_GW, firstq = as.numeric(summary(overprod)[2])/MW_GW, thirdq = as.numeric(summary(overprod)[5])/MW_GW)  %>% as.data.frame()
  # Proportion of Annual Load Met (%)
  irriggrouped <- irrig %>% group_by(Index, Year) %>% 
    summarise(perc_dev_annload = mean(perc_dev_annload, na.rm=TRUE)) %>% 
    group_by(Year) %>% 
    summarise(median = as.numeric(summary(perc_dev_annload)[3]), firstq = as.numeric(summary(perc_dev_annload)[2]), thirdq = as.numeric(summary(perc_dev_annload)[5])) %>% 
    as.data.frame()
  nonirriggrouped <- nonirrig %>% group_by(Index, Year) %>% 
    summarise(perc_dev_annload = mean(perc_dev_annload, na.rm=TRUE)) %>% 
    group_by(Year) %>% 
    summarise(median = as.numeric(summary(perc_dev_annload)[3]), firstq = as.numeric(summary(perc_dev_annload)[2]), thirdq = as.numeric(summary(perc_dev_annload)[5])) %>% 
    as.data.frame()
  
  # Calculate y limits -- rounds to nearest 0.5 GW above and below bounds
  max_y = max(max(irriggrouped$thirdq, na.rm = TRUE) , max(nonirriggrouped$thirdq, na.rm = TRUE), na.rm = TRUE) %>% plyr::round_any(50, f=ceiling) # %>% plyr::round_any(0.5, f=ceiling) -- if doing net elec (annual load budget)
  min_y = min(min(irriggrouped$firstq, na.rm = TRUE) , min(nonirriggrouped$firstq, na.rm = TRUE), na.rm = TRUE) %>% plyr::round_any(50, f=floor) # %>% plyr::round_any(0.5, f=ceiling) -- if doing net elec (annual load budget)
  break_rate = ((max_y - min_y) / 100) %>% plyr::round_any(200, f=ceiling) # ((max_y - min_y) / 6) %>% plyr::round_any(0.5, f=ceiling) # Result is 6 breaks -- if doing net elec (annual load budget)
  # Plot
  genMinusReq_plot <- ggplot(irriggrouped, aes(x=Year)) +
    geom_ribbon(data=irriggrouped, aes(ymin = firstq, ymax = thirdq), fill = "blue", alpha = alpha) +
    geom_ribbon(data=nonirriggrouped, aes(ymin = firstq, ymax = thirdq), fill = "goldenrod4", alpha = alpha) +
    geom_line(data=irriggrouped, aes(y=median), col = "blue", size = lwd) +
    geom_line(data=nonirriggrouped, aes(y=median), col = "goldenrod4", size = lwd) +
    ylab("Local Avg-Proportion of Annual Load (%)") + 
    geom_hline(yintercept=100, linetype = "dashed", color= "black", size = 1) +
    scale_x_continuous(breaks=seq(start_install_year, end_install_year, 1), limits = c(start_install_year-0.5, end_install_year+0.5), expand = c(0,0)) + 
    scale_y_continuous(breaks=seq(min_y, max_y, break_rate), limits = c(min_y, max_y), expand=c(0,0)) +
    theme_minimal() + 
    theme(legend.position = "none", plot.title = element_blank(), panel.border = element_rect(colour="black", fill=NA, size = 1), axis.text.x = element_blank(), 
          axis.title.x = element_blank(), plot.margin = margin_plot, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Plot with medians and quartiles (surplus and deficit)
  
  # Set plot vars
  pnt_size = 2.5
  lwd = 0.75
  # Create line plot with total surplus and deficit
  # Surplus DF
  surplusIrrig <- irrig[irrig$load=="Surplus", ] %>% group_by(Year) %>% summarize(over = sum(value)/MW_GW)
  surplusNonir <- nonirrig[nonirrig$load=="Surplus", ] %>% group_by(Year) %>% summarize(over = sum(value)/MW_GW)
  # Deficit DF
  deficitIrrig <- irrig[irrig$load=="Deficit", ] %>% group_by(Year) %>% summarize(under = sum(value)/MW_GW)
  deficitNonir <- nonirrig[nonirrig$load=="Deficit", ] %>% group_by(Year) %>% summarize(under = sum(value)/MW_GW)
  
  # Calculate y limits -- rounds to nearest 0.5 GW above and below bounds
  max_y = max(max(surplusIrrig$over, na.rm = TRUE) , max(surplusNonir$over, na.rm = TRUE), na.rm = TRUE) %>% plyr::round_any(10, f=ceiling)
  min_y = min(min(deficitIrrig, na.rm = TRUE) , min(surplusNonir$under, na.rm = TRUE), na.rm = TRUE) %>% plyr::round_any(10, f=floor)
  break_rate = ((max_y - min_y) / 10) %>% plyr::round_any(10, f=ceiling) # Result is 6 breaks
  
  # Plot
  SurplusDeficit_plot <- ggplot() +
    # Zero line
    geom_hline(yintercept = 0, size=0.75) +
    # Surplus
    geom_line(data = surplusIrrig, aes(x=Year, y=over), size=lwd, color="black", linetype=5) +
    geom_line(data = surplusNonir, aes(x=Year, y=over), size=lwd, color="black", linetype=5) +
    geom_point(data = surplusIrrig, aes(x=Year, y=over), size=pnt_size, color="blue", shape=19, alpha = 1) +
    geom_point(data = surplusIrrig, aes(x=Year, y=over), size=pnt_size, color="black", shape=21, alpha = 1) +
    geom_point(data = surplusNonir, aes(x=Year, y=over), size=pnt_size, color="darkgoldenrod4", shape=19, alpha = 1) +
    geom_point(data = surplusNonir, aes(x=Year, y=over), size=pnt_size, color="black", shape=21, alpha = 1) +
    # Deficit
    geom_line(data = deficitIrrig, aes(x=Year, y=under), size=lwd, color="black", linetype=5) +
    geom_line(data = deficitNonir, aes(x=Year, y=under), size=lwd, color="black", linetype=5) +
    geom_point(data = deficitIrrig, aes(x=Year, y=under), size=pnt_size, color="blue", shape=19, alpha = 0.5) +
    geom_point(data = deficitIrrig, aes(x=Year, y=under), size=pnt_size, color="black", shape=21, alpha = 1) +
    geom_point(data = deficitNonir, aes(x=Year, y=under), size=pnt_size, color="darkgoldenrod4", shape=19, alpha = 0.5) +
    geom_point(data = deficitNonir, aes(x=Year, y=under), size=pnt_size, color="black", shape=21, alpha = 1) +
    # Other graphing variables
    ylab("Regional Surplus and Deficit (GWh)") + 
    xlab("Year") +
    scale_y_continuous(breaks=seq(min_y, max_y, break_rate), limits = c(min_y-(break_rate/4), max_y*1.05), expand=c(0,0)) + # Min y should always be zero
    scale_x_continuous(breaks=seq(start_install_year, end_install_year, 2), limits = c(start_install_year-0.5, end_install_year+0.5), expand = c(0,0)) + 
    scale_fill_manual(values=c("black", "green")) +
    theme_minimal() + 
    theme(legend.position = "none", plot.title = element_blank(), panel.border = element_rect(colour="black", fill=NA, size = 1), plot.margin = margin_plot, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  # Plot together
  fig <- ggarrange(genMinusReq_plot, SurplusDeficit_plot, nrow = 2, ncol = 1, align = c("hv"), heights = c(5,5))
  ggsave(path = 'Outputs/Figures', filename = paste("Surplus_Deficits_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".pdf", sep = ""), width = 3.5, height = 6, units = "in", useDingbats=FALSE)
  return(fig)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get compile time
end_time <- Sys.time()
tot_time = (end_time - start_time) %>% format(format = "%H:%M:%S")# %>% round(2)
print(paste("Model compiled in ", tot_time, sep=""))

#_______________________________________________________________________________________________________________________________________ END FEWLS TOOL MODEL BUILD




