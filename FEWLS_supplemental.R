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
This is all supplemental code used in the publication *Agrisolar Co-Location: Balancing
Food and Energy Production Tradeoffs with Water Security*. Numerous supplemetnal
operations are perfomed here to either help set-up initial solar dataset (such as merging
Stid et al., 2022 solar PV data with Kruitwagen et al., 2021 solar PV data), generate values
for the publication itself derived from FEWLS_model results, and various other non-model-essential
intermediate products or depreicated products and code. If something is lost, its probabaly here, 
and then some. 

These are broad descritpions for processes peformed here: 
- A few quantitive summations to report in Agrisolar text
- Get average electricity rate (sensivity analysis, and utility scale water value)
- Get average rate of inflation for CPI and PPI
- Code to generate figure of individual array payback times
- Code to check if any crop types are missing crop dollar data (would need to updated FEWLS_cropRevenue.R)
- Code to generate Sankey diagram for fallowed proximal land to solar PV and associated data (and irrigation water use offset quantification)
- Code to generate figure for conference presentation on how we extrapolate water use to years outside of survey years
- Code for pre- and post-processing (from AcrPro) for Figure 2 in Agrisolar manuscript (study area)
- Code for looking at historical cultivated vs fallow land area in CCV between 2008 and 2021
- Code for acquiring associated utility name for each array
- Preliminary code for generating separate dataframes for non_irrig_annual_load_overpoduction and irrig_annual_load_overpoduction (Now performed in FEWLS_FarmElecBudgetExport.R
- Code to pre-process and stitch Stid et al. (2022) and Kruitwagen et al. (2021) Togegher
- Original and re-processed crop rotation selection for agriculturally co-located arrays between Kruitwagen and Stid datasets
- Check which arrays removed from Stid et al., 2022 by new ag-co-located definition (discussion in supplementary)
- Merging electricity rate model dataframe with in_solar dataframe (both Stid and Kruitwagen)
"

## ----------------------------------------------------------- ##
##                                                             ##
##            Supplemental Work, Data, and Plots               ##
##                                                             ##
## ----------------------------------------------------------- ##

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ A few paper numbers

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## MAIN TEXT NUMBERS --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# A few important paper numbers
tot_area = sum(in_solar_df$dir_a/1e6) 
tot_cap = sum(in_solar_df$Capacity) 
numKruitwagen =  in_solar_df[which(in_solar_df$Source=="Kruitwagen"), ] %>% nrow()
numStid =  in_solar_df[which(in_solar_df$Source=="Stid"), ] %>% nrow()
rangeCapacity = range(in_solar_df$Capacity)
percActiveAgland = tot_area / 38335 * 100 # number from CDL zonal stat of cultivated land (Data\Derived\GEE_ActiveAg_FallowLand_df)


# Get number of arrays associated with rotation
in_df <- in_solar_df
in_df$geometry=NULL
rotation_crop <- getCrop_rotation(in_df)
rotation <- rotation_crop %>% group_by(Index) %>% summarize(num = n(), Capacity = first(Capacity), dir_a = sum(dir_a))
rotation_num <- rotation[which(rotation$num>1), ] %>% nrow()
rotation_cap <- rotation$Capacity[which(rotation$num>1)] %>% sum()
rotation_dir_a <- rotation$dir_a[which(rotation$num>1)] %>% sum() / 1000000
# Get orchard and grape area
orchards <- Nass_Classifications$Crop[which(Nass_Classifications$Irrig_Crop=="Orchards")] %>% tolower()
berries <- Nass_Classifications$Crop[which(Nass_Classifications$Irrig_Crop=="Berry totals")] %>% tolower()
rotation_OrchardArea <- rotation_crop$dir_a[which(tolower(rotation_crop$Crp_n_1) %in% orchards)] %>% sum() / m2_ha

# Food, Energy, Water Real Life Comparisons
commFEW <- read.csv(paste("Outputs/Dataframes/ResourceModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
utilFEW <- read.csv(paste("Outputs/Dataframes/ResourceModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
# FOOD
commkcal <- commFEW %>% group_by(Year) %>% summarise(result = sum(kcal))
utilkcal <- utilFEW %>% group_by(Year) %>% summarise(result = sum(kcal))
kcal_yr = 2000 * 365 # 2500 Calories/day to Calories per year
people_yrLSfeed = (tail(commkcal$result, 1) + tail(utilkcal$result, 1)) / kcal_yr / system_lifespan
totKcalMWyr = (tail(commkcal$result, 1) + tail(utilkcal$result, 1)) / sum(in_solar_df$Capacity) / system_lifespan
# ENERGY
commEN <- commFEW %>% group_by(Year) %>% summarise(result = sum(energy_base))
utilEN <- utilFEW %>% group_by(Year) %>% summarise(result = sum(energy_base))
homesEN_yr = 10.6 / 1000 # MWh/yr to GWh/yr
homesPowered_yr = (tail(commEN$result, 1) + tail(utilEN$result, 1)) / homesEN_yr / system_lifespan
# ENERGY IRRIG
commENI <- commFEW %>% group_by(Year) %>% summarise(result = sum(irrEngy_fg))
utilENI <- utilFEW %>% group_by(Year) %>% summarise(result = sum(irrEngy_fg))
homesPowerd_yrIrrig = (tail(commENI$result, 1) + tail(utilENI$result, 1)) / homesEN_yr / system_lifespan
# WATER
commWA <- commFEW %>% group_by(Year) %>% summarise(result = sum(irrig_fg))
utilWA <- utilFEW %>% group_by(Year) %>% summarise(result = sum(irrig_fg))
oz_toliters_human_day = 80 / 33.814 # 80oz/day / 33.814 oz per liter
m3_yr_human = oz_toliters_human_day / 1000 * 365 #  80oz/day / liters/m3 * day/year = m3/yr
meter_orcharIrrig = (822.96 + 762)/2 / 1000 # Average orchard irrigDept Surveys (mm) / 1000mm/m
orchardAcres_yr = (tail(commWA$result, 1) + tail(utilWA$result, 1)) / meter_orcharIrrig / system_lifespan / m2_ha # m2 per acre
drinkWaterppl = ((tail(commWA$result, 1) + tail(utilWA$result, 1)) / m3_yr_human / system_lifespan) %>% signif(3)


# Food, Energy, Water perCrop comparisons
commFEW <- read.csv(paste("Outputs/Dataframes/ResourceModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
utilFEW <- read.csv(paste("Outputs/Dataframes/ResourceModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
# Add fars crops to both
nass = Nass_Classifications
nass$Crop = nass$Crop %>% tolower()
# Correct nass crops to modified FARS crops -- pulled from sankey diagram generation below
crp_grp1 = "wheat|sorghum|grain|rice"
crp_grp2 = "vegetable|tomatoes|potatoes|beans|corn"
crp_grp3 = "hay|pastureland"
crp_nme1 = "grain"
crp_nme2 = "vegetables"
crp_nme3 = "hay/pasture"
nass$sankey_crop <- tolower(nass$Irrig_Crop)
# Group further to simplify
nass$sankey_crop <- ifelse(grepl(crp_grp1,nass$sankey_crop), crp_nme1, nass$sankey_crop)
nass$sankey_crop <- ifelse(grepl(crp_grp2,nass$sankey_crop), crp_nme2, nass$sankey_crop)
nass$sankey_crop <- ifelse(grepl(crp_grp3,nass$sankey_crop), crp_nme3, nass$sankey_crop)
# Add custom for fallow idle and all other non-ag
nass$sankey_crop <- ifelse(nass$Value %in% c(61, 65), "fallow/idle", nass$sankey_crop)
nass$sankey_crop <- ifelse(is.na(nass$sankey_crop), "non-ag", nass$sankey_crop)
# Combine with outputs
commFEW$Crop <- commFEW$Crop  %>% tolower()
utilFEW$Crop <- utilFEW$Crop  %>% tolower()
# Really these are modified FARS crops
commFEW$FARS_crop <- nass$sankey_crop[match(commFEW$Crop, nass$Crop)] %>% tolower()
utilFEW$FARS_crop <- nass$sankey_crop[match(utilFEW$Crop, nass$Crop)] %>% tolower()
# Group and sum - USED IN PAPER
commFEW_crop <- commFEW %>% dplyr::select(!matches("Crop|Year|Index|Yr_inst", ignore.case = FALSE)) %>% group_by(FARS_crop) %>% summarise_all(sum) %>% mutate(across(where(is.numeric), ~ ./sum(.)))
utilFEW_crop <- utilFEW %>% dplyr::select(!matches("Crop|Year|Index|Yr_inst", ignore.case = FALSE)) %>% group_by(FARS_crop) %>% summarise_all(sum) %>%  mutate(across(where(is.numeric), ~ ./sum(.)))
# Get dataset total for kcal and irrig forgone - USED IN PAPER
FEW_crop <- rbind( commFEW %>% dplyr::select(matches("kcal|irrig_fg|FARS_crop")) , utilFEW %>% dplyr::select(matches("kcal|irrig_fg|FARS_crop")))
FEW_crop <- FEW_crop %>% group_by(FARS_crop) %>% summarise_all(sum) %>%  mutate(across(where(is.numeric), ~ ./sum(.)))


# Percent irrigated by scale
irrigArea = sum(in_solar_df$dir_a[which(in_solar_df$irrig==1)]/1000000)
datainfo = read.csv(paste("Outputs/Dataframes/Dataset_Information_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
percCommIrrig = datainfo$Commercial[which(datainfo$Info=="Previously Irrig")] / datainfo$Commercial[which(datainfo$Info=="Number")] * 100
percUtilIrrig = datainfo$Utility[which(datainfo$Info=="Previously Irrig")] / datainfo$Utility[which(datainfo$Info=="Number")] * 100
perIrrig_all = (datainfo$Commercial[which(datainfo$Info=="Previously Irrig")] + datainfo$Utility[which(datainfo$Info=="Previously Irrig")]) / nrow(in_solar_df) * 100


# Ratio of irrig water saved to O&M water used
commWA <- commFEW %>% group_by(Year) %>% summarise(result = sum(irrig_fg))
commOM <- commFEW %>% group_by(Year) %>% summarise(result = sum(oandmIr_base))
waterSaveIrrOandM = (tail(commWA$result, 1) / tail(commOM$result, 1)) 


# BIG NUMBER: Final budget / Food revenue loss
commecon <- read.csv(paste("Outputs/Dataframes/EconModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
commBudgSuperiority = commecon %>% group_by(Year) %>% summarise(budget = sum(Tot_Budg_base), food = sum(food_base))
finalEconComm_ratio = tail(commBudgSuperiority$budget, 1) / tail(commBudgSuperiority$food, 1)


# Get economic payback times (commercial, utility, and with half operational costs)
commecon <- read.csv(paste("Outputs/Dataframes/EconModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
utilecon <- read.csv(paste("Outputs/Dataframes/EconModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
# Get payback time vector with base methods
getPaybackTimes = function(df){
  df_temp <- df %>% dplyr::select(matches("Index|Year|Yr_inst|Tot_"))
  # For each index find base, best, and worst payback time and save out
  payback_df <- data.frame(base = as.numeric(), best = as.numeric(), worst = as.numeric(), Yr_inst = as.integer(), Index = as.integer())
  for(indx in unique(df_temp$Index)){
    payback_temp <- df_temp[which(df_temp$Index == indx), ]
    yr_inst = first(payback_temp$Yr_inst, na_rm = TRUE)
    payback_temp <- payback_temp[which(payback_temp$Year>=yr_inst), ]
    # Gets dataframe where row are only those with positive total econ values
    payback_base <- payback_temp %>% group_by(Year) %>% mutate(Inx = first(which(Tot_Budg_base > 0))) %>% filter(Inx == row_number()) %>% dplyr::select(-Inx) 
    payback_base <- as.numeric(min(payback_base$Year) - yr_inst)
    payback_best <- payback_temp %>% group_by(Year) %>% mutate(Inx = first(which(Tot_Budg_max > 0))) %>% filter(Inx == row_number()) %>% dplyr::select(-Inx) 
    payback_best <- as.numeric(min(payback_best$Year) - yr_inst)
    payback_worst <- payback_temp %>% group_by(Year) %>% mutate(Inx = first(which(Tot_Budg_min > 0))) %>% filter(Inx == row_number()) %>% dplyr::select(-Inx) 
    payback_worst <- as.numeric(min(payback_worst$Year) - yr_inst)
    payback_df <- rbind(payback_df, data.frame(base = payback_base, best = payback_best, worst = payback_worst, Yr_inst = yr_inst, Index = indx))
  }
  return(payback_df)
}
# Run for comm and util under standard conditions (regular operational costs)
comm_payback = getPaybackTimes(commecon)
util_payback = getPaybackTimes(utilecon)
# Find number of installations never reading payback
comm_payback_numNA <- comm_payback[is.infinite(comm_payback$base),] %>% nrow() %>% as.numeric()
util_payback_numNA <- util_payback[is.infinite(util_payback$base),] %>% nrow() %>% as.numeric()
# Get numeber of and remove installations with no data (inf or -inf)
comm_payback_numNA <- comm_payback[is.infinite(rowSums(comm_payback)),] %>% nrow() %>% as.numeric()
util_payback_numNA <- util_payback[is.infinite(rowSums(util_payback)),] %>% nrow() %>% as.numeric()
comm_payback <- comm_payback[!is.infinite(rowSums(comm_payback)),]
util_payback <- util_payback[!is.infinite(rowSums(util_payback)),]
# Get average payback times
comm_paybackTIME = colMeans(comm_payback[, c(1:3)])
util_paybackTIME = colMeans(util_payback[, c(1:3)])


# Run for comm and util under special conditions (half operational costs)
# First, remove half of all operational costs from each total budget column
reduce_factor = 0.5
commecon$Tot_Budg_base = commecon$Tot_Budg_base - (commecon$operation_base*reduce_factor)
commecon$Tot_Budg_max = commecon$Tot_Budg_max - (commecon$operation_max*reduce_factor)
commecon$Tot_Budg_min = commecon$Tot_Budg_min - (commecon$operation_min*reduce_factor)
# Rerun same analysis as above
comm_payback = getPaybackTimes(commecon)
util_payback = getPaybackTimes(utilecon)
# Get numeber of and remove installations with no data (inf or -inf)
comm_payback_numNA <- comm_payback[is.infinite(rowSums(comm_payback)),] %>% nrow() %>% as.numeric()
util_payback_numNA <- util_payback[is.infinite(rowSums(util_payback)),] %>% nrow() %>% as.numeric()
comm_payback <- comm_payback[!is.infinite(rowSums(comm_payback)),]
util_payback <- util_payback[!is.infinite(rowSums(util_payback)),]
# Get average payback times
comm_specialpaybackTIME = colMeans(comm_payback[, c(1:3)])
util_specialpaybackTIME = colMeans(util_payback[, c(1:3)])


# Proportion of California Energy Production 
# Food, Energy, Water Real Life Comparisons
gen2018 <- in_solar_df %>% dplyr::select(matches("Capacity|Gen_2018"))
gen2018$geometry = NULL
comm2018 = gen2018[which(gen2018$Capacity<capacity_threshold), ] %>% sum() / MW_GW
util2018 = gen2018[which(gen2018$Capacity>=capacity_threshold), ] %>% sum() / MW_GW
CA_2018totgen = 194842 # GWh
CA_2018solgen = 27265 # GWh
percCAtotGen = (comm2018 + util2018) / CA_2018totgen * 100
percCAsolGen = (comm2018 + util2018) / CA_2018solgen * 100


# Fallow/idle numbers below


# Proportion of dataset post SGMA
prop_afSGMA_comm <- in_solar_df[which(in_solar_df$Capacity<capacity_threshold & in_solar_df$Year > 2014), ] %>% nrow() / nrow(in_solar_df[which(in_solar_df$Capacity<capacity_threshold),])
prop_afSGMA_util <- in_solar_df[which(in_solar_df$Capacity>=capacity_threshold & in_solar_df$Year > 2014), ] %>% nrow() / nrow(in_solar_df[which(in_solar_df$Capacity>=capacity_threshold),])


# Value of irrigated land in CA per ha
irrigVal_ha = 14300 / acres_m2 * m2_ha


# Land Lease rates to USD ha-1 (from inputs, rates are in m2)
baseline_Lease = landlease_per_m2_avg * m2_ha
best_Lease = landlease_per_m2_max * m2_ha
worst_Lease = landlease_per_m2_min * m2_ha
otherWorst_Lease = 400 / acres_m2 * m2_ha
otherBest_Lease = 1200 / acres_m2 * m2_ha


# $/ac-ft to $/thousand-m3
USD_acft_thous_m3 = 40 / acft_m3 * 1000


# Range of NEM arrays used in results
nem_temp <- in_solar_df
nem_temp$geometry = NULL
nem_temp <- nem_temp[which(nem_temp$Capacity<capacity_threshold), ]
rangeNEMcap <- nem_temp %>% group_by(Class) %>% summarise(area = max(dir_a)/m2_ha)



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENATL NUMBERS -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Load percent met for commercial-scale farms (ALL)
loaddf = read.csv("Data/Derived/SolarFEWE_LoadDf_annual.csv")
loaddf = loaddf %>% group_by(Index) %>% summarise(percLoad = (mean(Gen) / mean(ann_operation_MWh) * 100), Year = mean(Year)) # means are irrelevant, all repeated numbers from how df was exported
loaddf_count = (loaddf[which(loaddf$percLoad>=100), ] %>% nrow()) / nrow(loaddf) * 100 # percentage of arrays meetin estimated load
# Load percent met for commercial-scale farms (Irrigated)
commIndices = in_solar_df$Index[which(in_solar_df$Capacity<capacity_threshold)]
# NON-Irrig
loaddf = read.csv("Data/Derived/SolarFEWE_LoadDf_annual.csv")
loaddf = loaddf[which(loaddf$irrig==0 & loaddf$Index %in% commIndices), ]
loaddf = loaddf %>% group_by(Index) %>% summarise(percLoad = (mean(Gen) / mean(ann_operation_MWh) * 100), Year = mean(Year)) # means are irrelevant, all repeated numbers from how df was exported
loaddf_med = median(loaddf$percLoad, na.rm=TRUE)
loaddf_count = (loaddf[which(loaddf$percLoad>=100), ] %>% nrow()) / nrow(loaddf) * 100
loaddf = loaddf %>% group_by(Year) %>% summarise(percLoad = median(percLoad))
# IRRIG
loaddf = read.csv("Data/Derived/SolarFEWE_LoadDf_annual.csv")
loaddf = loaddf[which(loaddf$irrig==1 & loaddf$Index %in% commIndices), ]
loaddf = loaddf %>% group_by(Index) %>% summarise(percLoad = (mean(Gen) / mean(ann_operation_MWh) * 100), Year = mean(Year)) # means are irrelevant, all repeated numbers from how df was exported
loaddf_med = median(loaddf$percLoad)
loaddf_count = (loaddf[which(loaddf$percLoad>=100), ] %>% nrow()) / nrow(loaddf) * 100
loaddf = loaddf %>% group_by(Year) %>% summarise(percLoad = median(percLoad))

# Actual GWh of surplus and deficit
load_df = read.csv("Data/Derived/SolarFEWE_LoadDf_annual.csv")
irrig <- load_df[which(load_df$irrig==1), ]
nonirrig <- load_df[which(load_df$irrig==0), ]
# Get start and end install year
start_install_year = in_solar_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
end_install_year  =  in_solar_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
# Get sums and medians for paper
total_irrig_surplus <- sum(irrig[irrig$load=="Surplus", ]$value) / MW_GW
total_nonir_surplus <- sum(nonirrig[nonirrig$load=="Surplus", ]$value) / MW_GW
total_irrig_deficit <- sum(irrig[irrig$load=="Deficit", ]$value) / MW_GW
total_nonir_deficit <- sum(nonirrig[nonirrig$load=="Deficit", ]$value) / MW_GW
median_percentOfLoad_irrig <- median(irrig$perc_dev_annload)
median_percentOfLoad_nonir <- median(nonirrig$perc_dev_annload)
# Year specific
#yr = 2018
#year_percentOfLoad_irrig <- median(irrig$perc_dev_annload[which(irrig$Year==yr)])
#year_percentOfLoad_nonir <- median(nonirrig$perc_dev_annload[which(nonirrig$Year==yr)])
#year_percOfLoad_irrig_group <- irrig %>% group_by(Year) %>% summarise(perc_dev_annload = median(perc_dev_annload))
#year_percOfLoad_nonir_group <- nonirrig %>% group_by(Year) %>% summarise(perc_dev_annload = median(perc_dev_annload))

# Some important values (mean loads for irrig and non irrig) in FEWLS_model.R in getResource_annFarmOperationReq() function



# Irrigation energy requiresment reporting
# Get df and get to kWh
irrigEn <- read.csv("Data/Derived/irrigEnergyReq_adjusted.csv")
irrigEn$GWh_m3 <- irrigEn$GWh_m3 * kWh_GWh * acft_m3
rangeIrrigEn <- range(irrigEn$GWh_m3)
meanIrrigEn <- mean(irrigEn$GWh_m3)



# Get total diff between mono- and muli-chrystalline silicon generation scenarios
commFEW <- read.csv(paste("Outputs/Dataframes/ResourceModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
utilFEW <- read.csv(paste("Outputs/Dataframes/ResourceModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
commEN <- commFEW %>% group_by(Year) %>% summarise(result_max = sum(energy_max), result_min = sum(energy_min))
utilEN <- utilFEW %>% group_by(Year) %>% summarise(result_max = sum(energy_max), result_min = sum(energy_min))
TWh_diff_scenario <- ( (tail(commEN$result_max, 1) + tail(utilEN$result_max, 1)) - (tail(commEN$result_min, 1) + tail(utilEN$result_min, 1)) ) / GW_TW
totKcalMWyr = (tail(commkcal$result, 1) + tail(utilkcal$result, 1)) / sum(in_solar_df$Capacity) / system_lifespan



# Predominant crops fallowed by solar -- Below in fallowed code



# Supplementary Discussion 3:  Interconnection fees
# Get interconneciotn fees
in_df <- in_solar_df[which(in_solar_df$Capacity<capacity_threshold), ]
in_df$geometry = NULL
yearsPV = in_df %>% group_by(Year) %>% summarise(num = n())
# NEM Interconnection Fee (one-time) -- $0 for NEM1
intconn_NEM2_high = 145 
intconn_NEM2_low = 75 
intconn_NEM2_base = 132
# Set interconnection variable
NEM1_period = c(0:2016)
NEM2_period = c(2017:2100)
intconn_num = in_df[which(in_df$Year %in% NEM2_period), ] %>% nrow() 
intconn_feeWorst = intconn_num * intconn_NEM2_high
# Everthing below here was kept from when this fee was included in the model
#in_df$interconnection_fee_base <- ifelse(in_df$Year %in% NEM1_period, 0, intconn_NEM2_base)
#in_df$interconnection_fee_high <- ifelse(in_df$Year %in% NEM1_period, 0, intconn_NEM2_high)
#in_df$interconnection_fee_low  <- ifelse(in_df$Year %in% NEM1_period, 0, intconn_NEM2_low)
# Group and projection
#in_df$geometry = NULL
#interconnection_fee_df <- in_df %>% group_by(Year) %>% summarise(interconnection_fee_base = sum(interconnection_fee_base, na.rm = TRUE),
#                                                                 interconnection_fee_high = sum(interconnection_fee_high, na.rm = TRUE),
#                                                                 interconnection_fee_low = sum(interconnection_fee_low, na.rm = TRUE))
#interconnection_fee_df[, c(2:4)] <- interconnection_fee_df[, c(2:4)] %>% cumsum()
# Get total fee associated with each array over time
#interconnection_fee_base <- c(interconnection_fee_df$interconnection_fee_base, rep(tail(interconnection_fee_df$interconnection_fee_base, 1), system_lifespan-1))
#interconnection_fee_best <- c(interconnection_fee_df$interconnection_fee_low, rep(tail(interconnection_fee_df$interconnection_fee_low, 1), system_lifespan-1))
#interconnection_fee_worst <- c(interconnection_fee_df$interconnection_fee_high, rep(tail(interconnection_fee_df$interconnection_fee_high, 1), system_lifespan-1))
# If missing years, include zeros for interconnection fees
#interconnection_fee_base <- c(interconnection_fee_base, rep(0, end_install_year - max(in_df$Year)))
#interconnection_fee_best <- c(interconnection_fee_best, rep(0, end_install_year - max(in_df$Year)))
#interconnection_fee_worst <- c(interconnection_fee_worst, rep(0, end_install_year - max(in_df$Year)))
# Create final dataframe
#interconnection_fee_df <- data.frame(Year = c(start_install_year:(end_install_year+system_lifespan-1)), 
#                                     base = interconnection_fee_base, 
#                                     best = interconnection_fee_best, 
#                                     worst = interconnection_fee_worst)



# Supplementary Discussion 2: Proportion of arrays larger than 1 MWpand single axis
propUtilSingleAxis = (in_solar_df[which(in_solar_df$Capacity>=capacity_threshold & in_solar_df$Class=="single_axis"), ] %>% nrow()) / (in_solar_df[which(in_solar_df$Capacity>=capacity_threshold), ] %>% nrow()) * 100



# Supplementary Discussion 3: Comparison of 5mw cutoff
# Shifted kcal
commFEW <- read.csv(paste("Outputs/Dataframes/ResourceModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
FEW <- read.csv(paste("Outputs/Dataframes/FEW_Resource_Impacts_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
FEW_5 <- read.csv(paste("Outputs/Dataframes/FEW_Resource_Impacts_", system_lifespan, "yrLS_5MW_", "sensitivity", ".csv", sep = ""))
baseline_Bkcal = FEW$Base[which(FEW$Resource=="Shifted Calories" & FEW$Scale=="Comm")]
baseline5_Bkcal = FEW_5$Base[which(FEW_5$Resource=="Shifted Calories" & FEW_5$Scale=="Comm")]

# Scale shifted arrays
scaleshifted_df = in_solar_df[which(in_solar_df$Index %in% commFEW_5$Index & !in_solar_df$Index %in% commFEW$Index), ]
scaleshifted_area = sum(scaleshifted_df$dir_a) / 1e6
scaleshifted_cap = sum(scaleshifted_df$Capacity)

# BIG NUMBER: Final budget / Food revenue loss for 5MW compared to baseline 1 MW
# Baseline
commecon <- read.csv(paste("Outputs/Dataframes/EconModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
commecon_5 <- read.csv(paste("Outputs/Dataframes/EconModel_CommArrayResults_", system_lifespan, "yrLS_", 5, "MW_", "sensitivity", ".csv", sep = ""))
# 1 MW
commBudgSuperiority = commecon %>% group_by(Year) %>% summarise(budget = sum(Tot_Budg_base), food = sum(food_base))
finalEconComm_ratio = tail(commBudgSuperiority$budget, 1) / tail(commBudgSuperiority$food, 1)
# 5 MW
commBudgSuperiority5 = commecon_5 %>% group_by(Year) %>% summarise(budget = sum(Tot_Budg_base, na.rm = TRUE), food = sum(food_base, na.rm = TRUE))
finalEconComm_ratio5 = tail(commBudgSuperiority5$budget, 1) / tail(commBudgSuperiority5$food, 1)




# SUPPLEMENTARY DISCSSION 3: Arrays with negative budgets
# Get economic payback times (commercial, utility, and with half operational costs)
commecon <- read.csv(paste("Outputs/Dataframes/EconModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
utilecon <- read.csv(paste("Outputs/Dataframes/EconModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
# Get payback time vector with base methods
getPaybackTimes = function(df){
  df_temp <- df %>% dplyr::select(matches("Index|Year|Yr_inst|Tot_"))
  # For each index find base, best, and worst payback time and save out
  payback_df <- data.frame(base = as.numeric(), best = as.numeric(), worst = as.numeric(), Yr_inst = as.integer(), Index = as.integer())
  for(indx in unique(df_temp$Index)){
    payback_temp <- df_temp[which(df_temp$Index == indx), ]
    yr_inst = first(payback_temp$Yr_inst, na_rm = TRUE)
    payback_temp <- payback_temp[which(payback_temp$Year>=yr_inst), ]
    # Gets dataframe where row are only those with positive total econ values
    payback_base <- payback_temp %>% group_by(Year) %>% mutate(Inx = first(which(Tot_Budg_base > 0))) %>% filter(Inx == row_number()) %>% dplyr::select(-Inx) 
    payback_base <- as.numeric(min(payback_base$Year) - yr_inst)
    payback_best <- payback_temp %>% group_by(Year) %>% mutate(Inx = first(which(Tot_Budg_max > 0))) %>% filter(Inx == row_number()) %>% dplyr::select(-Inx) 
    payback_best <- as.numeric(min(payback_best$Year) - yr_inst)
    payback_worst <- payback_temp %>% group_by(Year) %>% mutate(Inx = first(which(Tot_Budg_min > 0))) %>% filter(Inx == row_number()) %>% dplyr::select(-Inx) 
    payback_worst <- as.numeric(min(payback_worst$Year) - yr_inst)
    payback_df <- rbind(payback_df, data.frame(base = payback_base, best = payback_best, worst = payback_worst, Yr_inst = yr_inst, Index = indx))
  }
  return(payback_df)
}
# Run for comm and util under standard conditions (regular operational costs)
comm_payback = getPaybackTimes(commecon)
util_payback = getPaybackTimes(utilecon)
# Find number of installations never reading payback - USED IN SUPPLEMENTAL
comm_paybackNONPROFIT <- comm_payback[is.infinite(comm_payback$base),]
comm_paybackNONPROFITbest <- comm_payback[is.infinite(comm_payback$best),]
avgSize_nonProfit = in_solar_df$Capacity[which(in_solar_df$Index %in% comm_paybackNONPROFITbest$Index)] %>% mean()
# Utility-scale check
util_paybackNONPROFIT <- util_payback[is.infinite(util_payback$base),] 
util_paybackNONPROFITbest <- util_payback[is.infinite(util_payback$best),]



# Klise et al reported water use factors
wu_acft_thous_m3_best = 0.23 / acft_m3 * 1000
wu_acft_thous_m3_worst = 2.16 / acft_m3 * 1000
wu_acft_thous_m3_base = 1.20 / acft_m3 * 1000


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get average electricity rate (sensivity analysis, and utility scale water value)

# Adjust all econ variables to inflation
econ_df = in_solar_df[which(in_solar_df$Capacity<1), ] %>% dplyr::select(matches("Gen|Econ|Year"))
econ_df$geometry = NULL

# Adjust for inlfation to 2018 (real 2018 USD cost of electricity)
yr_list = unique(econ_df$Year)
rate_df <- data.frame(rate = numeric())
for(year in yr_list){
  # Get colname
  econ_colName = paste("Econ_", year, sep = "")
  # Get multiplyer
  infl_rate_mult <- inflation_rateCPI$infl_rate_mult[which(inflation_rateCPI$Year==year)]
  econ_df[,which(colnames(econ_df)==econ_colName)] <- econ_df[,which(colnames(econ_df)==econ_colName)] * infl_rate_mult
  # Get rates and average -- $/kWh
  gen_colName = paste("Gen_", year, sep = "")
  rate_df_temp <- (econ_df[,which(colnames(econ_df)==econ_colName)] / econ_df[,which(colnames(econ_df)==gen_colName)] / kW_MW) %>% mean(na.rm=TRUE)
  rate_df <- rbind(rate_df, rate_df_temp)
}

# Get rate
rate <- rate_df[,1] %>% mean()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get average rate of inflation 

# PPI
infl_rate = inflation_ratePPI %>% mutate(pct_change = (Annual/lag(Annual) - 1) * 100)
infl_ratePPI = infl_rate$pct_change %>% mean(na.rm=TRUE)
# CPI
infl_rate = inflation_rateCPI %>% mutate(pct_change = (Annual/lag(Annual) - 1) * 100)
infl_rateCPI = infl_rate$pct_change %>% mean(na.rm=TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Supplemental figure of individual payback time

# Get dfs
commecon <- read.csv(paste("Outputs/Dataframes/EconModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
utilecon <- read.csv(paste("Outputs/Dataframes/EconModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
# Get only operational years
getPaybackTime_fig = function(df, PVscale){
  df_temp <- df %>% dplyr::select(matches("Index|Year|Yr_inst|Tot_"))
  # For each index find base, best, and worst payback time and save out
  payback_df <- df_temp[which(df_temp$Index=="Get Empty Dataframe"), ]
  for(indx in unique(df_temp$Index)){
    payback_temp <- df_temp[which(df_temp$Index == indx), ]
    yr_inst = first(payback_temp$Yr_inst, na_rm = TRUE)
    payback_temp <- payback_temp[which(payback_temp$Year>=yr_inst), ]
    payback_temp$Year = c(1:nrow(payback_temp))
    payback_df = rbind(payback_df, payback_temp[which(payback_temp$Year<=system_lifespan), ]) # removes repeat cumulative values
  }
  
  # Plot variables
  # Plot baseline
  if(PVscale == "Commercial"){breaks = seq(-100,100,2)}else{breaks = seq(-100,100,1)}
  if(PVscale == "Commercial"){limits = c(-2,6)}else{limits = c(-0.5,4.5)}
  adjust = 1e6
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BASELINE
  # Get median and first and third quartiles (BASELINE)
  out_df = data.frame(Year = numeric(), mean = numeric(), first = numeric(), third = numeric(), ten = numeric(), ninety = numeric())
  for(year in c(1:system_lifespan)){
    #year_summary = summary(payback_df$Tot_Budg_base[which(payback_df$Year==year)])
    year_quantile = quantile(payback_df$Tot_Budg_base[which(payback_df$Year==year)], probs = c(0.10, 0.25, 0.5, 0.75, 0.90))
    out_df = rbind(out_df, data.frame(Year = year, mean = year_quantile[3], first = year_quantile[2], third = year_quantile[4], ten = year_quantile[1], ninety = year_quantile[5]))
  }
  basePayback_plot <- ggplot(out_df, aes(x=Year)) + 
    geom_ribbon(aes(ymin = ten/adjust, ymax = first/adjust), fill = "green", alpha = .25) +
    geom_ribbon(aes(ymin = third/adjust, ymax = ninety/adjust), fill = "green", alpha = .25) +
    geom_ribbon(aes(ymin = first/adjust, ymax = third/adjust), fill = "blue", alpha = .25) +
    geom_line(aes(y = mean/adjust), size = 1.5, color = "black") +
    geom_hline(yintercept = 0, size = 0.75, color = "black") +
    scale_x_continuous(name = "Years Since Installation", breaks=seq(0, 25, 5), limits = c(1, 25), expand = c(0, 0)) +
    scale_y_continuous(name = paste("Net-", PVscale, " Budget (Mil USD)", sep=""), breaks = breaks, limits = limits, expand = c(0, 0)) +
    theme(plot.title = element_blank(), legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(color = "black", fill = NA, size = border_size))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BEST
  # Get median and first and third quartiles (BASELINE)
  out_df = data.frame(Year = numeric(), mean = numeric(), first = numeric(), third = numeric(), ten = numeric(), ninety = numeric())
  for(year in c(1:system_lifespan)){
    #year_summary = summary(payback_df$Tot_Budg_base[which(payback_df$Year==year)])
    year_quantile = quantile(payback_df$Tot_Budg_max[which(payback_df$Year==year)], probs = c(0.10, 0.25, 0.5, 0.75, 0.90))
    out_df = rbind(out_df, data.frame(Year = year, mean = year_quantile[3], first = year_quantile[2], third = year_quantile[4], ten = year_quantile[1], ninety = year_quantile[5]))
  }
  bestPayback_plot <- ggplot(out_df, aes(x=Year)) + 
    geom_ribbon(aes(ymin = ten/adjust, ymax = first/adjust), fill = "green", alpha = .25) +
    geom_ribbon(aes(ymin = third/adjust, ymax = ninety/adjust), fill = "green", alpha = .25) +
    geom_ribbon(aes(ymin = first/adjust, ymax = third/adjust), fill = "blue", alpha = .25) +
    geom_line(aes(y = mean/adjust), size = 1.5, color = "black") +
    geom_hline(yintercept = 0, size = 0.75, color = "black") +
    scale_x_continuous(name = "Years Since Installation", breaks=seq(0, 25, 5), limits = c(1, 25), expand = c(0, 0)) +
    scale_y_continuous(name = paste("Net-", PVscale, " Budget (Mil USD)", sep=""), breaks = breaks, limits = limits, expand = c(0, 0)) +
    theme(plot.title = element_blank(), legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(color = "black", fill = NA, size = border_size))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ WORST
  # Get median and first and third quartiles (BASELINE)
  out_df = data.frame(Year = numeric(), mean = numeric(), first = numeric(), third = numeric(), ten = numeric(), ninety = numeric())
  for(year in c(1:system_lifespan)){
    #year_summary = summary(payback_df$Tot_Budg_base[which(payback_df$Year==year)])
    year_quantile = quantile(payback_df$Tot_Budg_min[which(payback_df$Year==year)], probs = c(0.10, 0.25, 0.5, 0.75, 0.90))
    out_df = rbind(out_df, data.frame(Year = year, mean = year_quantile[3], first = year_quantile[2], third = year_quantile[4], ten = year_quantile[1], ninety = year_quantile[5]))
  }
  worstPayback_plot <- ggplot(out_df, aes(x=Year)) + 
    geom_ribbon(aes(ymin = ten/adjust, ymax = first/adjust), fill = "green", alpha = .25) +
    geom_ribbon(aes(ymin = third/adjust, ymax = ninety/adjust), fill = "green", alpha = .25) +
    geom_ribbon(aes(ymin = first/adjust, ymax = third/adjust), fill = "blue", alpha = .25) +
    geom_line(aes(y = mean/adjust), size = 1.5, color = "black") +
    geom_hline(yintercept = 0, size = 0.75, color = "black") +
    scale_x_continuous(name = "Years Since Installation", breaks=seq(0, 25, 5), limits = c(1, 25), expand = c(0, 0)) +
    scale_y_continuous(name = paste("Net-", PVscale, " Budget (Mil USD)", sep=""), breaks = breaks, limits = limits, expand = c(0, 0)) +
    theme(plot.title = element_blank(), legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(color = "black", fill = NA, size = border_size))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BRING TOGETHER
  # Arrange & Save
  fig <- ggarrange(bestPayback_plot+theme(axis.title.y = element_blank()), basePayback_plot, worstPayback_plot+theme(axis.title.y = element_blank()), nrow = 3, ncol = 1, heights = c(1,1,1), align = c("hv")) 
  ggsave(path = 'Outputs/Figures', filename = paste("PaybackPeriod_", PVscale, "_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".pdf", sep=""), width = 3.5, height = 7, units = "in", useDingbats=FALSE)
  # Return fig
  return("Finished.")
}

# Get payback df for both util and comm
getPaybackTime_fig(commecon, "Commercial")
getPaybackTime_fig(utilecon, "Utility")

# Logic check
sum(commecon$Tot_Budg_base[which(commecon$Year==2042)]/1e9)
sum(utilecon$Tot_Budg_base[which(utilecon$Year==2042)]/1e9)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Supplemental figure of scenario descriptsion

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RESOURCE

# Set to three significant figures for base/best/worst
df_food = data.frame(Resoruce = "Displaced Food (kcal)", 
                     Baseline = "Low-yielding land (6.5% yield deficit) was displaced with no change in average yields through time.", 
                     Best = "Extremely low-yielding land (25% yield deficit) was displaced with decreases (-0.1%/yr) in average yields through time.", 
                     Worst = "Average-yielding land (0% yield deficit) was displaced with increases (0.8%/yr) average yields through time.")
df_energy = data.frame(Resoruce = "Generated Energy (GWh)", 
                       Baseline = "Installed PV modules were a mix of mono-Si and muli-Si with the avearge reported efficiency for a given installation year (Supplementary Fig. 4).", 
                       Best = "Installed PV modules were entirely mono-Si extrapolated efficiency for a given installation year (Supp. Fig. 4).", 
                       Worst = "Installed PV modules were a entirely muli-Si with extrapolated efficiency for a given installation year (Supp. Fig. 4).")
df_water = data.frame(Resoruce = "Saved Irrigation Water (mil-m3)", 
                      Baseline = "Average irrigation depth prediction between FRIS and IWMS, and avearge annual precipitation during the addition phase projected forward.", 
                      Best = "Maximum irrigation depth prediction between FRIS and IWMS, and highest offset irrigation water use due to estimted precipitaiton of driest year during the addition phase projected forward.", 
                      Worst = "Minimum irrigation depth prediction between FRIS and IWMS, and lowest offset irrigation water use due to estimted precipitaiton of wettest year during the addition phase projected forward.")
df_irrigEn = data.frame(Resoruce = "Offset Irrigation Energy (GWh)", 
                        Baseline = "Same as Saved Irrigation Water baseline, with county-level irrigation energy requirement estimation.", 
                        Best = "Same as Saved Irrigation Water best-case, with county-level irrigation energy requirement estimation.", 
                        Worst = "Same as Saved Irrigation Water worst-case, with county-level irrigation energy requirement estimation.")
df_oandmwu = data.frame(Resoruce = "O&M Water Use (mil-m3)", 
                        Baseline = "Average O&M water use for dry cooled systems (0.97 thousand-m3/MW/yr) from Klise et al. (2013).", 
                        Best = "Reproted low O&M water use for dry cooled systems (0.19 thousand-m3/MW/yr) from Klise et al. (2013).", 
                        Worst = "Reproted high O&M water use for dry cooled systems (1.75 thousand-m3/MW/yr) from Klise et al. (2013).")
df = rbind(df_food, df_energy, df_water, df_irrigEn, df_oandmwu)
# Clean pub datatable
library(gt)
df_export <- df %>% gt()
gtsave(df_export, paste("Outputs/Dataframes/Resource_Scenario_Descriptions.png", sep = ""), vwidth = 1000, vheight = 1800)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ECON

# Set to three significant figures for base/best/worst
df_nem = data.frame(Resoruce = "Cash Flow Saved (Offset Load) and Earned (Surplus Generation NEM), Commercial-Scale", 
                    Baseline = "Same scenario as baseline Generated Energy (Supp. Table 18), with utility-specific energy charge rates and Annual Energy Outlook projected changes in electricity prices.", 
                    Best = "Same scenario as best-case Generated Energy (Supp. Table 18), with utility-specific energy charge rates and Annual Energy Outlook projected changes in electricity prices.", 
                    Worst = "Same scenario as worst-case Generated Energy (Supp. Table 18), with utility-specific energy charge rates and Annual Energy Outlook projected changes in electricity prices.")
df_landlease = data.frame(Resoruce = "Cash Flow from Annual Land Lease Contract, Utility-Scale", 
                          Baseline = "Average consultant and industry reported land lease contract rate ($2,450/ha/yr).", 
                          Best = "High consultant and industry reported land lease contract rate ($4,950/ha/yr).", 
                          Worst = "Low consultant and industry reported land lease contract rate ($750/ha/yr).")
df_food = data.frame(Resoruce = "Forgone Crop Cash Flow", 
                     Baseline = "Same scenario as baseline Displaced Food (Supp. Table 18), with food prices scaling with projected changes in electricity prices.", 
                     Best = "Same scenario as best-case Displaced Food (Supp. Table 18), with food prices scaling with projected changes in electricity prices.", 
                     Worst = "Same scenario as worst-case Displaced Food (Supp. Table 18), with food prices scaling with projected changes in electricity prices.")
df_operation = data.frame(Resoruce = "Saved Cash Flow from Business as Usual Farm Operation", 
                          Baseline = "Same scenario as baseline Forgone Crop Cash Flow, with costs represented in USD/kg to produce.", 
                          Best = "Same scenario as best-case Forgone Crop Cash Flow, with costs represented in USD/kg to produce.", 
                          Worst = "Same scenario as worst-case Forgone Crop Cash Flow, with costs represented in USD/kg to produce.")
df_water = data.frame(Resoruce = "Cash Flow from Change in Water Budget", 
                      Baseline = "Baseline scenario for Saved Irrigation Water minus baseline scenario for O&M Water Use (Supp. Table 18), with constant irrigation energy requirement cost (Cash flow saved and earned). with a CCV-wide water right rate.", 
                      Best = "Best-case scenario for Saved Irrigation Water minus best-case scenario for O&M Water Use (Supp. Table 18), with constant irrigation energy requirement cost (Cash flow saved and earned). with a CCV-wide water right rate.", 
                      Worst = "Worst-case scenario for Saved Irrigation Water minus worst-case scenario for O&M Water Use (Supp. Table 18), with constant irrigation energy requirement cost (Cash flow saved and earned). with a CCV-wide water right rate.")
df_oandm= data.frame(Resoruce = "Cash Flow to O&M", 
                     Baseline = "Average O&M costs from NREL Cost Benchmark. Projections of fixed O&M costs under the 'moderate' Annual Technology Basline scenario", 
                     Best = "Lowest O&M costs from NREL Cost Benchmark. Projections of fixed O&M costs under the 'advanced' Annual Technology Basline scenario", 
                     Worst = "Highest O&M costs from NREL Cost Benchmark. Projections of fixed O&M costs under the 'conservative' Annual Technology Basline scenario")
df_installation = data.frame(Resoruce = "Initial Installation Cost", 
                             Baseline = "Median reported installation costs to PV system owner from Barbose et al. (2021) including an assumed 30% Solar ITC reduction.", 
                             Best = "20th percentile reported installation costs to PV system owner from Barbose et al. (2021) including an assumed 30% Solar ITC reduction.", 
                             Worst = "80th percentile reported installation costs to PV system owner from Barbose et al. (2021) including an assumed 30% Solar ITC reduction.")
df = rbind(df_nem, df_landlease, df_food, df_operation, df_water, df_oandm, df_installation)
# Clean pub datatable
library(gt)
df_export <- df %>% gt()
gtsave(df_export, paste("Outputs/Dataframes/Economic_Scenario_Descriptions.png", sep = ""), vwidth = 1000, vheight = 1800) # not sure why these heights and widths work



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Supplemental figure of efficiency LM, and export for econ model

# ----------------------------------------------- #
# Get linear model of mono vs multi for scenarios #  -- this is done within the model, based on updated mono share and efficiency data, but exported here
# ----------------------------------------------- #
# Set indf
in_df <- in_solar_df
# Non user variables 
start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
end_install_year  = in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)
# SETUP Generated Variables -- no input necessary
install_period = end_install_year - start_install_year
model_length = install_period + system_lifespan 
# Call eff_df
energy_efficiency_df_loc = 'Data/Downloaded/efficiency_yearly_LBerkeley.csv'
eff_df <- read.csv(energy_efficiency_df_loc)
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
write.csv(eff_df_moml, "Data/Derived/efficiency_LBNL_08to18_monopoly.csv", row.names = FALSE) # to econ model

# Plot eff
eff_plot <- ggplot(eff_df_moml, aes(x=Year)) + 
  geom_line(aes(y=efficiency_mono*100), size = 1, color = "blue") +
  geom_line(aes(y=efficiency_multi*100), size = 1, color = "red") +
  geom_line(aes(y=efficiency_reported*100), size = 1, color = "black") +
  xlab("Year") + 
  scale_x_continuous(breaks=seq(2008, 2018, 2), expand = c(0, 0)) +#, limits = c(2007, 2044)) +
  scale_y_continuous(name = "Annual PV Efficiency (%)", breaks=seq(10,20,1), limits = c(12, 19), expand = c(0, 0)) +
  theme(plot.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y.left = element_text(color="black"), axis.title.y.right = element_text(color = "red"))
eff_plot
ggsave("eff_monopoly_plot.pdf", width = 3.5, height = 2.5, units = "in")
#rm(coeff, Energy_df, gt_convert)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Check if any crop types are missing crop dollar data

# Check crop type 
util <- read.csv(paste("Outputs/Dataframes/EconModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
comm <- read.csv(paste("Outputs/Dataframes/EconModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
U_na <- util[which(is.na(util$food_base)), ]
C_na <- comm[which(is.na(comm$food_base)), ]
df <- data.frame(Year = integer, food_base = numeric(), food_min = numeric(), food_max = numeric(), Index = integer())
for(i in c(1:nrow(in_solar_df))){
  df_temp <- getEconomic_food(in_solar_df[i, ])
  df_temp$Index <- in_solar_df[i, ]$Index
  df <- rbind(df, df_temp)
  print(i)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Extra figures and code for Solar-FEWE paper

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIGURE IN PAPER 

# NOTE THIS IS CURRENTLY SET TO RUN FOR COMMERICAL SCALE

# list all export files
file_list <- list.files("S:/Users/stidjaco/R_files/Chapter_2/Data/Derived/AdjacentLCC_PostInstall", pattern = "*csv", full.names = TRUE)
read_list <- lapply(file_list, read.csv)
solar <- do.call("rbind", read_list)
solar <- solar[which(solar$Index %in% in_solar_df$Index), ]

# Post LC year
solar$Capacity = in_solar_df$Capacity[match(solar$Index, in_solar_df$Index)]
solar$dir_a = in_solar_df$dir_a[match(solar$Index, in_solar_df$Index)]
solar$irrig = in_solar_df$irrig[match(solar$Index, in_solar_df$Index)]

# Subset for irrigated utility
#solar <- solar[which(solar$Capacity<capacity_threshold), ] # COMMERCIAL SCALE
solar <- solar[which(solar$Capacity>=capacity_threshold), ] # UTILITY SCALE
solar <- solar[which(solar$irrig==1), ]

# Get cdl crop
nass = Nass_Classifications
# Correct nass crops to modified FARS crops
crp_grp1 = "wheat|sorghum|grain|rice"
crp_grp2 = "vegetable|tomatoes|potatoes|beans|corn|lettuce|cotton|crops, other"
crp_grp3 = "hay|pastureland"
crp_grp4 = "berry totals|orchards"
crp_nme1 = "grain"
crp_nme2 = "veg/cotton/other"
crp_nme3 = "hay/pasture"
crp_nme4 = "orchards/grapes"
nass$sankey_crop <- tolower(nass$Irrig_Crop)
# Group further to simplify
nass$sankey_crop <- ifelse(grepl(crp_grp1,nass$sankey_crop), crp_nme1, nass$sankey_crop)
nass$sankey_crop <- ifelse(grepl(crp_grp2,nass$sankey_crop), crp_nme2, nass$sankey_crop)
nass$sankey_crop <- ifelse(grepl(crp_grp3,nass$sankey_crop), crp_nme3, nass$sankey_crop)
nass$sankey_crop <- ifelse(grepl(crp_grp4,nass$sankey_crop), crp_nme4, nass$sankey_crop)
# Add custom for fallow idle and all other non-ag
nass$sankey_crop <- ifelse(nass$Value %in% c(61, 65), "fallow/idle", nass$sankey_crop)
nass$sankey_crop <- ifelse(is.na(nass$sankey_crop), "non-ag", nass$sankey_crop)
# Match solar with 
solar$pre_crop <- nass$sankey_crop[match(solar$cdl_be1, nass$Value)] %>% tolower()
solar$post_crop <- nass$sankey_crop[match(solar$cdl_af1, nass$Value)] %>% tolower()
solar <- solar[which(!is.na(solar$pre_crop) & !is.na(solar$post_crop)),]

# Remove those that were already fallowed or non ag surrounding before or after
solar <- solar[which(solar$pre_crop!="non-ag" & solar$post_crop!="non-ag"), ]

# ~~~~~~~~~~~~~~ Quantify amount of water associated with these arrays AND surrounding cropland -- USED IN PAPER
# Get arrays that shift to fallow/idle
solar_toFallow <- solar[which(solar$pre_crop!="fallow/idle" & solar$post_crop=="fallow/idle"), ]
solar_toFallow_Years <- solar_toFallow %>% group_by(Year) %>% summarize(num = n())
in_df <- in_solar_df[which(in_solar_df$Index %in% solar_toFallow$Index), ]
df <- data.frame(Index = integer(), irrig_offset = numeric())
for(i in c(1:nrow(in_df))){
  in_df_temp <- in_df[i, ]
  irrig_offset <- getResource_irrigWater(in_df_temp)
  df_temp <- data.frame(Index = in_df_temp$Index, irrig_offset = tail(irrig_offset$irrig_fg, 1))
  df <- rbind(df, df_temp)
}
# Get sum for array offset
tot_toFallow_offset = sum(df$irrig_offset/1000000)
tot_toFallow_capacity = sum(solar_toFallow$Capacity)
tot_toFallow_dir_a = sum(solar_toFallow$dir_a)/1000000
prop = tot_toFallow_dir_a / (sum(in_solar_df$dir_a[which(in_solar_df$Capacity>capacity_threshold)])/1000000) * 100

# ~~~~ GET FARM SHIFT TO FALLOW -- USED IN PAPER -- theorectical if proximal cropland was fallowed

# Get arrays that shift to fallow/idle
in_df <- in_solar_df[which(in_solar_df$Index %in% solar_toFallow$Index), ]
df <- data.frame(Index = integer(), irrig_offset = numeric(), farm_size = numeric())
for(i in c(1:nrow(in_df))){
  in_df_temp <- in_df[i, ]
  irrig_offset <- getResource_irrigWaterFarmFallow(in_df_temp)
  df_temp <- data.frame(Index = in_df_temp$Index, irrig_offset = tail(irrig_offset$irrig_fg, 1), farm_size = tail(irrig_offset$farm_size, 1))
  df <- rbind(df, df_temp)
}
# Get sum for array offset
tot_FARMtoFallow_offset = sum(df$irrig_offset/1000000)
tot_FARMtoFallow_offset_yr = tot_FARMtoFallow_offset / system_lifespan
tot_FARMtoFallow_dir_a = sum(df$farm_size)/1000000

# ~~~~~~~~~~~~~~

# Solar grouping
solar <- solar %>% group_by(pre_crop, post_crop) %>% summarise(num = n(), Capacity = sum(Capacity), dir_a = sum(dir_a)/1000000)

# Group by pre-crop -- USED IN PAPER
solar_cropprops <- solar %>% group_by(pre_crop) %>% summarise(area = sum(dir_a))
solar_cropprops <- solar_cropprops[which(solar_cropprops$pre_crop!="fallow/idle"), ]
solar_cropprops$prop <- solar_cropprops$area / sum(solar_cropprops$area) * 100

# Get newly idle arrays
newidle_num <- sum(solar$num[which(solar$pre_crop!="fallow/idle" & solar$post_crop=="fallow/idle")]) %>% as.numeric()
newidle_area <- (sum(solar$dir_a[which(solar$pre_crop!="fallow/idle" & solar$post_crop=="fallow/idle")])) %>% as.numeric() 
#newidle_area_prop <- (sum(solar$dir_a[which(solar$pre_crop!="fallow/idle" & solar$post_crop=="fallow/idle")])) / sum(solar$dir_a) * 100
newidle_capacity <- sum(solar$Capacity[which(solar$pre_crop!="fallow/idle" & solar$post_crop=="fallow/idle")]) %>% as.numeric()

# Get prop of each crop type - USED IN PAPER
toFallow_cropprop = solar[which(solar$pre_crop!="fallow/idle" & solar$post_crop=="fallow/idle"), ] 
toFallow_cropprop$prop = toFallow_cropprop$d

# Temp for simplification
#solar_plot <- solar[which(solar$Capacity>10), ]
solar_plot = solar

# Plot cols
#crop_plot_cols_sankey = c("fallow/idle" = "black", 
#                          "grain" = "peru", 
#                          "hay/pasture" = "khaki2",
#                          "orchards/grapes" = "purple3", # Modified from other plots for simplicity
#                          "veg/cotton/other" = "darkorange1") 
crop_plot_cols_sankey = c("fallow/idle" = "black", 
                          "grain" = "peru", 
                          "hay/pasture" = "khaki2",
                          "orchards/grapes" = "mediumvioletred", # Modified from other plots for simplicity
                          "veg/cotton/other" = "seagreen") # Modified from other plots for simplicity
# Create Sankey
library(ggalluvial)
ggplot(data = solar_plot, aes(axis1 = pre_crop, axis2 = post_crop, y = dir_a)) +
  scale_x_discrete(limits = c("Proximal Cropland Pre-Solar", "Proximal Cropland Post-Solar"), expand = c(0.15, .05)) +
  ylab("Utility Solar: Irrigated Area Converted (km2)") +
  geom_alluvium(aes(fill = pre_crop),  alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  #scale_y_continuous(breaks=seq(0, 30, 3), limits = c(-0.5,30), expand = c(0, 0)) +
  scale_y_continuous(breaks=seq(0, 3, 0.1), limits = c(0,2.3), expand = c(0, 0)) +
  scale_fill_manual(values = crop_plot_cols_sankey, name = "Crop") +
  theme_minimal() +
  theme(#panel.grid.major = element_blank(),
    #panel.border = element_rect(colour="black", fill=NA, size = 1),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
#ggsave("Sankey_PrePost_SolarCropland.pdf", width = 7, height = 5, useDingbats=FALSE)
ggsave("Sankey_PrePost_SolarCropland_comm.pdf", width = 7, height = 5, useDingbats=FALSE)

# Get datatable for this too
# Set to three significant figures for base/best/worst
solar$Capacity <- solar$Capacity %>% round(0) #signif(3) %>% format(scientific = F)
solar$dir_a <- solar$dir_a %>% round(2) #signif(3) %>% format(scientific = F)
solar <- solar %>% arrange(desc(solar$dir_a))
write.csv(solar, "Outputs/Dataframes/shift_fallow_df_comm.csv")
#write.csv(solar, "Outputs/Dataframes/shift_fallow_dfcsv")

# Clean pub datatable
library(gt)
solar$Transition = paste(solar$pre_crop, solar$post_crop, sep = " to ")
solar <- solar %>% relocate(Transition)
names(solar)[names(solar) == "dir_a"] <- "Area"
solar$post_crop = NULL
solar$pre_crop = NULL
df_export <- solar %>% gt()
gtsave(df_export, "Outputs/Dataframes/shift_fallow_df.png")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Graphs for idle/fallow land and active ag land
idle <- read.csv("/Data/Derived/cv_fallow_idle.csv")
#idle21 <- read.csv("Data/Derived/cv_fallow_idle_thru2021.csv")
#idle <- rbind(idle, idle21)
# PLOT 
fallowidle_plot <- ggplot(idle) +
  stat_smooth(aes(x=Year, y=idle_ar), method = "lm", col = "black") +
  geom_point(aes(x=Year, y=idle_ar), size=4, color="brown", shape=19) +
  theme_minimal() +
  #geom_vline(aes(xintercept=2018)) +
  xlab("Year") + 
  ylab("Fallow/Idle Area (km2)") + 
  scale_x_continuous(breaks=seq(start_install_year-2, end_install_year+2, 2), limits = c(start_install_year-0.5, end_install_year+0.5), expand = c(0,0)) + 
  scale_y_continuous(breaks=seq(0, 8000, 1000), limits = c(3300, 7500), expand=c(0,0)) +
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_rect(colour="black", fill=NA, size = 1),axis.title=element_text(size=14,face="bold"))
fallowidle_plot
ggsave("fallowidle_plot.pdf", width = 3.5, height = 3.5, units = "in")

# Active ag area
active_ag <- read.csv(paste(wd, "/Data/Derived/cv_active_ag_area.csv", sep = ""))
#active_ag$ag_ar = active_ag$ag_ar + idle$idle_ar
# PLOT 
activeag_plot <- ggplot(active_ag) +
  stat_smooth(aes(x=Year, y=ag_ar), method = "lm", col = "black") +
  geom_point(aes(x=Year, y=ag_ar), size=4, color="blue", shape=19) +
  theme_minimal() +
  xlab("Year") + 
  ylab("Active Ag Area (km2)") + 
  scale_x_continuous(breaks=seq(start_install_year, end_install_year+3, 2), limits = c(start_install_year-0.5, end_install_year+3+0.5), expand = c(0,0)) + 
  scale_y_continuous(breaks=seq(30000, 41000, 1000), limits = c(36500, 40700), expand=c(0,0)) +
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_rect(colour="black", fill=NA, size = 1), axis.title=element_text(size=14,face="bold"))
activeag_plot

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Graph for Sustain Valencia
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

# PREPARE
ppt_df_year <- ppt_df_year[which(ppt_df_year$Year %in% USGS_yrs), ]
x = ppt_df_year$ppt[order(ppt_df_year$ppt)]
y = rev(ppt_df_year$ppt[order(ppt_df_year$ppt)])
y = y * runif(length(y),0.9,1.1)
model = lm(y~x)
# PLOT 
precipWu_plot <- ggplot() +
  stat_smooth(aes(x=x, y=y), method = "lm", col = "black") +
  geom_point(aes(x=x, y=y), size=4, color="blue", shape=19) +
  theme_minimal() +
  xlab("County Annual Precipitation") + 
  ylab("County Irrigation Water Use") + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_rect(colour="black", fill=NA, size = 1), axis.text.x = element_blank(), 
        axis.text.y = element_blank(), axis.title=element_text(size=14,face="bold"))
ggsave("Precip_irrigWU_concept.png", width = 3.5, height = 3.5, units = "in")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ NEW H3 plot for cv area


# Gett plotting stuff 
# Mapping
library("ggplot2")
library("ggspatial")
library("ggpubr")
library("rnaturalearth")
library("rnaturalearthdata")
library("maps")
library("smoothr")
library("tigris")
library("plotly")
library("vioplot")
library("RColorBrewer")
library("hexbin")
library("reshape2")

# Plotting crop groups
crp_grp1 = "wheat|sorghum|grain|rice"
crp_grp2 = "vegetable|tomatoes|potatoes|beans|corn"
crp_grp3 = "hay|pastureland"
crp_nme1 = "grain"
crp_nme2 = "vegetables"
crp_nme3 = "hay/pasture"
# Set margins
margin_plot <- unit(c(0,0,0,0),"cm")

# Get in_df
in_df <- in_solar_df
in_df <- st_transform(in_df, crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
# Get CV
CV_Boundary <- st_read("S:/Users/stidjaco/R_files/Chapter_1/Data/Downloaded/CV_alluvial_bnd/Alluvial_Bnd.shp") %>% st_as_sf() %>% st_transform(st_crs(in_df))
CV_Boundary <- fill_holes(CV_Boundary, threshold = units::set_units(50000000, 'km^2'))
# Run get croprotaion on in_df
in_df_save <- in_df
in_df$geometry = NULL
in_df <- in_df %>% getCrop_rotation()

# Add fars crop calsses
nass = Nass_Classifications
nass$Crop <- tolower(nass$Crop)
in_df$Crp_n_1 <- in_df$Crp_n_1 %>% tolower()
in_df$FARS_crop <- nass$Irrig_Crop[match(in_df$Crp_n_1, nass$Crop)] %>% tolower()
# Group further to simplify
in_df$FARS_crop <- ifelse(grepl(crp_grp1,in_df$FARS_crop), crp_nme1, in_df$FARS_crop)
in_df$FARS_crop <- ifelse(grepl(crp_grp2,in_df$FARS_crop), crp_nme2, in_df$FARS_crop)
in_df$FARS_crop <- ifelse(grepl(crp_grp3,in_df$FARS_crop), crp_nme3, in_df$FARS_crop)
in_df$Year <- as.character(in_df$Year)

# Turn fars crop into number
in_df$FARS_crop_num <- ifelse(in_df$FARS_crop=="berry totals", 1,
                              ifelse(in_df$FARS_crop=="cotton", 2,
                                     ifelse(in_df$FARS_crop=="crops, other", 3,
                                            ifelse(in_df$FARS_crop=="grain", 4, 
                                                   ifelse(in_df$FARS_crop=="hay/pasture", 5, 
                                                          ifelse(in_df$FARS_crop=="orchards", 6,
                                                                 ifelse(in_df$FARS_crop=="vegetables", 7, 8)))))))

# Get mode of fars crop
getmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]}
in_df <- in_df[!is.na(in_df$FARS_crop), ] %>% group_by(Index) %>% dplyr::summarise(FARS_crop = getmode(FARS_crop_num))

# Return geometry
in_df_save$FARS_crop <- in_df$FARS_crop[match(in_df_save$Index, in_df$Index)]

# Write to Arcpro
st_write(in_df_save, "fars_crop_fig_H3prep.shp")

'Arc Process: 
IN ArcGIS Online: 
- https://www.esri.com/arcgis-blog/products/arcgis-online/analytics/use-h3-hexagons-for-spatial-analysis-in-arcgis-online/ 
- Resolution: https://h3geo.org/docs/core-library/restable/#average-area-in-km2 
- Map Viewer, Analysis, Generate Tesselation, H3 Hexagon, Extent (CV), only intersecting, select resolution (7 is small, 6 is good for CV), hit run and wait

IN ARCPRO (OPTION 1): 
- Go to ArcPro and add data from "My Content"
- *Add Data* fars_crop_fig_H3prep.shp
- *Project* fars_crop_fig_H3prep.shp to Projected CS --> Continental --> NA --> NAD 1983 Contiguous USA Albers
- *Feature to Raster* FARS_crop field with result cell size to 30m raster
- *Zontal Statistics* in Spatial Analyst toolbox, Tesselation to Input raster or feature zone data, and raster resulting from last operation to raster, select "Majority"
- *Project Raster* to PrWGS 1984 --> GCS --> World --> WGS 1984 (or keep NAD1983 for visuals -- this is done for the paper)
- Change colors (below)
- Create layout with basemap of CV alluvium, overlain by Zonal Stats raster, overlain by H3 polygons empty fill

IN ARCPRO (OPTION 2): -- This is what we chose for the paper
- Go to ArcPro and add data from "My Content"
- *Add Data* fars_crop_fig_H3prep.shp
- *Project* fars_crop_fig_H3prep.shp to Projected CS --> Continental --> NA --> NAD 1983 Contiguous USA Albers
- *Feature to Raster* FARS_crop field with result cell size to 30m raster
- *Zontal Statistics* in Spatial Analyst toolbox, Tesselation to Input raster or feature zone data, and raster resulting from last operation to raster, select "Majority"
- Create polygon centroid layer using Feature to Point (Data Management) on H3 polygons
- Append raster value to centroid layer using Extract Multi Values to Points (Spatial Analyst)
- Join centroid layer with raster value to original H3 polygon layer.
- Change colors (below)
- Create layout with joined H3 shapes to export

Colors to Change to (Crop, ID, #HEX): 
- Grapes: 1 -- #7d26cd
- Cotton: 2 -- #b0e2ff
- Other crops: 3 -- #ff7f00
- Grain: 4 -- #cd853f
- Haypasture: 5 -- #eee685
- Orchard: 6 -- #ee0000
- Vegetables: 7 -- #9acd32

All done in ArcPro and Affinity' 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Old pixel map

# Get orchards
orchards <- Nass_Classifications[which(Nass_Classifications$Irrig_Crop=="Orchards"), ]
berries <- Nass_Classifications[which(Nass_Classifications$Irrig_Crop=="Berry totals"), ]

# Gett plotting stuff 
# Mapping
library("ggplot2")
library("ggspatial")
library("ggpubr")
library("rnaturalearth")
library("rnaturalearthdata")
library("maps")
library("smoothr")
library("tigris")
library("plotly")
library("vioplot")
library("RColorBrewer")
library("hexbin")
library("reshape2")

# Plotting crop groups
crp_grp1 = "wheat|sorghum|grain|rice"
crp_grp2 = "vegetable|tomatoes|potatoes|beans|corn"
crp_grp3 = "hay|pastureland"
crp_nme1 = "grain"
crp_nme2 = "vegetables"
crp_nme3 = "hay/pasture"
# Set margins
margin_plot <- unit(c(0,0,0,0),"cm")

# Get centroids
in_df <- in_solar_df
centroids <- st_centroid(in_df)
centroids <- st_transform(centroids, crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
coords <- as.data.frame(st_coordinates(centroids$geometry))
centroids$Lat <- coords$Y
centroids$Long <- coords$X
# Get CV
CV_Boundary <- st_read("S:/Users/stidjaco/R_files/Chapter_1/Data/Downloaded/CV_alluvial_bnd/Alluvial_Bnd.shp") %>% st_as_sf() %>% st_transform(st_crs(centroids))
CV_Boundary <- fill_holes(CV_Boundary, threshold = units::set_units(50000000, 'km^2'))
# Run get croprotaion on centroids
centroids_save <- centroids
centroids$geometry = NULL
centroids <- centroids %>% getCrop_rotation()

# Add fars crop calsses
nass = Nass_Classifications
nass$Crop <- tolower(nass$Crop)
centroids$Crp_n_1 <- centroids$Crp_n_1 %>% tolower()
centroids$FARS_crop <- nass$Irrig_Crop[match(centroids$Crp_n_1, nass$Crop)] %>% tolower()
# Group further to simplify
centroids$FARS_crop <- ifelse(grepl(crp_grp1,centroids$FARS_crop), crp_nme1, centroids$FARS_crop)
centroids$FARS_crop <- ifelse(grepl(crp_grp2,centroids$FARS_crop), crp_nme2, centroids$FARS_crop)
centroids$FARS_crop <- ifelse(grepl(crp_grp3,centroids$FARS_crop), crp_nme3, centroids$FARS_crop)
centroids$Year <- as.character(centroids$Year)

# Turn fars crop into number
centroids$FARS_crop_num <- ifelse(centroids$FARS_crop=="berry totals", 1,
                                  ifelse(centroids$FARS_crop=="cotton", 2,
                                         ifelse(centroids$FARS_crop=="crops, other", 3,
                                                ifelse(centroids$FARS_crop=="grain", 4, 
                                                       ifelse(centroids$FARS_crop=="hay/pasture", 5, 
                                                              ifelse(centroids$FARS_crop=="orchards", 6,
                                                                     ifelse(centroids$FARS_crop=="vegetables", 7, 8)))))))

# Get mode of fars crop
getmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]}
centroids <- centroids[!is.na(centroids$FARS_crop), ] %>% group_by(Index) %>% dplyr::summarise(FARS_crop = getmode(FARS_crop_num), Long = min(Long), Lat = min(Lat))

# Write to Arcpro
write.csv(centroids, "fars_crop_fig_centroids.csv")

'Arc Process: 
- *Add Data* fars_crop_fig_centroids.csv
- *Display XY Coordinates* fars_crop_fig_centroids.csv
- *Project* fars_crop_fig_centroids.csv to Projected CS --> Continental --> NA --> NAD 1983 Contiguous USA Albers
- *Feature to Raster* FARS_crop field with result cell size to 4km raster
- *Project Raster* to PrWGS 1984 --> GCS --> World --> WGS 1984
- *Export Raster* to call file below.'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Call in raster 
r <- raster("S:/Users/stidjaco/R_files/Chapter_2/CV_pv_crops_fig.tif")
r_spdf <- as(r, "SpatialPixelsDataFrame")
r_df <- as.data.frame(r_spdf)
colnames(r_df) <- c("value", "x", "y")

# For some reason, nodatavalue is filling as "2147483647", so set to NA if below
r_df$value = ifelse(r_df$value < 15, r_df$value, NA)

# Cast CV_boundary
CV_Boundary <- st_transform(CV_Boundary, crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') # %>% st_buffer(1000)

# New scale fill
crop_ras_cols = c("red2", "lightskyblue1", "yellowgreen", "peru", "khaki2","purple3", "darkorange1")
#Derive desired break/legend colours from gradient of selected brewer palette
crop_ras_cols <- colorRampPalette(crop_ras_cols, space="rgb")(7)

# Centroid plot with legend
pv_crop_plot <- ggplot() +
  geom_sf(data = CV_Boundary, fill="black") +
  geom_raster(data = r_df, aes(x=x, y=y, fill=value)) + 
  scale_fill_gradientn(colours=crop_ras_cols, na.value = "transparent") + 
  annotation_scale(location = "bl", width_hint = 0.75) +
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.25, "in"), pad_y = unit(0.7, "in"), style = north_arrow_fancy_orienteering) + 
  theme(rect = element_blank()) +
  theme(legend.title=element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank()) +
  theme(legend.title=element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none")
# Save
pv_crop_plot
ggsave("pv_crop_plot.pdf", width = 3.5, height = 5, units = "in", dpi = 300)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Utilities

# California Electric Utilitys 
ungroupPV_df <- in_solar_df

# Data source: https://gis.data.ca.gov/datasets/b95ca182aa254c3db8ad4d92bd32a73c_0?geometry=-161.087%2C31.071%2C-77.459%2C43.276
sf_use_s2(FALSE)
CA_utility <- st_read('S:/Users/stidjaco/R_files/Chapter_2/Data/Downloaded/CA_UtilService_shp/California_Electric_Utility_Service_Areas.shp') %>%
  st_transform(crs = st_crs(4326)) %>%
  st_buffer(dist = 0) %>%
  st_make_valid()
#st_intersection(CV_Boundary)
#st_difference() # overlapping?
CA_utility$Utility <- CA_utility$Utility %>% as.character()

# Get bounds for Sid
CA_utility$lon_min = NA
CA_utility$lon_max = NA
CA_utility$lat_min = NA
CA_utility$lat_max = NA
for (i in 1:nrow(CA_utility)){
  CA_utility[i, ]$lon_min <- st_bbox(CA_utility[i, ]$geometry)[1] %>% as.numeric()
  CA_utility[i, ]$lat_min <- st_bbox(CA_utility[i, ]$geometry)[2] %>% as.numeric()
  CA_utility[i, ]$lon_max <- st_bbox(CA_utility[i, ]$geometry)[3] %>% as.numeric()
  CA_utility[i, ]$lat_max <- st_bbox(CA_utility[i, ]$geometry)[4] %>% as.numeric()
}

# Intersect and Group shapes by utility 
# Had to do for loop because of st_intersection multipolygon error
df <- data.frame(ID=character(), Utility=character())
namesdf <- names(df)
for(i in 1:nrow(ungroupPV_df)){
  centroid <- st_centroid(ungroupPV_df[i, ])
  Utility <- st_intersection(CA_utility, centroid)
  idata <- t(data.frame(c(0, 0)))
  idata[1] <- i
  idata[2] <- Utility$Utility
  idata <- data.frame(idata)
  row.names(idata) <- i
  names(idata) <- namesdf
  names(df) <- namesdf
  df <- rbind(df, idata)
  print(paste('Completed ', i))
}
ungroupPV_df$Utility <- df$Utility
rm(df, i, namesdf, Utility, idata, centroid)

# Get arrays in certain utility
#utility_oi <- ungroupPV_df[which(ungroupPV_df$Utility == "Merced Irrigation District"), ]

# Match lat and long min and max
#ungroupPV_df$lon_min <- CA_utility$lon_min[match(ungroupPV_df$Utility, CA_utility$Utility)]
#ungroupPV_df$lon_max <- CA_utility$lon_max[match(ungroupPV_df$Utility, CA_utility$Utility)]
#ungroupPV_df$lat_min <- CA_utility$lat_min[match(ungroupPV_df$Utility, CA_utility$Utility)]
#ungroupPV_df$lat_max <- CA_utility$lat_max[match(ungroupPV_df$Utility, CA_utility$Utility)]

# Gruop by utility
CA_utility <- ungroupPV_df %>% 
  dplyr::group_by(Utility) %>% 
  summarise(Yr_inst = toString(unique(Year)), num = n(), Capacity = sum(Capacity))
CA_utility$geometry = NULL

# Get prop of capacity
CA_utility$Cap_Prop <- CA_utility$Capacity / sum(CA_utility$Capacity, na.rm = TRUE)

# Remove lat and long min and maxes from ungroupPVdf
ungroupPV_df$lon_min <- NULL
ungroupPV_df$lon_max <- NULL
ungroupPV_df$lat_min <- NULL
ungroupPV_df$lat_max <- NULL
ungroupPV_df$geometry = NULL

# Order CA_utility
CA_utility <- CA_utility[order(-CA_utility$Capacity), ]

# Export
#write.csv(ungroupPV_df, "CA_Utility_Service_PV.csv")
write.csv(CA_utility, "CA_Utility_ServiceFigure.csv")

# Prep for pretty export
library(gt)
CA_utility$Yr_inst = NULL
CA_utility$Cap_Prop = NULL
CA_utility$Capacity = CA_utility$Capacity %>% round(0)
df_export <- CA_utility %>% gt()
gtsave(df_export,"Outputs/Dataframes/CA_Utility_ServiceFigure.png")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test <- in_solar_df[which(in_solar_df$Capacity<1 & in_solar_df$irrig==1), ] 
test$geometry = NULL
test <- test %>% group_by(Year) %>% summarise(avg_size = median(dir_a))
plot(test$Year, test$avg_size)

# Plot vertical irrig vs non irrig operational energy requirements and overproduction
margin_plot <- unit(c(0,0,0,0),"cm")
irrig_over <- ggplot(irrig) + 
  geom_bar(aes(fill=load, y=value, x=Year), position=position_fill(reverse = TRUE), stat="identity") +
  ylab("Prop Overproduced") +
  scale_x_continuous(breaks=seq(start_install_year, end_install_year, 2), limits = c(start_install_year-0.5, end_install_year+0.5), expand = c(0,0)) + 
  scale_fill_manual(values=c("black", "green")) +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_rect(colour="black", fill=NA, size = 1), plot.margin = margin_plot) 
irrig_annCont <- ggplot(irrig) +
  geom_boxplot(aes(y=perc_dev_annload, x=as.character(Year)), width=0.5, outlier.shape=NA, fill = "blue") +
  geom_hline(yintercept=100) +
  ylab("% Annual Load") +
  scale_fill_manual(values=c("black", "green")) +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_rect(colour="black", fill=NA, size = 1), axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = margin_plot) + 
  coord_cartesian(ylim = c(0, 800))
irrig_plot <- ggarrange(irrig_annCont, irrig_over, nrow = 2, ncol = 1, align = c("hv"), heights = c(1,3))
ggsave(path = 'Outputs/Figures', filename = "irrig_annual_load_overpoduction.pdf", width = 3.5, height = 3.5, units = "in")

# Plot vertical irrig vs non irrig operational energy requirements and overproduction
margin_plot <- unit(c(0,0,0,0),"cm")
nonirrig_over <- ggplot(nonirrig) + 
  geom_bar(aes(fill=load, y=value, x=Year), position=position_fill(reverse = TRUE), stat="identity") +
  ylab("Prop Overproduced") +
  scale_x_continuous(breaks=seq(start_install_year, end_install_year, 2), limits = c(start_install_year-0.5, end_install_year+0.5), expand = c(0,0)) + 
  scale_fill_manual(values=c("black", "green")) +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_rect(colour="black", fill=NA, size = 1), plot.margin = margin_plot) 
nonirrig_annCont <- ggplot(nonirrig) +
  geom_boxplot(aes(y=perc_dev_annload, x=as.character(Year)), width=0.5, outlier.shape=NA, fill = "darkgoldenrod4") +
  geom_hline(yintercept=100) +
  ylab("% Annual Load") +
  scale_fill_manual(values=c("black", "green")) +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_rect(colour="black", fill=NA, size = 1), axis.text.x = element_blank(), axis.title.x = element_blank(), plot.margin = margin_plot) + 
  coord_cartesian(ylim = c(0, 2300))
nonirrig_plot <- ggarrange(nonirrig_annCont, nonirrig_over, nrow = 2, ncol = 1, align = c("hv"), heights = c(1,3))
ggsave(path = 'Outputs/Figures', filename = "non_irrig_annual_load_overpoduction.pdf", width = 3.5, height = 3.5, units = "in")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Merging and preprocessing Solar datasets

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~ Code to pre-process and stitch Stid et al. (2022) and Kruitwagen et al. (2021) Togegher~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
### Stid et al OMITTTED arrays addition ###
# Get missing arrays from Kruitwagen
####################
# Kruitwagen Solar #
####################

# Attribute info: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03957-7/MediaObjects/41586_2021_3957_MOESM1_ESM.pdf
kruit_pv <- read_sf("S:/Users/stidjaco/R_files/PhD_Ch1/Data/Downloaded/Kruitwagen_etal/predicted_set.geojson") %>% st_as_sf()
kruit_pv$ind_comm_10km = NULL
kruit_pv$wdpa_10km = NULL

# Subset US
PV_US <- kruit_pv[which(kruit_pv$iso.3166.2=="US-CA"), ]
PV_US <- st_transform(PV_US, crs = '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
# Subset CA and Central Valley
CV_Boundary <- st_cast(st_read("S:/Users/stidjaco/R_files/Chapter_1/Data/Downloaded/CV_alluvial_bnd/Alluvial_Bnd.shp"), "MULTIPOLYGON") %>% st_cast("POLYGON")
CV_Boundary <- fill_holes(CV_Boundary, threshold = units::set_units(50000000, km^2))
CV_Boundary <- as_Spatial(CV_Boundary)
CV_Boundary <- spTransform(CV_Boundary, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% st_as_sf() 
PV_US$CV <- st_intersects(PV_US, CV_Boundary) %>% as.numeric()
PV_CV <- PV_US[which(PV_US$CV==1),]
rm(PV_US)

# Get only necessary attributes
PV_CV <- select(PV_CV, matches("unique_id|area|capacity_mw"))
PV_CV$area_error <- NULL
colnames(PV_CV) <- c("Index", "dir_a", "cap_mw", "geometry")
PV_CV$Source <- 'Kruitwagen'

####################
#    Stid Solar    #
####################

# Call in our dataset
ungroupPV_df <- st_read('S:/Users/stidjaco/R_files/PhD_Ch1/Data/Downloaded/Stid_etal/Arrays/PV_ID_CV.shp')

# Subset to normalize
ungroupPV_df <- select(ungroupPV_df, matches("Index|Tot_a|TPVPp"))
colnames(ungroupPV_df)[colnames(ungroupPV_df)=='Tot_a'] <- 'dir_a'
colnames(ungroupPV_df)[colnames(ungroupPV_df)=='TPVPp'] <- 'cap_mw'
ungroupPV_df$Source <- 'Stid'

# Intersection between their dataset and ours
PV_CV_overlap <- st_intersection(PV_CV, ungroupPV_df)
'%!in%' <- Negate('%in%')
PV_stid_omit <- PV_CV[which(PV_CV$Index %!in% PV_CV_overlap$Index),]
rm(coords, CV_Boundary, PV_CV_overlap, kruit_pv, PV_CV, ungroupPV_df)

# Write Stid omitted arrays to shp to run in GEE renewables script
#st_write(PV_stid_omit, "PV_stid_omit.shp")
# Perform GEE Solar opeartions, manual yod interp, export and pull for Sid's model run, and pull back in and merge with our dataset -- POST initial GEE solar assessment
####################
#     Post-GEE     #
####################

# list all export files
file_list <- list.files("S:/Users/stidjaco/R_files/Chapter_2/Data/Derived/Kruit_Stid_Omit/GEE_Solar_Analysis", pattern = "*csv", full.names = TRUE)
read_list <- lapply(file_list, read.csv)
Kruit_Stid_Omit <- do.call("rbind", read_list)
Kruit_Stid_Omit <- Kruit_Stid_Omit %>% distinct(Index, .keep_all = TRUE) # Had to be done due to gee error with random export sampling 
stidomit <- Kruit_Stid_Omit

# Remove rooftop, missing panels, and non-ag, and normalize for Electricity model export
rooftop_num <- Kruit_Stid_Omit[which(Kruit_Stid_Omit$rooftop==1), ] %>% nrow() %>% as.numeric()
Kruit_Stid_Omit <- Kruit_Stid_Omit[which(Kruit_Stid_Omit$rooftop==0), ]
panels_missing <- Kruit_Stid_Omit[which(Kruit_Stid_Omit$mount==-9999), ] %>% nrow() %>% as.numeric()
Kruit_Stid_Omit <- Kruit_Stid_Omit[which(Kruit_Stid_Omit$mount!=-9999), ]
non_ag <- c(61, 63, 64, 65, 81, 82, 83, 87, 88, 92, 111, 112, 121, 122, 123, 124, 131, 141, 142, 143, 152, 190, 195) # Same classes use in chapter I
non_ag_kruit <- Kruit_Stid_Omit[which(Kruit_Stid_Omit$cdl_1yr %in% non_ag), ] %>% nrow() %>% as.numeric()
Kruit_Stid_Omit <- Kruit_Stid_Omit[which(!Kruit_Stid_Omit$cdl_1yr %in% non_ag), ]
Kruit_Stid_Omit <- Kruit_Stid_Omit %>% select(matches("mount|pnl_a|Yr_inst|Index|cap_mw"))
names(Kruit_Stid_Omit) <- c("Index", "Yr_inst", "pnl_a", "Class", "TPVPp")

# Give kruit_stid_omit lat long
shp <- st_read('S:/Users/stidjaco/R_files/PhD_Ch1/Data/Derived/CONUS_Renewables_PreGEE/CONUS_Renewables.shp')
shp <- shp[which(shp$Index %in% Kruit_Stid_Omit$Index), ] %>% select(matches("Index"))
Kruit_Stid_Omit <- merge(shp, Kruit_Stid_Omit, by=c('Index'))
coords <- as.data.frame(st_coordinates(st_centroid(Kruit_Stid_Omit$geometry)))
Kruit_Stid_Omit$Lat <- coords$Y
Kruit_Stid_Omit$Long <- coords$X
Kruit_Stid_Omit$geometry = NULL

# Get NA install years
missing_year <- Kruit_Stid_Omit[which(Kruit_Stid_Omit$Yr_inst<0), ]
# TEMP GET THE ALREADY CREATED YOD MANUAL INTERP SO I DONT HAVE TO DO AS MUCH WORK
yodvalid <- read.csv("S:/Users/stidjaco/R_files/Chapter_2/Data/Derived/Kruit_Stid_Omit/SAVE_kruitstidomit_completeYOD.csv")
missing_year$Yr_inst <- yodvalid$Yr_inst[match(missing_year$Index, as.numeric(yodvalid$Index))]

# Export to ArcPRO for YOD interp
#write.csv(missing_year, "Kruit_Stid_Omit_preYOD_08to18.csv")
rm(file_list, read_list, shp, coords)

# Post ArcPro YOD Validation
############################################################################################## Post Validation

#Call in validated yods
# 33238 was interesting array which shifted panel distribution between 2017 anbd 2018 and switched from double axis to fixed
NON_solar <- c(33116, 33179, 33400, 32140, 32229, 32469, 33061, 33262, 
               31891, 32853, 32929, 33063, 33170, 33313, 33424, 32541, 33262, 33370, 33339, 
               33134) # rooftop
#yodvalid <- read.csv("S:/Users/stidjaco/R_files/Chapter_2/Data/Derived/Kruit_Stid_Omit/Kruit_Stid_Omit_postYOD_08to18.csv")
#Kruit_Stid_Omit$Yr_inst <- ifelse(Kruit_Stid_Omit$Yr_inst<0, yodvalid$Yr_inst[match(Kruit_Stid_Omit$Index, as.numeric(yodvalid$Index))], Kruit_Stid_Omit$Yr_inst) 
Kruit_Stid_Omit <- Kruit_Stid_Omit[which(!Kruit_Stid_Omit$Index %in% NON_solar), ]
Kruit_Stid_Omit$Yr_inst <- ifelse(Kruit_Stid_Omit$Yr_inst<2008, 2008, 
                                  ifelse(Kruit_Stid_Omit$Yr_inst>2018, 2018, Kruit_Stid_Omit$Yr_inst)) # two of them were manually interp 2003 and 2007
rm(yodvalid)

# Test a manual removal based on capacity to panel area ratio
#Kruit_Stid_Omit$PArat_rm <- Kruit_Stid_Omit$pnl_a / (Kruit_Stid_Omit$TPVPp*1000)
#test_if_comission <- Kruit_Stid_Omit[which(Kruit_Stid_Omit$PArat_rm<1), ]
#write.csv(test_if_comission, "panel_area_capacity_comissiontest.csv") # Add these comissions to non_solar and re-run

# After checking if all YOD and Pnl_a present, and removing comissions, Now export to Sid's Model

# Export for Sid
#write.csv(Kruit_Stid_Omit, "Stidetal_OmitKruit_PV_ID.csv")

# Call in omitted dataset and our dataset
Kruit_Stid_Omit <- read.csv('S:/Users/stidjaco/R_files/Chapter_2/Data/Derived/Kruit_Stid_Omit/pvlib/Stid_KruitOMIT_wValidDualAxis_results.csv')
shp <- st_read('S:/Users/stidjaco/R_files/PhD_Ch1/Data/Derived/CONUS_Renewables_PreGEE/CONUS_Renewables.shp')
shp <- shp[which(shp$Index %in% Kruit_Stid_Omit$Index), ] %>% select(matches("Index"))
Kruit_Stid_Omit <- merge(shp, Kruit_Stid_Omit, by=c('Index'))
rm(shp)

# Also write a shp file for GEE CDL 5-year landcover classes
# Group with shp dataframe, get centroids, remove geometry for export for electricity model
#shp <- st_read('S:/Users/stidjaco/R_files/PhD_Ch1/Data/Derived/CONUS_Renewables_PreGEE/CONUS_Renewables.shp')
#shp <- shp[which(shp$Index %in% Kruit_Stid_Omit$Index), ] %>% select(matches("Index"))
#Kruit_Stid_Omit <- merge(shp, Kruit_Stid_Omit, by=c('Index'))
st_write(Kruit_Stid_Omit, 'Stidetal_OmitKruit_PV_ID_cropRot.shp')

##### FINAL CALL IN ########

# PREPARE KRUITWAGEN OMITTED ARRAYS: Electricity generation adjustments on kruitwagen omitted arrays -- Post pvlib modeling
# Call in generation model resutls (kWh)
generation_kWh <- read.csv('S:/Users/stidjaco/R_files/Chapter_2/Data/Derived/Kruit_Stid_Omit/pvlib_prep_and_results/Stid_KruitOMIT_wValidDualAxis_results.csv')
generation_kWh <- select(generation_kWh, matches("PV.generation|Index"))

# Convert to MW for easier export
generation_kWh[,2:12] <- generation_kWh[,2:12]/1000
is.num <- sapply(generation_kWh, is.numeric)
generation_kWh[,is.num] <- lapply(generation_kWh[,is.num], round, 4)
names(generation_kWh) <- c("Index", "Gen_2008", "Gen_2009", "Gen_2010", "Gen_2011", "Gen_2012", "Gen_2013", "Gen_2014", 
                           "Gen_2015", "Gen_2016", "Gen_2017", "Gen_2018")

# Call in omitted dataset again to retain attributes
Kruit_Stid_Omit <- read.csv('S:/Users/stidjaco/R_files/Chapter_2/Data/Derived/Kruit_Stid_Omit/pvlib_prep_and_results/Stid_KruitOMIT_wValidDualAxis_results.csv')
Kruit_Stid_Omit <- select(Kruit_Stid_Omit, matches("Index|Yr_inst|pnl_a|Class|TPVPp"))
Kruit_Stid_Omit <- merge(Kruit_Stid_Omit, generation_kWh, by=c('Index'))
Kruit_Stid_Omit$Index <- Kruit_Stid_Omit$Index %>% as.character() %>% as.integer()

# Get Kruit shapes as well to add
shp <- st_read('S:/Users/stidjaco/R_files/PhD_Ch1/Data/Derived/CONUS_Renewables_PreGEE/CONUS_Renewables.shp')
shp <- shp[which(shp$Index %in% Kruit_Stid_Omit$Index), ] %>% select(matches("Index|dir_a|Source"))
Kruit_Stid_Omit <- merge(shp, Kruit_Stid_Omit, by=c('Index'))
rm(generation_kWh, shp, is.num)

################ Incorporate All losses into Sids Results

# Adjust electricity generation to account for all loss (preinverter derate and losses from inverter efficiency)
# Vector of inverter derate % from NREL 2020 Cost Benchmark
soiling_loss <- 3 # %

# Inverter derate
inv_derate <- 100 - c(90, 90, 90, 90.1, 90.2, 90.3, 90.4, 90.5, 90.5, 90.5, 90.5)

# Inverter efficiency
inv_eff <- 100 - c(94, 94, 94, 94.8, 95.6, 96.4, 97.2, 98, 98, 98, 98)

# Generation efficiency accounting for all losses - soiling
gen_eff <- 1 - ((inv_derate + inv_eff - soiling_loss) / 100)
rm(soiling_loss, inv_derate, inv_eff)

# Get all generations
gen_df <- Kruit_Stid_Omit[, grepl("Gen|Yr_inst", names(Kruit_Stid_Omit))]
gen_df$geometry = NULL

# Multiply all generations by gen_eff for each respective year of installation
for(j in c(2:ncol(gen_df))){
  for(i in c(1:nrow(gen_df))){
    if(gen_df[i,]$Yr_inst==2008){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[1]
    }
    if(gen_df[i,]$Yr_inst==2009){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[2]
    }
    if(gen_df[i,]$Yr_inst==2010){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[3]
    }
    if(gen_df[i,]$Yr_inst==2011){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[4]
    }
    if(gen_df[i,]$Yr_inst==2012){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[5]
    }
    if(gen_df[i,]$Yr_inst==2013){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[6]
    }
    if(gen_df[i,]$Yr_inst==2014){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[7]
    }
    if(gen_df[i,]$Yr_inst==2015){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[8]
    }
    if(gen_df[i,]$Yr_inst==2016){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[9]
    }
    if(gen_df[i,]$Yr_inst==2017){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[10]
    }
    if(gen_df[i,]$Yr_inst==2018){
      gen_df[i,j] <- gen_df[i,j] * gen_eff[11]
    }
  }
}

# Merge dataframes
Kruit_Stid_Omit <- Kruit_Stid_Omit[, !grepl("Gen|Yr_inst", names(Kruit_Stid_Omit))]
Kruit_Stid_Omit <- merge.data.frame(Kruit_Stid_Omit, gen_df, by="row.names")
Kruit_Stid_Omit$Row.names=NULL
rm(gen_eff, i, j, gen_df)

# ----------------------------------------------------------------------------------------- 

# Post GEE CDL 5-year prior assessment to get crop rotation for Kruitwagen arrays
# Read in post gee analysis, group by index
crop_rotation <- read.csv("S:/Users/stidjaco/R_files/Chapter_2/Data/Derived/Kruit_Stid_Omit/Crop_Rotation/kruit_stid_yrPriorcrops.csv")
crop_rotation <- crop_rotation %>% group_by(Index) %>% summarise(Yr_inst = first(Yr_inst), sys_C_1 = cdl_yr[which(yr_prior==1)], sys_C_2 = cdl_yr[which(yr_prior==2)], sys_C_3 = cdl_yr[which(yr_prior==3)], sys_C_4 = cdl_yr[which(yr_prior==4)], sys_C_5 = cdl_yr[which(yr_prior==5)])

# ----------------------------------------------------------------------------------------- 

# Combine into final DF with Stid et al and export
# Get Kruit crop rotation
crop_rotation$Yr_inst = NULL
Kruit_Stid_Omit <- merge(crop_rotation, Kruit_Stid_Omit, by=c('Index'))
Kruit_Stid_Omit <- Kruit_Stid_Omit %>% st_as_sf()
rm(crop_rotation)

# Add nass crop types
# Get crop names from nass
nass_temp <- read.csv('S://Users/stidjaco/R_files/Chapter_2/Data/Derived/Nass_Classifications.csv')
in_df = Kruit_Stid_Omit
in_df$Crp_n_1 <- nass_temp$Crop[match(in_df$sys_C_1, nass_temp$Value)] 
in_df$Crp_n_2 <- nass_temp$Crop[match(in_df$sys_C_2, nass_temp$Value)] 
in_df$Crp_n_3 <- nass_temp$Crop[match(in_df$sys_C_3, nass_temp$Value)]
in_df$Crp_n_4 <- nass_temp$Crop[match(in_df$sys_C_4, nass_temp$Value)] 
in_df$Crp_n_5 <- nass_temp$Crop[match(in_df$sys_C_5, nass_temp$Value)] 
Kruit_Stid_Omit = in_df
rm(in_df, nass_temp)

# Call in and prep our dataset for co-located arrays
ungroupPV_df <- st_read('S:/Users/stidjaco/R_files/Chapter_2/Data/Downloaded/PV_Products/PV_ID_CV.shp')

# Grab ag-co-located arrays from our dataset
'%!in%' <- function(x,y)!('%in%'(x,y))

# Remove non croptypes
rm_non_ag <- function(df){
  # Subset by non ag year prior
  #df <- df[which(!df$Pnl_cnt==1), ]
  non_ag <- c(61, 63, 64, 65, 81, 82, 83, 87, 88, 92, 111, 112, 121, 122, 123, 124, 131, 141, 142, 143, 152, 190, 195)
  df <- df[df$sys_C_1 %!in% non_ag, ]
  
  # Remove all non ag
  df$Crp_n_1 <- ifelse(df$sys_C_1 %in% non_ag, NA, as.character(df$Crp_n_1))
  df$Crp_n_2 <- ifelse(df$sys_C_2 %in% non_ag, NA, as.character(df$Crp_n_2))
  df$Crp_n_3 <- ifelse(df$sys_C_3 %in% non_ag, NA, as.character(df$Crp_n_3))
  df$Crp_n_4 <- ifelse(df$sys_C_4 %in% non_ag, NA, as.character(df$Crp_n_4))
  df$Crp_n_5 <- ifelse(df$sys_C_5 %in% non_ag, NA, as.character(df$Crp_n_5))
  return(df)
}

# Apply to both ungroupPV_df and util data
crop_df <- rm_non_ag(ungroupPV_df)
Kruit_Stid_Omit <- rm_non_ag(Kruit_Stid_Omit)
rm('%!in%', rm_non_ag, ungroupPV_df)

# Get columns of interest and normalize
crop_df <- select(crop_df, matches("Index|Yr_inst|Class|Tot_a|Pnl_a|Gen|sys_c_|TPVPp|Crp_n_"))
crop_df$Source = "Stid"
crop_df$Class <- ifelse(crop_df$Class=="Si_Single_E/W", "single_axis",
                        ifelse(crop_df$Class=="Si_Fixed_S", "fixed_axis", NA))
names(crop_df)[which(names(crop_df) == "Pnl_a")] = "pnl_a"
names(crop_df)[which(names(crop_df) == "Tot_a")] = "dir_a"

# Merge df's 
crop_df_final <- rbind(crop_df, Kruit_Stid_Omit)
names(crop_df_final)[which(names(crop_df_final) == "TPVPp")] = "Capacity"
st_write(crop_df_final, 'CV_CoLocated_Solar_df.shp')
rm(list=ls())


# Run lanID irrig export from GEE, call in and group and make final export
pre_irrig <- st_read('S:/Users/stidjaco/R_files/Chapter_2/Data/Derived/PV_ID_Complete/pre_irrig/CV_CoLocated_Solar_df.shp')
irrig_df <- read.csv('S:/Users/stidjaco/R_files/Chapter_2/Data/Derived/PV_ID_Complete/pre_irrig/PV_ID_CV_irrig.csv')
post_irrig <- merge(pre_irrig, irrig_df, by='Index')
names(post_irrig)[which(names(post_irrig) == "Yr_inst")] = "Year"


# Ensure all columns are correct type
post_irrig$Index <- post_irrig$Index %>% as.character() %>% as.integer()
post_irrig$Class <- post_irrig$Class %>% as.character()
post_irrig$dir_a <- post_irrig$dir_a %>% as.character() %>% as.numeric()
post_irrig$pnl_a <- post_irrig$pnl_a %>% as.character() %>% as.numeric()
post_irrig$sys_C_1 <- post_irrig$sys_C_1 %>% as.character() %>% as.integer()
post_irrig$sys_C_2 <- post_irrig$sys_C_2 %>% as.character() %>% as.integer()
post_irrig$sys_C_3 <- post_irrig$sys_C_3 %>% as.character() %>% as.integer()
post_irrig$sys_C_4 <- post_irrig$sys_C_4 %>% as.character() %>% as.integer()
post_irrig$sys_C_5 <- post_irrig$sys_C_5 %>% as.character() %>% as.integer()
post_irrig$Crp_n_1 <- post_irrig$Crp_n_1 %>% as.character()
post_irrig$Crp_n_2 <- post_irrig$Crp_n_2 %>% as.character()
post_irrig$Crp_n_3 <- post_irrig$Crp_n_3 %>% as.character()
post_irrig$Crp_n_4 <- post_irrig$Crp_n_4 %>% as.character()
post_irrig$Crp_n_5 <- post_irrig$Crp_n_5 %>% as.character()
post_irrig$Capacity <- post_irrig$Capacity %>% as.character() %>% as.numeric()
post_irrig$Year <- post_irrig$Year %>% as.character() %>% as.integer()
post_irrig$irrig <- post_irrig$irrig %>% as.character() %>% as.integer()
post_irrig$Source <- post_irrig$Source %>% as.character()
post_irrig$Gen_2008 <- post_irrig$Gen_2008 %>% as.character() %>% as.numeric()
post_irrig$Gen_2009 <- post_irrig$Gen_2009 %>% as.character() %>% as.numeric()
post_irrig$Gen_2010 <- post_irrig$Gen_2010 %>% as.character() %>% as.numeric()
post_irrig$Gen_2011 <- post_irrig$Gen_2011 %>% as.character() %>% as.numeric()
post_irrig$Gen_2012 <- post_irrig$Gen_2012 %>% as.character() %>% as.numeric()
post_irrig$Gen_2013 <- post_irrig$Gen_2013 %>% as.character() %>% as.numeric()
post_irrig$Gen_2014 <- post_irrig$Gen_2014 %>% as.character() %>% as.numeric()
post_irrig$Gen_2015 <- post_irrig$Gen_2015 %>% as.character() %>% as.numeric()
post_irrig$Gen_2016 <- post_irrig$Gen_2016 %>% as.character() %>% as.numeric()
post_irrig$Gen_2017 <- post_irrig$Gen_2017 %>% as.character() %>% as.numeric()
post_irrig$Gen_2018 <- post_irrig$Gen_2018 %>% as.character() %>% as.numeric()
post_irrig <- post_irrig %>% st_as_sf()
st_write(post_irrig, 'CV_CoLocated_Solar_df.shp')
rm(list=ls())



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get croptypes from LC_RP_MultiYear (redone [check bottom of doc] because of maxRaw error in GEE-- see LC_RP_MultiYear.js)

# Remove all prior crop type information (only need shape and )
in_solar_df_temp <- st_read('Data/Derived/PV_ID_Complete/CV_CoLocated_Solar_df.shp', stringsAsFactors = FALSE, quiet = TRUE) %>% st_as_sf()
in_df <- in_solar_df_temp %>% dplyr::select(matches("Index|Year"))
#st_write(in_df_temp, "CV_Solar_LCprep_rnd2.shp")
# Post GEE CDL 5-year prior assessment to get crop rotation for Kruitwagen arrays
# Read in post gee analysis, group by index
crop_rotation <- read.csv("Data/Derived/Kruit_Stid_Omit/Crop_Rotation/kruit_stid_yrPriorcrops2.csv")
crop_rotation <- crop_rotation %>% group_by(Index) %>% summarise(Year = first(Year), sys_C_1 = cdl_yr[which(yr_prior==1)], sys_C_2 = cdl_yr[which(yr_prior==2)], sys_C_3 = cdl_yr[which(yr_prior==3)], sys_C_4 = cdl_yr[which(yr_prior==4)], sys_C_5 = cdl_yr[which(yr_prior==5)])
# Add nass crop types
in_df <- crop_rotation
# Get crop names from nass
nass_temp <- read.csv('S://Users/stidjaco/R_files/Chapter_2/Data/Derived/Nass_Classifications.csv')
in_df$Crp_n_1 <- nass_temp$Crop[match(in_df$sys_C_1, nass_temp$Value)] 
in_df$Crp_n_2 <- nass_temp$Crop[match(in_df$sys_C_2, nass_temp$Value)] 
in_df$Crp_n_3 <- nass_temp$Crop[match(in_df$sys_C_3, nass_temp$Value)]
in_df$Crp_n_4 <- nass_temp$Crop[match(in_df$sys_C_4, nass_temp$Value)] 
in_df$Crp_n_5 <- nass_temp$Crop[match(in_df$sys_C_5, nass_temp$Value)] 
# Grab ag-co-located arrays from our dataset
'%!in%' <- function(x,y)!('%in%'(x,y))
# Remove non croptypes
rm_non_ag <- function(df){
  # Subset by non ag year prior
  #df <- df[which(!df$Pnl_cnt==1), ]
  non_ag <- c(61, 63, 64, 65, 81, 82, 83, 87, 88, 92, 111, 112, 121, 122, 123, 124, 131, 141, 142, 143, 152, 190, 195)
  df <- df[df$sys_C_1 %!in% non_ag, ]
  # Remove all non ag
  df$Crp_n_1 <- ifelse(df$sys_C_1 %in% non_ag, NA, as.character(df$Crp_n_1))
  df$Crp_n_2 <- ifelse(df$sys_C_2 %in% non_ag, NA, as.character(df$Crp_n_2))
  df$Crp_n_3 <- ifelse(df$sys_C_3 %in% non_ag, NA, as.character(df$Crp_n_3))
  df$Crp_n_4 <- ifelse(df$sys_C_4 %in% non_ag, NA, as.character(df$Crp_n_4))
  df$Crp_n_5 <- ifelse(df$sys_C_5 %in% non_ag, NA, as.character(df$Crp_n_5))
  return(df)}
# Apply to both ungroupPV_df and util data
crop_df <- rm_non_ag(in_df)
# Make sure newly omitted arrays are removed
in_solar_df_temp <- in_solar_df_temp[which(in_solar_df_temp$Index %in% crop_df$Index), ]
# Match new crop types
in_solar_df_temp$sys_C_1 <- crop_df$sys_C_1[match(in_solar_df_temp$Index, crop_df$Index)] 
in_solar_df_temp$sys_C_2 <- crop_df$sys_C_2[match(in_solar_df_temp$Index, crop_df$Index)] 
in_solar_df_temp$sys_C_3 <- crop_df$sys_C_3[match(in_solar_df_temp$Index, crop_df$Index)] 
in_solar_df_temp$sys_C_4 <- crop_df$sys_C_4[match(in_solar_df_temp$Index, crop_df$Index)] 
in_solar_df_temp$sys_C_5 <- crop_df$sys_C_5[match(in_solar_df_temp$Index, crop_df$Index)] 
in_solar_df_temp$Crp_n_1 <- crop_df$Crp_n_1[match(in_solar_df_temp$Index, crop_df$Index)]
in_solar_df_temp$Crp_n_2 <- crop_df$Crp_n_2[match(in_solar_df_temp$Index, crop_df$Index)]
in_solar_df_temp$Crp_n_3 <- crop_df$Crp_n_3[match(in_solar_df_temp$Index, crop_df$Index)]
in_solar_df_temp$Crp_n_4 <- crop_df$Crp_n_4[match(in_solar_df_temp$Index, crop_df$Index)]
in_solar_df_temp$Crp_n_5 <- crop_df$Crp_n_5[match(in_solar_df_temp$Index, crop_df$Index)]
st_write(in_solar_df_temp, "CV_Agrisolar_df.shp")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Check which arrays removed from Stid et al., 2022

# From the above code for "Get croptypes from LC_RP_MultiYear" (re-do crop rotation analysis to fix CDL error)
in_df_missingStid <- in_df[which(!in_df$Index %in% crop_df$Index), ]
insolartemp <- in_solar_df[which(!in_solar_df$Index %in% crop_df$Index), ]
croprottemp <- crop_rotation[which(!crop_rotation$Index %in% crop_df$Index), ]
# Three are removed because of maxRaw issue, and 7 are remvoed because of condition that active ag land had to be the year prior to installation, not just in crop history

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get total loss dataframe for Sid to use in electricity rate analysis

# Get all losses 
soiling_loss <- 3 # %
# Inverter derate
inv_derate <- 100 - c(90, 90, 90, 90.1, 90.2, 90.3, 90.4, 90.5, 90.5, 90.5, 90.5)
# Inverter efficiency
inv_eff <- 100 - c(94, 94, 94, 94.8, 95.6, 96.4, 97.2, 98, 98, 98, 98)

# Generation efficiency accounting for all losses - soiling
eff_loss <- ((inv_derate + inv_eff + soiling_loss))
df_export <- data.frame(Year = c(2008:2018), derating_factor = eff_loss)
write.csv(df_export, "Data/Derived/Derating_Loss_Total_08to18.csv", row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pull in historic utility dataset from Sids model and merge with dataset

# Call in cv dataset - Delete current CV_Agrisolar_df first
in_solar_df <- st_read('All_Backup/CV_Agrisolar_df_backup/CV_Agrisolar_df.shp', stringsAsFactors = FALSE, quiet = TRUE) %>% st_as_sf()

# Call in raw data output from AgrisolarNEM_Returns.py
# For each install year dataframe, columns are class (1 or 0 for fixed or single respectively), Index, Econ_20## for install year, Econ_20## for install year +1 and so on
file_list <- list.files("pvlib_Analysis/FEWLS_pvlibResults", pattern = "*csv", full.names = TRUE)
for(file in file_list){
  file_name <- gsub('.csv','', file)
  file_name <- gsub(".*FEWLS_pvlibResults/", "", file_name)
  assign(file_name, read.csv(file, header = FALSE)) # Header = false is important because econ model outputs have no columns, so read.csv uses the first row
}

# For each install year file, remove NULL, and normalize so all files have same columns, and rbind into one df
names <- list.files("pvlib_Analysis/FEWLS_pvlibResults", pattern = "*csv", full.names = FALSE)
names <- gsub('.csv','', names)
start_year = 2008
end_year = 2018
df <- data.frame(matrix(nrow = 0, ncol = (length(c(start_year:end_year))+1)))
for(name in names){
  df_temp = get(name) # gets data from global environment
  yr_inst <- gsub(".*economic_value_", "", name) %>% as.numeric() # get installation year for dataset
  df_temp <- df_temp[which(df_temp[,2]>0), ] # removes empty export matrix resulting from AgrisolarNEM_Returns.py methods 
  econ_colNames <- paste("Econ_", c(yr_inst:end_year) %>% as.character(), sep = "") # get column names to something readable
  names(df_temp) <- c("Class_bin", "Index", econ_colNames)
  missing_econ_colNames <- c(start_year:end_year)[which(!c(start_year:end_year) %in% c(yr_inst:end_year))] # get missing model years because prior to installation
  if( length(missing_econ_colNames)>0 ){ # if model years exist for years prior to current df installation year, add empty columns
    df_temp[,c( (ncol(df_temp)+1) : (ncol(df_temp)+length(missing_econ_colNames)))] <- as.numeric(0)
    names(df_temp) <- c("Class_bin", "Index", econ_colNames, paste("Econ_", missing_econ_colNames %>% as.character(), sep = ""))
  }
  df_temp$Class_bin = NULL # at the moment, unneccessary column
  df <- rbind(df, df_temp)
}

# Fill df with NA for non-commercial scale installations
util_rows = in_solar_df[which(!in_solar_df$Index %in% df$Index), ] %>% nrow()
df_util <- data.frame(matrix(nrow = util_rows, ncol = (length(c(start_year:end_year))+1)))
names(df_util) <- c("Index",  paste("Econ_", c(start_year:end_year), sep = ""))
df_util$Index = in_solar_df$Index[which(!in_solar_df$Index %in% df$Index)]
df <- rbind(df, df_util)

# Call in and merge
in_solar_df <- merge(in_solar_df, df, by = "Index")
st_write(in_solar_df, "CV_Agrisolar_df.shp")

#names(economic_value_2018) <- c("X1.0", "X2.0", "X3.0")
#economic_value_2018 <- add_row(economic_value_2018, X1.0 = 1, X2.0 = 2, X3.0 = 60275.4600469082)
#economic_value_2018 <- economic_value_2018[-c(64), ]
#write.csv(economic_value_2018, "Econ_Analysis/FEWLS_EconResults/economic_value_2018.csv", row.names = FALSE)

##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// TEMPORARY
"
# Get crop unique df
crops <- unique(c(in_solar_df$Crp_n_1, in_solar_df$Crp_n_2, in_solar_df$Crp_n_3, in_solar_df$Crp_n_4, in_solar_df$Crp_n_5))
in_df <- in_solar_df
in_df$geometry = NULL
crp1_index <- in_df %>% group_by(Crp_n_1) %>% summarise(Index = first(Index))
crops_miss <- setdiff(crops, crp1_index$Crp_n_1)
crp2_index <- in_df[which(in_df$Crp_n_2 %in% crops_miss), ] %>% group_by(Crp_n_2) %>% summarise(Index = first(Index))
crops_miss <- setdiff(crops_miss, crp2_index$Crp_n_2)
crp3_index <- in_df[which(in_df$Crp_n_3 %in% crops_miss), ] %>% group_by(Crp_n_3) %>% summarise(Index = first(Index))
crops_miss <- setdiff(crops_miss, crp3_index$Crp_n_3)
crp4_index <- in_df[which(in_df$Crp_n_4 %in% crops_miss), ] %>% group_by(Crp_n_4) %>% summarise(Index = first(Index)) # covered all crops
crop_indexes <- c(crp1_index$Index, crp2_index$Index, crp3_index$Index, crp4_index$Index)
in_solar_df <- in_solar_df[which(in_solar_df$Index %in% crop_indexes), ]
"
"THIS WILL BE DONE IN SUPPELEMENTAL CODE ONCE SIDS ANALYSIS IS DONE SO THAT INPUT DATASET CONTAINS ECON VARIABLES"
##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// TEMPORARY

