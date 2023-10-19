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
These are the necesssary input variables for the FEWLS tool. Specified are several essential
dataframes (delineated with loc) as well as conversion factors and constant variables. When
chaning dataframes, ensure attributes contain the same nomenclature and units. When chaning 
constant variables, ensure constant units. This is, for instance, where the solar PV 
lifespan can be altered from 25-years, to 15-years, or to 50-years to acquire new resource 
and economic end-values under different scenarios. Final exported dataframes and figures
contain information on the lifespan used and the capacity threshold used to delineate 
utility- vs. commercial-scale practices. The FEWLS_model.R script automatically calls and 
compiles these inputs. At the moment, this model is limited to installation years between 
2008 and 2018 and for California's Central Valley. Although, the only limiting resource is 
water use and irrigation maps.
"

# Set global settings (warnings to false, figs to real dpi)
options(warn = -1)
#showtext_auto()
#showtext_opts(dpi = 300)

# Input file locations: 
# Spatio-temporal solar dataset with installation years, removed crop types 5-years prior, capacity, presence of irrigation year prior, and historical 
# utility data. Example with required attributes and names found in "Data/Example" 
historical_solar_price_df_loc = "Data/Derived/final_PV_ID_CV_utility_results.csv"
inflation_ratesPPI_loc = 'Data/Downloaded/PPI.csv' # Not CPI
inflation_ratesCPI_loc = 'Data/Downloaded/CPI.csv' # Not CPI
county_irrig_energy_vars_loc = 'Data/Derived/irrigEnergyReq_adjusted.csv'
NASS_classes_loc = 'Data/Derived/Nass_Classifications.csv'
food_yield_loc = 'Data/Downloaded/USDA_yield_08to18_state.csv'
FoodData_food_loc = 'Data/Downloaded/FoodData/food.csv'
FoodData_nutrient_loc = 'Data/Downloaded/FoodData/food_nutrient.csv'
food_dollarValue_loc = 'Data/Downloaded/crop_values_CA_05_18.csv'
energy_efficiency_df_loc = 'Data/Downloaded/efficiency_yearly_LBerkeley.csv'
USGS_wu_loc = 'Data/Downloaded/IrrigWU_CONUScounty85to15.csv'
conus_county_annu_precip_loc = 'Data/Derived/county_annu_gridMET.csv'
water_irrig_perCrop_state_loc = 'Data/Downloaded/USDA_FARS_1318.csv'
water_irrig_perCrop_conus_loc = 'Data/Downloaded/USDA_FARS_1318_CONUS.csv'
agland_operation_conus_loc = 'Data/Downloaded/CONUS_aglandarea_operations.csv'
NREL_ATB_future_loc = 'Data/Downloaded/NREL_ATB_historical_all.csv'
NREL_ATB_historical_loc = 'Data/Downloaded/NREL_ATB_Historical_Trends_data.csv'

# Constants and Conversions (conversion are written in the "from_to" format)
# FOOD 
m2_acre = 4046.86 # m2/acre
LB_kg = 1 / 2.205 # 2.205 lbs per kg 
BU_kg = 60 * LB_kg # 60lbs of wheat in one bushel * lbs per kg -- Wheat: https://www.cfd.coop/go/doc/f/CMDT_Tables_for_Weights_and_Measurement_for_Crops.pdf -- Eggplant: https://extension.missouri.edu/publications/g6369#:~:text=Eggplants%20can%20yield%20500%20to,pounds%20per%20bushel)%20per%20acre. -- 33lbs????
TONS_kg = 2000 * LB_kg # 2000 lbs per short ton / lbs per kg -- LOGIC
CWT_kg = 100 * LB_kg # 100 lbs per hundred weight / lbs per kg -- LOGIC
BOXES_kg = 55 * LB_kg # 55 lbs per box / lbs per kg -- Oranges: https://coststudyfiles.ucdavis.edu/uploads/cs_public/13/0d/130df154-60ed-4096-842a-3e1e64da87e7/orangevs2005.pdf 
BOX_kg = 80 * LB_kg # 80 lbs per box / lbs per kg -- CA oranges, grapefruit, lemons, tangerines, mandarines: https://www.nass.usda.gov/Publications/Todays_Reports/reports/cfrt0820.pdf
BARRELS_kg = 100 * LB_kg # 100 lbs per barrel / lbs per kg -- Cranberry: https://www.nass.usda.gov/Statistics_by_State/New_Jersey/Publications/Cranberry_Statistics/2014%20CRANSUM.pdf
acres_m2 = 4046.86 # m2/acre

# Energy
joules_GWh = 3.6e12 # 3.6e12 joules per GWh
watt_GW = 1000000000 # watts per gigawatt
kWh_GWh = 1000000 # kW per GW and kWh per GWh
kW_MW = 1000 # kW per MW
MW_GW = 1000 # MW per GW and MWh per GWh
kW_GW = 1000000 # kilawatts per gigawatt
GW_TW = 1000 # GW per TW

# WATER
mm_acft_acre = 304.8 # acrefeet per acre (depth in feet) to milimieters
mm3_year_Mgal_day = 365*1000000/325851.43188913*1233000000000 # Days per year * gal per mgal / gal per acrefoot * mm3 per acrefoot = mm3/yr
acft_m3 = 1233.4818375475 # acrefeet per m3
mm_m = 1000 # millimeters per 
m3_millionm3 = 1000000 # cubic meters per million cubic meters
gal_acft = 325851.42857143 # gallons per acre-foot

# Variable Inputs

# NOTES: All input vectors must be the same length as the model length. All economic values here are nominal. The model will adjust for inflation depending on the desired USD year. Also note variables which need to be temporal vectors based on system lifetime and start install year

#////////////////////////#
#________________________#
#    Setup Variables     #
#________________________#
#\\\\\\\\\\\\\\\\\\\\\\\\#
# REQUIRED INPUTS ____________________________________

# Model variables
system_lifespan = 25 # system lifespan in years
capacity_threshold = 5 # MW capacity threshold -- cutoff for Net Metering
use_single_elecRate = FALSE # Check 'Energy Variables' for this rate (elec_price). This is TRUE only for running sensitivity analysis, or if no economic modeling has been done (in which case this is average $/kWh expected returns)

# Nominal discount rate for Discount Cash Flow (DCF) assessment 
discount_rate_nominal = 0.05 # DECIMAL -- 5% nominal discount rate Kelley et al. 2010 (https://doi.org/10.1016/j.rser.2010.07.061) "On the feasibility of solar powered irrigation", many others with similar rates discussed in agrisolar manuscript
inflation_rate_future = 0.03 # DECIMAL -- 3% expected inflation rate, historical average (2000 to 2022) PPI and CPI (3.4% and 2.4%) and roughly the same as Vartiainen et al. (2020) and Liu et al. (2014) respectively

# Inflation adjustment year, non-agricultural CDL crops, and State variables, although, model is currently only fully operable in CA. 
USD_infl_adj_year = 2018 # desired inflation adjustment year (ie. if 2018, all resulting economic values are in 2018 USD)
cdl_non_ag = c(0, 61, 63, 64, 65, 81, 82, 83, 87, 88, 92, 111, 112, 121, 122, 123, 124, 131, 141, 142, 143, 152, 190, 195) # CDL ag classes to exclude in analysis
state = "California" # state which dataset is contained
statefp = 6 # state fp which dataset is contained


#////////////////////////#
#________________________#
#     Food Variables     #
#________________________#
#\\\\\\\\\\\\\\\\\\\\\\\\#
# Resource Impact Variables & Scenarios ____________________________________
# Installing solar onlow yield land
yield_deficit_base = 0.065 # 6.5% new cropland yield deficit, similar to our value 4% https://www.nature.com/articles/s41467-020-18045-z & stid et al., 2022 -- Fixed axis marginal land use
yield_deficit_best = 0.25 # 25% less productive than avg # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0089501, https://www.science.org/doi/10.1126/sciadv.aav9318
yield_deficit_worst = 0.0 # worst case scenario is if solar is installed on land with average yields (would likely never be higher than average yields)

# Potential changes in yield through time
yield_time_base = 0.0 # no change in average yield -- https://link.springer.com/article/10.1007/s10584-011-0305-4 -- some crops show no change in yield
yield_time_worst= 0.008 #  0.8%/yr increase in yield from Deepack et al., 2013 (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0066428)
yield_time_best = -0.04/41 # Lee et al., 2011 (https://link.springer.com/article/10.1007/s10584-011-0305-4) 4% decrease in wheat yield between 2009 and 2050, similar to Lobell et al., 2011 (https://link.springer.com/article/10.1007/s10584-011-0303-6) 5% decrease in grape yield between 2000 and 2010

# Crop Production Adjustment for Non-direct Human Consumption (Grazing and Silage) ____________________________________
feedsilage_cdl_crops = c("Alfalfa", "Corn", "Grassland/Pasture", "Clover/Wildflowers") # https://corn.ucanr.edu/About_California_Corn/ -- Most corn in CA is used for silage # https://ucanr.edu/sites/Tulare_County/files/74080.pdf -- slide 8 (top feed crops)
calorie_conversion_efficiency_beef = 0.029 # https://iopscience.iop.org/article/10.1088/1748-9326/11/10/105002 -- Table 1 (2.9%)
calorie_conversion_efficiency_dairy = 0.17 # https://iopscience.iop.org/article/10.1088/1748-9326/11/10/105002 -- Table 1 (17.0%)
regional_dairy_beef_ratio = 1740 / 670 # https://www.cdfa.ca.gov/statistics/PDFs/2018-2019AgReportnass.pdf ---- 1720000 / 680000 # https://www.nass.usda.gov/Quick_Stats/Ag_Overview/stateOverview.php?state=CALIFORNIA

# Crop Production Adjustment for Non-direct Human Consumption (Grazing and Silage) ____________________________________
seedoil_cdl_crops = c("Canola", "Sunflower", "Safflower", "Flaxseed", "Cotton")
seedoil_mass_prop = 0.61 # https://www.cdfa.ca.gov/statistics/PDFs/2018-2019AgReportnass.pdf -- 61% from pg 13
seedoil_feed_prop = 0.95 # https://ccgga.org/cotton-information/ca-cotton-facts/#:~:text=2011%20production%20is%20estimated%20in,is%20crushed%20for%20the%20oil. -- 95% used to feed dairy cattle


#////////////////////////#
#________________________#
#    Energy Variables    #
#________________________#
#\\\\\\\\\\\\\\\\\\\\\\\\#

# Resource Impact Variables & Scenarios ____________________________________
irrig_energymix_prop = 0.7 # https://agwaterstewards.org/practices/water_energy/#:~:text=Agriculture%20in%20California%20also%20consumes,Southern%20and%20Western%20Central%20Valley.

# Electricity price ($/kWh) -- Only when performing sensitivity analysis "getEconomic_energy_SingleRate()" function -- can also be used if electric utility rate information is less available to approximate NEM returns
# This can be manually set, but is calculated as the average inflation adjusted (2018 USD) rate of all utilities ($/kWh) in FEWLS_supplemental.R

# Only when performing sensitivity analysis "getEconomic_energy_SingleRate()" function -- can also be used if electric utility rate information is less available to approximate NEM returns
# Electricity price ($/kWh) used in the economic analysis for the offset irrigation-water use electricity requirements
elec_price = 0.1722 # $/kWh

# Degradation rates -- currently dependent on mount # New degradation rates based on median rates delineated by mount tech from Jordan et al., 2022 (https://doi.org/10.1002/pip.3566), updated from DC Jordan et al. prior work in 2016
efficiency_degradation_rate = 0.0075 # DECIMAL NOT PERCENT -- Although Jordan et al., 2022 reports -0.68%/yr median fixed-axis & 0.0075 # -0.76%/yr median single_axis, they report no statistic diff, with average of -0.75%/yr

# Economic Impact Variables & Scenarios ____________________________________
# End install year (2018) to end model year (2042) from AEO 2019 Table 8 "End-Use Prices: Commercial" 2018 USD (line 58) -- checked AEO 2023, rate of change nearly identical. Growth rate vector is therefore, an annual real % change. Direct from datatable line 58. 
real_EnergyPriceAEO = c(10.835227,10.747263,10.456611,10.315327,10.318168,10.31889,10.397223,10.520648,10.587178,10.585861,10.600105,10.575869,10.602732,10.608257,10.651743,10.661316,10.665332,10.656591,10.668245,10.67905,
                        10.65926,10.622561,10.556154,10.577965,10.536174)
real_EnergyPriceAEO = data.frame(n = c(1:length(real_EnergyPriceAEO)), r = real_EnergyPriceAEO) %>% mutate(pct_change = (r/lag(r) - 1)) # results in decimal percent change, not % -- ie, dont multiply by 100
# Rate of growth in energy prices -- code includes auto extension of this data incase it does not encompass full lifesopan (extends vector by vector average). 
# Also, does not matter initial real year USD adjusted to, since any inflation adjustment would just be a multiplier and we are calculating a rate of change. 
energy_price_growth_rate = real_EnergyPriceAEO$pct_change %>% replace(is.na(.), 0) %>% as.numeric() # gives 2018 NA rate of change a 0% rate of change, end result in % annual change IN DECIMAL


#////////////////////////#
#________________________#
#    Water Variables     #
#________________________#
#\\\\\\\\\\\\\\\\\\\\\\\\#
# Resource Impact Variables & Scenarios ____________________________________
lowOandM_wu = 0.23 * acft_m3 # low for dry cooled systems https://www.osti.gov/servlets/purl/1090206 - MW
highOandM_wu = 2.16 * acft_m3 # high for dry cooled systems https://www.osti.gov/servlets/purl/1090206 - MW
baseOandM_wu = 1.195 * acft_m3 # average of high and low https://www.osti.gov/servlets/purl/1090206 - MW

# Economic Impact Variables & Scenarios ____________________________________
waterRightDollar_m3 = 40 / acft_m3 # $ per acrefeet  / acrefeet per m3 = $/m3
waterRightDollar_year = 2018 # inflation adjustment


#//////////////////////////////#
#______________________________#
# Utility Land Lease Variables #
#______________________________#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# Economic Impact Variables & Scenarios ____________________________________
landlease_price_growth_rate = 0 # (rep(0, install_period), rep(0.02, model_length-install_period)) # Rates should be near estimates of future inflation with average figures between 1.5% to 2.5% annually -- https://strategicsolargroup.com/what-is-the-average-solar-farm-lease-rate/
landlease_per_m2_avg = 1000 / acres_m2 # https://www.vantrumpreport.com/what-you-need-to-know-big-money-leasing-farmland-to-solar-operators/
landlease_per_m2_min = 300  / acres_m2 # https://www.vantrumpreport.com/what-you-need-to-know-big-money-leasing-farmland-to-solar-operators/
landlease_per_m2_max = 2000 / acres_m2 # https://www.vantrumpreport.com/what-you-need-to-know-big-money-leasing-farmland-to-solar-operators/
# https://www.landmarkdividend.com/solar-farm-land-lease-rates-2/ -- good general land lease method source
# https://www.ysgsolar.com/blog/how-much-money-can-solar-farm-make-2022-ysg-solar -- $21,000 - $41,000 per acre???
# https://farmoffice.osu.edu/sites/aglaw/files/site-library/Farmland_Owner%27s_Guide_to_Solar_Leasing.pdf -- interesting discussion on developer paying land owner farmer
# https://www.blm.gov/sites/blm.gov/files/docs/2022-05/MS-2806%20rel%202-307%20Chapter%206.pdf -- table E & D; Land "rent" from BLM??? -- PUBLIC LAND $75/acre -- also gives out to 2050


#////////////////////////#
#________________________#
# O&M Economic variable  #
#________________________#
#\\\\\\\\\\\\\\\\\\\\\\\\#
# Economic Impact Variables & Scenarios ____________________________________
NREL_ATB_yr = 2022 # Year of the Annual Technology Baseline to adjust for inflation for O&M Costs


#/////////////////////////////#
#_____________________________#
# Installation Cost VARIABLES #
#_____________________________#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# Economic Impact Variables & Scenarios ____________________________________
# SET UP: Commercial CAPEX installation costs from NREL Annual Technology Baseline: https://atb.nrel.gov/electricity/2022/commercial_pv -- in 2019 USD, needs to be adjusted
inst_USD_yr = 2019

# Solar Investment Tax Credit (30% up through 2023)
solar_ITC <- c(0.30) # 30% solar ITC between 2008 and 2018


#/////////////////////////////#
#_____________________________#
#    Plot Input Variables     #
#_____________________________#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# Plot colors
food_col = "wheat4"
energy_col = "orange"
water_col = "blue"
CO2_col = "black"
phase_alpha = 0.2
crop_plot_cols = c("berry totals" = "purple3", 
                   "cotton" = "lightskyblue1", 
                   "crops, other" = "darkorange1", 
                   "grain" = "peru", 
                   "hay/pasture" = "khaki2",
                   "orchards" = "red2", 
                   "vegetables" = "yellowgreen")

# Border size for plots
border_size = 0.75

# Factor to scale kcal by
kcal_adj= 1e9 # Convert to billions of kcal

# Factor to scale USD by 
econ_scale = 1e9 # Billions of 2018 USD

# Plotting crop groups
crp_grp1 = "wheat|sorghum|grain|rice"
crp_grp2 = "vegetable|tomatoes|potatoes|beans|corn"
crp_grp3 = "hay|pastureland"
crp_nme1 = "grain"
crp_nme2 = "vegetables"
crp_nme3 = "hay/pasture"

# Set margins
legend_margin_plot <- unit(c(0,-0.5,0,-0.5),"cm")
margin_plot <- unit(c(0,0,0,0),"cm")
right_margin_plot <- unit(c(0,-0.5,0,0),"cm")
solo_margin_plot <- unit(c(0.5, 0.5, 0.5, 0.5),"cm")

# Old Color Profile
#crop_plot_cols = c("berry totals" = "red2", "cotton" = "lightskyblue1", "crops, other" = "yellowgreen", "grain" = "peru", "hay/pasture" = "khaki2","orchards" = "purple3", "vegetables" = "darkorange1")
