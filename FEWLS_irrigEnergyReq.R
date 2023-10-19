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
This is an optional file (downloaed folder comes with defaults) that generates a dataframe
of irrigation energy requirements (GWh/m3) at the county level from Methods derived in 
McCarthy et al. (2021; 2022). Only needs to be run if cleanFEWLS() is run, or if new DTW
product is available. 
"
## ----------------------------------------------------------- ##
##                                                             ##
##       Build irrigation energy requirement estimation        ##
##                                                             ##
## ----------------------------------------------------------- ##

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# Packages for Processing 
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(sf)
library(lwgeom)
library(tidyverse)
library(base)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(miscTools)
library(lubridate)
library(tigris)
library(stringr)
library(smoothr)

# METHOD
# IrrigEnergyReq (kWh) = Mass(m3/1000 kg_per_m3) * Lift (m) * g (9.8016 m/s^2) * E (77% pump efficiency)
# Lift (m) = Lwt (distance to GW) + Lpr (lift distance due to pressurization of water and pipe losses) + Lcd&wd (lift from cone of depression and drawdown from well efficiency)
# Lcd&wd = Q (pumprate) / 4p*pi*Transmissivity [-0.5772 - ln(r2^S / 4*Transmissivity*timepumped)] * 1.5 (accoutngin fro well drag of well efficinty 50%)

# ------------------- #
#  Neccessary Inputs  #
# ------------------- #
state = "California"
statefp = 6
county_irrig_energy_vars_loc = 'Data/Downloaded/County_IrrigEnergyInputsMcCarthy.csv'
irrig_energy_vars <- read.csv(county_irrig_energy_vars_loc, stringsAsFactors = FALSE)

# ------------------- #
# New SGMA County DTW # -- From: https://sgma.water.ca.gov/webgis/?appid=SGMADataViewer#gwlevels
# ------------------- #

# For now, specifically Central Valley, but give option for CONUS
dtw_useCONUS = FALSE

# Call in depth to water raster -- if California, use SGMA product, else use CONUS product from Zell and Sanford
if(dtw_useCONUS == FALSE){
  dtw <- raster("Data/Derived/DTW_raster/Spring2018_DTW.tif") * 0.3048 # Convert feet to meters
} else {
  dtw <- raster("Data/Derived/Depth_To_Water/Zell_Sanford_2020/conus_MF6_SS_Unconfined_250_dtw.tif") 
}

# Pull in counties shape file
counties <- counties(state, TRUE) %>% as("Spatial") %>% spTransform(crs(dtw)) %>% st_as_sf()
counties$COUNTYFP <- gsub("(?<![0-9])0+", "", counties$COUNTYFP, perl = TRUE)

# Only use county shapes and regions which intersect CV
CV_Boundary <- st_read("Data/Downloaded/CV_alluvial_bnd/Alluvial_Bnd.shp") %>% as("Spatial") %>% spTransform(crs(dtw))
CV_Boundary <- fill_holes(CV_Boundary, threshold = units::set_units(50000000, 'km^2')) %>% st_as_sf()

# Get intersection
cv_counties <- st_intersection(counties, CV_Boundary)
cv_counties <- cv_counties %>% dplyr::select(matches("COUNTYFP"))

# Extract new depth to water from counties
cv_counties$dtw <- raster::extract(dtw, cv_counties, fun=mean) %>% as.numeric()

# ------------------------------------------------------ #
# Get Zell and Sanford dtw from bens results and compare # 
# ------------------------------------------------------ #

# MS Ben Mcarthy's County irrigation energy use dataframe results
irrigEnergyVars <- irrig_energy_vars[which(irrig_energy_vars$STATE_FIPS==statefp), ]
irrigEnergyVars$FIPS <- irrigEnergyVars$FIPS %>% as.character()
irrigEnergyVars$COUNTYFP <- substring(irrigEnergyVars$FIPS, 2)
irrigEnergyVars$COUNTYFP <- gsub("(?<![0-9])0+", "", irrigEnergyVars$COUNTYFP, perl = TRUE)
irrigEnergyVars <- irrigEnergyVars[which(irrigEnergyVars$COUNTYFP %in% cv_counties$COUNTYFP), ]

# Check for sanity depth to water reported vs new product
dtw_zands <- data.frame(COUNTYFP = as.numeric(irrigEnergyVars$COUNTYFP), dtw = irrigEnergyVars$DTW)
dtw_sgma <- data.frame(COUNTYFP = as.numeric(cv_counties$COUNTYFP), dtw = cv_counties$dtw)
dtw <- merge(dtw_sgma, dtw_zands, by=c("COUNTYFP"))
names(dtw) <- c("COUNTYFP", "SGMA_dtw", "ZandS_dtw")
dtw$diff <- dtw$SGMA_dtw - dtw$ZandS_dtw

# ------------------------------------------------- #
# Recalculate irrig energy requirement with new dtw # 
# ------------------------------------------------- #

# Get new dtw 
irrigEnergyVars$DTW <- dtw_sgma$dtw[match(irrigEnergyVars$COUNTYFP, dtw_sgma$COUNTYFP)]

# Calculate new lift
pump_loss = 10 # 10psi
pump_eff = 0.77 # pumping efficiency
well_radius = 0.1524 # 6 inches
psi_m_of_head = 1/1.41 # Pressure per m of head
well_eff_correction = 1.5 # Assuming a 50% well efficiency
irrigEnergyVars$psi_lift <- (irrigEnergyVars$PSI+pump_loss) * psi_m_of_head
irrigEnergyVars$Liftcdwd <- (irrigEnergyVars$PUMP_RATE / (4 * pi * irrigEnergyVars$TRANS)) * (-0.5772 - log((well_radius^2 * irrigEnergyVars$spec_yield_lay3) / (4 * irrigEnergyVars$TRANS * irrigEnergyVars$irrigDays)))
irrigEnergyVars$Total_Lift <- (irrigEnergyVars$DTW + irrigEnergyVars$psi_lift + irrigEnergyVars$Liftcdwd) * well_eff_correction

# Calcualte irrig energy based on new lift
kg_m3 = 1000
g = 9.8016 # m/s^2
irrigEnergyVars$Total_Energy <- (irrigEnergyVars$GW_Total*kg_m3 * irrigEnergyVars$Total_Lift * g) + irrigEnergyVars$SW_Energy

# Export with GWh_m3
joules_GWh = 3.6e12 # 3.6e12 joules per GWh
irrigEnergyVars$GWh_m3 <- (irrigEnergyVars$Total_Energy/joules_GWh) / (irrigEnergyVars$GW_Total + irrigEnergyVars$SW_Total)
write.csv(irrigEnergyVars, "Data/Derived/irrigEnergyReq_adjusted.csv", row.names=FALSE)