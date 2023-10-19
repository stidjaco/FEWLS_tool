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
This is the execute file for the FEWLS tool. If using defaults, the tool can be run directly 
from this file with no alteration. The only needed input is *in_solar_df.shp* with the 
necessary attributes. No model functions are called here. Instead model scripts for resource 
(FEWLS_ResourceModelResults.R) and ecnomic (FEWLS_EconModelResults.R) sub-models are called 
and run. Make sure to specify the working directory that contains the FEWLS_tool folder.

Additionally, FEWLS_FarmElecBudgExport.R is called here. This is a non-essential script to 
run every model run, but is needed for the AgrisolarGenNEM.py scripts. It is, therefore, 
currently commented out. cleanFEWLS() is also a non-essential function, but clears model
memory of certain input dataframes (described in FEWLS_model.R)

Notes for future updates: 
1) Split Econ and Resource model results scripts into dataframe and figure generation. Also, 
   Split Econ and Resource model results scripts into commercial- and utility-scale: Currently,
   model only runs if both commericial- and utility-scale arrays are present. 
  
2) Create default model function where defaults can be changed in an automated fashion for 
   output (ie runFEWLS() vs. runFEWLS(system_lifespan = 35, discount_rate = 15%))
"
#_______________________________________________________________________________________________________________________________________ Model Setup

# Set working directory
wd <- 'S://Users/stidjaco/R_files/FEWLS_tool'
setwd(wd)

# Packages for Processing 
library(sf)
library(lwgeom)
library(tidyverse)
library(base)
library(dplyr)
library(stringr)
library(showtext)

# Model Name: Give it a name! This appends the model name to all exported files -- Capacity and System Lifespan are inherently included in model export names already
model_name = "baseline"

# Pull in solar dataset -- stringsAsFactors must be false, must be sf object
# Must contain a shape file with the attributes: geometry, Year, Capacity, dir_a, irrig, Crp_yr_(1-5) 
in_solar_df <- st_read('CV_Agrisolar_df/CV_Agrisolar_df.shp', stringsAsFactors = FALSE, quiet = TRUE) %>% st_as_sf()

# Compile model inputs, model, and model result outputs
source("FEWLS_model.R")

# Run models and get outputs
source("FEWLS_ResourceModelResults.R")
source("FEWLS_EconModelResults.R")

# "Start from Scratch" -- Function to clear FEWLS memory of generated data, model checks if these exist and forgos running them if they do
#cleanFEWLS()

# Get farm-level electricity budget and econ inputs (Agrisolar Manuscript) # -- Not directly used in scripts, should move to supp. Outputs three dataframes and fig described in README
source("FEWLS_FarmElecBudgExport.R") # Does not need to be run every model run
