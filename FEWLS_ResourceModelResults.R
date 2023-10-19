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
This file executes all model functions derived in *FEWLS_model.R* for resource implications
(food, energy, water) of lifespan agrisolar co-location for both commerical- and utility- 
scale installations. Exports dataframe containing all cumulative resource values for every 
array during every year of the lifespan (ResourceModel_[Comm|Util]ArrayResults_##yrLS_##MW.csv),
a dataframe of grouped cumulative FEW resource reuslts (FEW_Resource_Impacts_##yrLS_##MW.csv), 
along with a .png of the same dataframe, and a FEW impact figure.

Additionally, getPlot_Array_Tech_Dist() is called here. This is a non-essential script to run 
every model run, but is valuable for visualizing the capacity distribution of your dataset and 
is figure 2 in the Agrisolar Manuscript. 
"

# --------------------------- #
# Get dataset info and figure # 
# --------------------------- #
# Get start time
start_time <- Sys.time()

# Separate by threshold
in_df <- in_solar_df
Commercial_df <- in_df[which(in_df$Capacity<=capacity_threshold), ]
Utility_df <- in_df[which(in_df$Capacity>capacity_threshold), ]

# Get dataset info -- Does not need to be run every model
getData_info(in_solar_df)
#getPlot_Array_Tech_Dist(in_solar_df)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                       DATAFRAMES                                          #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


####### Get input dataset info (Commercial) #######
commercial <- data.frame(Year = integer(), kcal = numeric(), kcal_min = numeric(), kcal_max = numeric(), energy_base = numeric(), energy_max = numeric(), energy_min = numeric(), 
                         irrEngy_fg = numeric(), IrrEngy_fg_dry = numeric(), irrEngy_fg_wet = numeric(), irrig_fg = numeric(), irrig_fg_dry = numeric(), irrig_fg_wet = numeric(), 
                         oandm_base = numeric(), oandm_min = numeric(), oandm_max = numeric(), oandmIr_base = numeric(), oandmIr_min = numeric(), oandmIr_max = numeric(), 
                         Index = integer(), Crop = character(), Yr_inst = integer(), Capacity = numeric())
# For every commerical installation, run the economic model
for(i in c(1:nrow(Commercial_df))){
  # Get array
  in_df <- Commercial_df[i, ]
  # Source non user inputs
  source("FEWLS_NonUserInputs.R", local = parent.frame())
  # Run commercial resource model
  commercial_temp <- getResource_CommUtil(in_df)
  # Append to commercial df
  commercial <- rbind(commercial, commercial_temp) 
  print(paste(round(i/nrow(Commercial_df)*100, 2), "% done")) } #  -- printing progress reduces efficiency
# For every econ-variable and every installation, extend vector to include full range of model
for(indx in unique(commercial$Index)){
  # Subset array econ model and get total model range
  array <- commercial[which(commercial$Index==indx), ]
  model_range <- c(range(commercial$Year)[1]:range(commercial$Year)[2])
  # Save decomission year and final econ values
  array_decomYr <- tail(array$Year, 1) %>% as.numeric()
  array_endVals <- tail(array, 1) %>% dplyr::select(-matches("Year|Yr_inst|Crop|Index|Capacity"))
  # If vectors differ, add rows to dataframe with missing years
  if(setequal(array$Year, model_range)==FALSE){
    missing_years <- setdiff(model_range, array$Year)
    for(year in missing_years){
      array <- add_row(array, Year=year, Index=first(array$Index), Crop=first(array$Crop), Yr_inst=first(array$Yr_inst), Capacity = first(array$Capacity))
    }
    # Order by year for cumsum operations
    array <- array[order(array$Year), ]
    array[is.na(array)] <- 0
  }
  # If new year-row is after decommission, retain final operational year economic value for each variable (because this is a cumulative analysis)
  array[which(array$Year>array_decomYr), -which(names(array) %in% c("Year","Yr_inst","Crop","Index", "Capacity"))] <- array_endVals
  # Return expanded array data to main df
  commercial <- commercial[which(!commercial$Index==indx), ]
  commercial <- rbind(commercial, array)
}
# Save df for crop or index export
rownames(commercial) <- c(1:nrow(commercial))
commercial <- commercial %>% as.data.frame()
commercial_save <- commercial
# Group by year
commercial <- commercial %>% dplyr::select(-matches("Crop|Index|Yr_inst|Capacity")) %>% group_by(Year) %>% summarise_all(sum, na.rm=TRUE)



####### Get input dataset info (Utility) #######
# Ensure in_solar_df returns to in_df
in_df <- in_solar_df
# Build new dataframe for utility output
utility <- data.frame(Year = integer(), kcal = numeric(), kcal_min = numeric(), kcal_max = numeric(), energy_base = numeric(), energy_max = numeric(), energy_min = numeric(), 
                      irrEngy_fg = numeric(), IrrEngy_fg_dry = numeric(), irrEngy_fg_wet = numeric(), irrig_fg = numeric(), irrig_fg_dry = numeric(), irrig_fg_wet = numeric(), 
                      oandm_base = numeric(), oandm_min = numeric(), oandm_max = numeric(), oandmIr_base = numeric(), oandmIr_min = numeric(), oandmIr_max = numeric(), 
                      Index = integer(), Crop = character(), Yr_inst = integer(), Capacity = numeric())
# For every utility scale installation, run the economic model 
for(i in c(1:nrow(Utility_df))){
  # Get array
  in_df <-  Utility_df[i, ]
  # Source non user inputs
  source("FEWLS_NonUserInputs.R", local = parent.frame())
  # Run commercial econ model
  utility_temp <- getResource_CommUtil(in_df)
  # Append to utility df
  utility <- rbind(utility, utility_temp) 
  print(paste(round(i/nrow(Utility_df)*100, 2), "% done")) #  -- printing progress reduces efficiency
}
# For every econ-variable and every installation, extend vector to include full range of model
for(indx in unique(utility$Index)){
  # Subset array econ model and get total model range
  array <- utility[which(utility$Index==indx), ]
  model_range <- c(range(utility$Year)[1]:range(utility$Year)[2])
  # Save decomission year and final econ values
  array_decomYr <- tail(array$Year, 1) %>% as.numeric()
  array_endVals <- tail(array, 1) %>% dplyr::select(-matches("Year|Yr_inst|Crop|Index|Capacity"))
  # If vectors differ, add rows to dataframe with missing years
  if(setequal(array$Year, model_range)==FALSE){
    missing_years <- setdiff(model_range, array$Year)
    for(year in missing_years){
      array <- add_row(array, Year=year, Index=first(array$Index), Crop=first(array$Crop), Yr_inst=first(array$Yr_inst), Capacity = first(array$Capacity))
    }
    # Order by year for cumsum operations
    array <- array[order(array$Year), ]
    array[is.na(array)] <- 0
  }
  # If new year-row is after decommission, retain final operational year economic value for each variable (because this is a cumulative analysis)
  array[which(array$Year>array_decomYr), -which(names(array) %in% c("Year","Yr_inst","Crop","Index", "Capacity"))] <- array_endVals
  # Return expanded array data to main df
  utility <- utility[which(!utility$Index==indx), ]
  utility <- rbind(utility, array)
}
# Save df for crop or index export
rownames(utility) <- c(1:nrow(utility))
utility <- utility %>% as.data.frame()
utility_save <- utility
# Group by year
utility <- utility %>% dplyr::select(-matches("Crop|Index|Yr_inst|Capacity")) %>% group_by(Year) %>% summarise_all(sum, na.rm=TRUE)

# Export Commercial and Utility Dataframes -- Raw model output, every year and every array 
write.csv(commercial_save, paste("Outputs/Dataframes/ResourceModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""), row.names=FALSE)
write.csv(utility_save, paste("Outputs/Dataframes/ResourceModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""), row.names=FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Prepare export dataframes for manuscript 


# Get resource model results into clean dataframe 
# ------------------------ #
# Get resource projections # 
# ------------------------ #

# FOOD
commercial_kcal_df <- commercial %>% dplyr::select(matches("kcal"))
utility_kcal_df <- utility %>% dplyr::select(matches("kcal"))

# ENERGY PRODUCED AND CONSERVED FROM IRRIGATION
commercial_energy_df <- commercial %>% dplyr::select(matches("energy"))
utility_energy_df <- utility %>% dplyr::select(matches("energy"))
commercial_irrEngy_df <- commercial %>% dplyr::select(matches("irrEngy"))
utility_irrEngy_df <- utility %>% dplyr::select(matches("irrEngy"))

# WATER USED FOR O&M AND CONSERVED FROM IRRIGATION
commercial_irrig_df <- commercial %>% dplyr::select(matches("irrig_fg"))
utility_irrig_df <- utility %>% dplyr::select(matches("irrig_fg"))
commercial_OandM_df <- commercial %>% dplyr::select(matches("oandm_"))
utility_OandM_df <- utility %>% dplyr::select(matches("oandm_"))
# For only irrigated arrays
commercial_OandM_df_irrig <- commercial %>% dplyr::select(matches("oandmIr_"))
utility_OandM_df_irrig <- utility %>% dplyr::select(matches("oandmIr_"))

# Create dataframe to export
bil = 1e9
rnd = 1 # decimals to round to
df <- data.frame(Resource = character(), Scale = character(), Base = numeric(), Best = numeric(), Worst = numeric(), Unit = character())
df <- add_row(df, Resource = "Shifted Calories", Scale = "Comm", Base = round(tail(commercial_kcal_df$kcal, 1)/bil, rnd), Best = round(tail(commercial_kcal_df$kcal_min, 1)/bil, rnd), Worst = round(tail(commercial_kcal_df$kcal_max, 1)/bil, rnd), Unit = "Bil Kcal")
df <- add_row(df, Resource = "Shifted Calories", Scale = "Utility", Base = round(tail(utility_kcal_df$kcal, 1)/bil, rnd), Best = round(tail(utility_kcal_df$kcal_min, 1)/bil, rnd), Worst = round(tail(utility_kcal_df$kcal_max, 1)/bil, rnd), Unit = "Bil Kcal")
df <- add_row(df, Resource = "Elect. Produced", Scale = "Comm", Base = round(tail(commercial_energy_df$energy_base, 1), rnd), Best = round(tail(commercial_energy_df$energy_max, 1), rnd), Worst = round(tail(commercial_energy_df$energy_min, 1), rnd), Unit = "GWh")
df <- add_row(df, Resource = "Elect. Produced", Scale = "Utility", Base = round(tail(utility_energy_df$energy_base, 1), rnd), Best = round(tail(utility_energy_df$energy_max, 1), rnd), Worst = round(tail(utility_energy_df$energy_min, 1), rnd), Unit = "GWh")
df <- add_row(df, Resource = "Irrig. Elect. Saved", Scale = "Comm", Base = round(tail(commercial_irrEngy_df$irrEngy_fg, 1), rnd), Best = round(tail(commercial_irrEngy_df$irrEngy_fg_dry, 1), rnd), Worst = round(tail(commercial_irrEngy_df$irrEngy_fg_wet, 1), rnd), Unit = "*GWh")
df <- add_row(df, Resource = "Irrig. Elect. Saved", Scale = "Utility", Base = round(tail(utility_irrEngy_df$irrEngy_fg, 1), rnd), Best = round(tail(utility_irrEngy_df$irrEngy_fg_dry, 1), rnd), Worst = round(tail(utility_irrEngy_df$irrEngy_fg_wet, 1), rnd), Unit = "*GWh")
df <- add_row(df, Resource = "Irrig. Water Saved", Scale = "Comm", Base = round(tail(commercial_irrig_df$irrig_fg, 1)/m3_millionm3, rnd), Best = round(tail(commercial_irrig_df$irrig_fg_dry, 1)/m3_millionm3, rnd), Worst = round(tail(commercial_irrig_df$irrig_fg_wet, 1)/m3_millionm3, rnd), Unit = "*Mil m3")
df <- add_row(df, Resource = "Irrig. Water Saved", Scale = "Utility", Base = round(tail(utility_irrig_df$irrig_fg, 1)/m3_millionm3, rnd), Best = round(tail(utility_irrig_df$irrig_fg_dry, 1)/m3_millionm3, rnd), Worst = round(tail(utility_irrig_df$irrig_fg_wet, 1)/m3_millionm3, rnd), Unit = "*Mil m3")
df <- add_row(df, Resource = "O&M Water Used", Scale = "Comm", Base = round(tail(commercial_OandM_df$oandm_base, 1)/m3_millionm3, rnd), Best = round(tail(commercial_OandM_df$oandm_min, 1)/m3_millionm3, rnd), Worst = round(tail(commercial_OandM_df$oandm_max, 1)/m3_millionm3, rnd), Unit = "*Mil m3")
df <- add_row(df, Resource = "O&M Water Used", Scale = "Utility", Base = round(tail(utility_OandM_df$oandm_base, 1)/m3_millionm3, rnd), Best = round(tail(utility_OandM_df$oandm_min, 1)/m3_millionm3, rnd), Worst = round(tail(utility_OandM_df$oandm_max, 1)/m3_millionm3, rnd), Unit = "*Mil m3")

# Get per MW and per year values
m2_ha = 10000
df$Per_MW <- ifelse(df$Scale=="Comm", df$Base / sum(Commercial_df$Capacity), df$Base / sum(Utility_df$Capacity)) %>% round(2)
df$Per_ha <- ifelse(df$Scale=="Comm", df$Base / sum(Commercial_df$dir_a/m2_ha), df$Base / sum(Utility_df$dir_a/m2_ha)) %>% round(2)
df$Per_MW_Yr <- ifelse(df$Scale=="Comm", df$Base / sum(Commercial_df$Capacity) / system_lifespan, df$Base / sum(Utility_df$Capacity) / system_lifespan) %>% round(6)

# Correct irrigated vs non-irrigated arrays
df$Per_MW <- ifelse(df$Resource %in% c("Irrig. Elect. Saved", "Irrig. Water Saved") & df$Scale=="Comm", df$Base / sum(Commercial_df$Capacity[which(Commercial_df$irrig==1)]), df$Per_MW) %>% round(2)
df$Per_MW <- ifelse(df$Resource %in% c("Irrig. Elect. Saved", "Irrig. Water Saved") & df$Scale=="Util", df$Base / sum(Utility_df$Capacity[which(Utility_df$irrig==1)]), df$Per_MW) %>% round(2)
df$Per_ha <- ifelse(df$Resource %in% c("Irrig. Elect. Saved", "Irrig. Water Saved") & df$Scale=="Comm", df$Base / sum(Commercial_df$dir_a[which(Commercial_df$irrig==1)]/m2_ha), df$Per_ha) %>% round(2)
df$Per_ha <- ifelse(df$Resource %in% c("Irrig. Elect. Saved", "Irrig. Water Saved") & df$Scale=="Util", df$Base / sum(Utility_df$dir_a[which(Utility_df$irrig==1)]/m2_ha), df$Per_ha) %>% round(2)
df$Per_MW_Yr <- ifelse(df$Resource %in% c("Irrig. Elect. Saved", "Irrig. Water Saved") & df$Scale=="Comm", df$Base / sum(Commercial_df$Capacity[which(Commercial_df$irrig==1)]) / system_lifespan, df$Per_MW_Yr) %>% round(6)
df$Per_MW_Yr <- ifelse(df$Resource %in% c("Irrig. Elect. Saved", "Irrig. Water Saved") & df$Scale=="Util", df$Base / sum(Utility_df$Capacity[which(Utility_df$irrig==1)]) / system_lifespan, df$Per_MW_Yr) %>% round(6)
# Special correction for O&M water use: All scenarios report whole dataset, for perMW, Yr, and Ha, use only irrigated arrays
df$Per_MW <- ifelse(df$Resource %in% c("O&M Water Used") & df$Scale=="Comm", round(tail(commercial_OandM_df_irrig$oandmIr_base, 1)/m3_millionm3, rnd) / sum(Commercial_df$Capacity[which(Commercial_df$irrig==1)]), df$Per_MW) %>% round(2)
df$Per_MW <- ifelse(df$Resource %in% c("O&M Water Used") & df$Scale=="Util", round(tail(utility_OandM_df_irrig$oandmIr_base, 1)/m3_millionm3, rnd) / sum(Utility_df$Capacity[which(Utility_df$irrig==1)]), df$Per_MW) %>% round(2)
df$Per_ha <- ifelse(df$Resource %in% c("O&M Water Used") & df$Scale=="Comm", round(tail(commercial_OandM_df_irrig$oandmIr_base, 1)/m3_millionm3, rnd) / sum(Commercial_df$dir_a[which(Commercial_df$irrig==1)]/m2_ha), df$Per_ha) %>% round(2)
df$Per_ha <- ifelse(df$Resource %in% c("O&M Water Used") & df$Scale=="Util", round(tail(utility_OandM_df_irrig$oandmIr_base, 1)/m3_millionm3, rnd) / sum(Utility_df$dir_a[which(Utility_df$irrig==1)]/m2_ha), df$Per_ha) %>% round(2)
df$Per_MW_Yr <- ifelse(df$Resource %in% c("O&M Water Used") & df$Scale=="Comm", round(tail(commercial_OandM_df_irrig$oandmIr_base, 1)/m3_millionm3, rnd) / sum(Commercial_df$Capacity[which(Commercial_df$irrig==1)]) / system_lifespan, df$Per_MW_Yr) %>% round(6)
df$Per_MW_Yr <- ifelse(df$Resource %in% c("O&M Water Used") & df$Scale=="Util", round(tail(utility_OandM_df_irrig$oandmIr_base, 1)/m3_millionm3, rnd) / sum(Utility_df$Capacity[which(Utility_df$irrig==1)]) / system_lifespan, df$Per_MW_Yr) %>% round(6)
units <- df$Unit 
df$Unit = NULL
df$Unit = units

# Set to three significant figures for base/best/worst
df$Base <- df$Base %>% signif(3) %>% format(scientific = F)
df$Best <- df$Best %>% signif(3) %>% format(scientific = F)
df$Worst <- df$Worst %>% signif(3) %>% format(scientific = F)
df$Per_MW <- df$Per_MW %>% signif(3) %>% format(scientific = F)
df$Per_ha <- df$Per_ha %>% signif(3) %>% format(scientific = F)
df$Per_MW_Yr <- df$Per_MW_Yr %>% round(4) %>% signif(3) %>% format(scientific = F)
names(df)[which(names(df)=="Per_MW_Yr")] <- "Unit/MW/Yr"
names(df)[which(names(df)=="Per_ha")] <- "/ha"
names(df)[which(names(df)=="Per_MW")] <- "/MW"
write.csv(df, paste("Outputs/Dataframes/FEW_Resource_Impacts_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""), row.names=FALSE)

# Clean pub datatable
library(gt)
df_export <- df %>% gt()
gtsave(df_export, paste("Outputs/Dataframes/FEW_Resource_Impacts_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".png", sep = ""))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                        FIGURES                                            #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Prep work for figure generation

# Ensure export quality, call in dataframes from above
commercial_save <- read.csv(paste("Outputs/Dataframes/ResourceModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
utility_save <- read.csv(paste("Outputs/Dataframes/ResourceModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
# Group by year
commercial <- commercial_save %>% dplyr::select(-matches("Crop|Index|Yr_inst|Capacity")) %>% group_by(Year) %>% summarise_all(sum, na.rm=TRUE)
utility <- utility_save %>% dplyr::select(-matches("Crop|Index|Yr_inst|Capacity")) %>% group_by(Year) %>% summarise_all(sum, na.rm=TRUE)

# Return in_df to in_solar_df
in_df <- in_solar_df

# Rectangles for phases of operational LCA
start_add <-c(range(in_df$Year)) %>% min() %>% as.numeric()
start_constant <- c(range(in_df$Year)) %>% max() %>% as.numeric()
start_removal <- start_add + system_lifespan - 1

# Get distribution by crop
getResource_CropDistFig = function(df, resource){
  # Get kcal df
  crop_df <- df
  nass = Nass_Classifications
  nass$Crop <- tolower(nass$Crop)
  crop_df$FARS_crop <- nass$Irrig_Crop[match(tolower(crop_df$Crop), nass$Crop)] %>% tolower()
  # Group further to simplify
  crop_df$FARS_crop <- ifelse(grepl(crp_grp1,crop_df$FARS_crop), crp_nme1, crop_df$FARS_crop)
  crop_df$FARS_crop <- ifelse(grepl(crp_grp2,crop_df$FARS_crop), crp_nme2, crop_df$FARS_crop)
  crop_df$FARS_crop <- ifelse(grepl(crp_grp3,crop_df$FARS_crop), crp_nme3, crop_df$FARS_crop)
  # Select for resource of interest
  crop_df <- crop_df %>% dplyr::select(matches(paste("Year|FARS_crop|", resource, sep = "")))
  return(crop_df)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COMMERCIAL KCAL PLOTTING
# Get max and min y's 
max_kcalC = (commercial %>% dplyr::select(matches("kcal"))/kcal_adj) %>% max() %>% signif(digits = 2)
max_kcalU = (utility %>% dplyr::select(matches("kcal"))/kcal_adj) %>% max() %>% signif(digits = 2)
max_kcal = max(max_kcalC, max_kcalU)
kcal_breaks <- plyr::round_any(max_kcal/7, 100, ceiling) # 8 is the desired number of y axis labels
accuracy = 0.01

# Kcal project plot -- BOTH SCALES
kcal_proj_plot <- ggplot(data=commercial, aes(x=Year)) + 
  annotate(geom = "rect", xmin = start_add, xmax = start_constant, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray75") +
  annotate(geom = "rect", xmin = start_constant, xmax = start_removal, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray50") +
  annotate(geom = "rect", xmin = start_removal, xmax = start_constant+system_lifespan-1, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray25") +
  geom_ribbon(aes(ymin = kcal_max/kcal_adj*-1, ymax = kcal_min/kcal_adj*-1), fill = food_col, alpha = 0.5) +
  geom_line(aes(y = kcal/kcal_adj*-1), color = food_col, size = 1) +
  geom_ribbon(data=utility, aes(ymin = kcal_max/kcal_adj*-1, ymax = kcal_min/kcal_adj*-1), fill = food_col, alpha = 0.5) +
  geom_line(data=utility, aes(y = kcal/kcal_adj*-1), color = food_col, size = 1,  linetype="dashed") +
  ylab(paste("kcal (", "billion",")", sep="")) +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  scale_y_continuous(breaks=seq(max_kcal*-1, 0, kcal_breaks), limits = c(max_kcal*-1*1.025,0), expand = c(0, 0)) +
  theme(plot.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.text.x = element_blank(), 
        plot.margin=margin_plot, panel.border = element_rect(color = "black", fill = NA, size = border_size))

# Get kcal crop distribution -- Commercial
crop_kcal_df <- getResource_CropDistFig(commercial_save, 'kcal')
crop_kcal_df <- crop_kcal_df %>% group_by(Year, FARS_crop) %>% summarise(kcal = sum(kcal))
graph_name = "Commercial"
commercial_kcal_plot <- ggplot(crop_kcal_df, aes(x=Year, y=kcal, fill = FARS_crop)) +
  geom_area(position = "fill", colour = "black", size = .2, alpha = 1, linetype=ifelse(graph_name=="Commercial", "solid", "dashed")) +
  ylab("Fg kcal (% of total)") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  scale_fill_manual(values = crop_plot_cols, name = "kcal (%)") + 
  theme_bw() + theme(axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust=0.5, face="bold")) +
  theme(axis.title.y = element_blank(), axis.text.x = element_blank(), plot.margin=margin_plot, panel.border = element_rect(color = "black", fill = NA, size = border_size))

# Get kcal crop distribution -- Utility
crop_kcal_df <- getResource_CropDistFig(utility_save, 'kcal')
crop_kcal_df <- crop_kcal_df %>% group_by(Year, FARS_crop) %>% summarise(kcal = sum(kcal))
graph_name = "Utility"
utility_kcal_plot <- ggplot(crop_kcal_df, aes(x=Year, y=kcal, fill = FARS_crop)) +
  geom_area(position = "fill", colour = "black", size = .2, alpha = 1, linetype=ifelse(graph_name=="Commercial", "solid", "dashed")) +
  ylab("Fg kcal (% of total)") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  scale_fill_manual(values = crop_plot_cols, name = "kcal (%)") + 
  theme_bw() + theme(axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_blank(), plot.margin=right_margin_plot, panel.border = element_rect(color = "black", fill = NA, size = border_size))

# Get legend to save and remove from plots
legend_save <- get_legend(commercial_kcal_plot) %>% as_ggplot() + theme(plot.margin=legend_margin_plot)
commercial_kcal_plot <- commercial_kcal_plot + theme(legend.position = "none")

# Save plot out to list
Food_fig <- list(kcal_proj_plot, commercial_kcal_plot, utility_kcal_plot, legend_save)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COMMERCIAL ENERGY PLOTTING
# Get max and min y's 
energy_adj = 1000 # MWh to GWh
max_EnergyC = commercial %>% dplyr::select(matches("energy")) %>% max() 
max_EnergyU = utility %>% dplyr::select(matches("energy")) %>% max()
max_Energy = ceiling(max(max_EnergyC, max_EnergyU) / energy_adj) %>% signif(2) %>% plyr::round_any(50, ceiling)
max_IrEnergyC = commercial %>% dplyr::select(matches("irrEngy")) %>% max() 
max_IrEnergyU = utility %>% dplyr::select(matches("irrEngy")) %>% max()
max_IrEnergy = (max(max_IrEnergyC, max_IrEnergyU) / energy_adj) %>% signif(2)
# Set second axis coefficient
ecoeff <- plyr::round_any((max_Energy/max_IrEnergy), 100, f=ceiling)/2

# Kcal project plot -- BOTH SCALES
energy_proj_plot <- ggplot(data=commercial, aes(x=Year)) + 
  annotate(geom = "rect", xmin = start_add, xmax = start_constant, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray75") +
  annotate(geom = "rect", xmin = start_constant, xmax = start_removal, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray50") +
  annotate(geom = "rect", xmin = start_removal, xmax = start_constant+system_lifespan-1, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray25") +
  # Electricity produced
  geom_ribbon(data=commercial, aes(ymin = energy_min/energy_adj, ymax = energy_max/energy_adj), fill = energy_col, alpha = 0.5) + 
  geom_line(data=commercial, aes(y = energy_base/energy_adj), color = energy_col, size = 1) +
  geom_ribbon(data=utility, aes(ymin = energy_min/energy_adj, ymax = energy_max/energy_adj), fill = energy_col, alpha = 0.5) +
  geom_line(data=utility, aes(y = energy_base/energy_adj), color = energy_col, size = 1, linetype="dashed") +
  # Irrigation energy forgone
  geom_ribbon(data = commercial, aes(ymin = irrEngy_fg_wet*ecoeff/energy_adj, ymax = irrEngy_fg_dry*ecoeff/energy_adj), fill = water_col, alpha = 0.5) +
  geom_line(data = commercial, aes(y=irrEngy_fg*ecoeff/energy_adj), size = 1, color = water_col) +
  geom_ribbon(data = utility, aes(ymin = irrEngy_fg_wet*ecoeff/energy_adj, ymax = irrEngy_fg_dry*ecoeff/energy_adj), fill = water_col, alpha = 0.5) +
  geom_line(data = utility, aes(y=irrEngy_fg*ecoeff/energy_adj), size = 1, color = water_col, linetype="dashed") +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  # Graph design
  scale_y_continuous(name = "Electricity (TWh)", sec.axis = sec_axis(~./ecoeff, name="Elec Saved (TWh)"), 
                     breaks=seq(0, max_Energy, 25), limits = c(0, max_Energy*1.025,0), expand = c(0, 0)) +
  theme(plot.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y.right = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.text.y.left = element_text(colour = "orange"), axis.text.y.right = element_text(colour = "blue", angle = 45), plot.margin = margin_plot, panel.border = element_rect(color = "black", fill = NA, size = border_size))

# Get kcal crop distribution -- Commercial
crop_energy_df <- getResource_CropDistFig(commercial_save, 'energy_base')
crop_energy_df <- crop_energy_df %>% group_by(Year, FARS_crop) %>% summarise(energy = sum(energy_base))
graph_name = "Commercial"
commercial_energy_plot <- ggplot(crop_energy_df, aes(x=Year, y=energy, fill = FARS_crop)) +
  geom_area(position = "fill", colour = "black", size = .2, alpha = 1, linetype=ifelse(graph_name=="Commercial", "solid", "dashed")) +
  ylab("Fg kcal (% of total)") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  scale_fill_manual(values = crop_plot_cols, name = "kcal (%)") + 
  theme_bw() + theme(axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust=0.5, face="bold")) +
  theme(axis.title.y = element_blank(), axis.text.x = element_blank(), plot.margin=margin_plot, panel.border = element_rect(color = "black", fill = NA, size = border_size))

# Get kcal crop distribution -- Utility
crop_energy_df <- getResource_CropDistFig(utility_save, 'energy_base')
crop_energy_df <- crop_energy_df %>% group_by(Year, FARS_crop) %>% summarise(energy = sum(energy_base))
graph_name = "Utility"
utility_energy_plot <- ggplot(crop_energy_df, aes(x=Year, y=energy, fill = FARS_crop)) +
  geom_area(position = "fill", colour = "black", size = .2, alpha = 1, linetype=ifelse(graph_name=="Commercial", "solid", "dashed")) +
  ylab("Fg kcal (% of total)") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  scale_fill_manual(values = crop_plot_cols, name = "kcal (%)") + 
  theme_bw() + theme(axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_blank(), plot.margin=right_margin_plot, panel.border = element_rect(color = "black", fill = NA, size = border_size))

# Get legend to save and remove from plots
legend_save <- get_legend(commercial_energy_plot) %>% as_ggplot() + theme(plot.margin=legend_margin_plot)
commercial_energy_plot <- commercial_energy_plot + theme(legend.position = "none")

# Save plot out to list
Energy_fig <- list(energy_proj_plot, commercial_energy_plot, utility_energy_plot, legend_save)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COMMERCIAL WATER PLOTTING
# Get max and min y's 
max_wuC <- commercial %>% dplyr::select(matches("oandm_")) %>% max() 
max_wuU <- utility %>% dplyr::select(matches("oandm_")) %>% max() 
max_wu <- (max(max_wuC, max_wuU) / m3_millionm3) %>% signif(2) %>% signif(2) %>% plyr::round_any(200, f=ceiling)
max_wu_offsetC <- commercial %>% dplyr::select(matches("irrig")) %>% max() 
max_wu_offsetU <- utility %>% dplyr::select(matches("irrig")) %>% max() 
max_wu_offset <- (max(max_wu_offsetC, max_wu_offsetU) / m3_millionm3) %>% signif(2) %>% plyr::round_any(200, f=ceiling)

# Kcal project plot -- BOTH SCALES
water_proj_plot <- ggplot(data=commercial, aes(x=Year)) + 
  annotate(geom = "rect", xmin = start_add, xmax = start_constant, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray75") +
  annotate(geom = "rect", xmin = start_constant, xmax = start_removal, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray50") +
  annotate(geom = "rect", xmin = start_removal, xmax = start_constant+system_lifespan-1, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray25") +
  # Electricity produced
  geom_ribbon(data=commercial, aes(ymin = irrig_fg_dry/m3_millionm3*-1, ymax = irrig_fg_wet/m3_millionm3*-1), fill = water_col, alpha = 0.5) +
  geom_line(data=commercial, aes(y = irrig_fg/m3_millionm3*-1), color = water_col, size = 1) +
  geom_ribbon(data=utility, aes(ymin = irrig_fg_dry/m3_millionm3*-1, ymax = irrig_fg_wet/m3_millionm3*-1), fill = water_col, alpha = 0.5) +
  geom_line(data=utility, aes(y = irrig_fg/m3_millionm3*-1), color = water_col, size = 1,  linetype="dashed") +
  # Irrigation energy forgone
  geom_ribbon(data = commercial, aes(ymin = oandm_min/m3_millionm3, ymax = oandm_max/m3_millionm3), fill = "black", alpha = 0.5) +
  geom_line(data = commercial, aes(y = oandm_base/m3_millionm3), color = "black", size = 1) +
  geom_ribbon(data = utility, aes(ymin = oandm_min/m3_millionm3, ymax = oandm_max/m3_millionm3), fill = "black", alpha = 0.5) +
  geom_line(data = utility, aes(y = oandm_base/m3_millionm3), color = "black", size = 1,  linetype="dashed") +
  # Graph design
  ylab("Water Use (million m3)") +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  scale_y_continuous(breaks=seq(max_wu_offset*-1, max_wu, 200), limits = c(max_wu_offset*-1*1.025,max_wu*1.025), expand = c(0, 0)) +
  theme(plot.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), 
        axis.text.x = element_text(angle=45, hjust=1), plot.margin=margin_plot, panel.border = element_rect(color = "black", fill = NA, size = border_size))

# Get kcal crop distribution -- Commercial
crop_water_df <- getResource_CropDistFig(commercial_save, 'irrig_fg')
crop_water_df <- crop_water_df %>% group_by(Year, FARS_crop) %>% summarise(irrig_fg = sum(irrig_fg))
graph_name = "Commercial"
commercial_water_plot <- ggplot(crop_water_df, aes(x=Year, y=irrig_fg, fill = FARS_crop)) +
  geom_area(position = "fill", colour = "black", size = .2, alpha = 1, linetype=ifelse(graph_name=="Commercial", "solid", "dashed")) +
  ylab("Fg kcal (% of total)") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  scale_fill_manual(values = crop_plot_cols, name = "kcal (%)") + 
  theme_bw() + theme(axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust=0.5, face="bold")) +
  theme(axis.title.y = element_blank(), axis.text.x = element_blank(), plot.margin=margin_plot, panel.border = element_rect(color = "black", fill = NA, size = border_size))

# Get kcal crop distribution -- Utility
crop_water_df <- getResource_CropDistFig(utility_save, 'irrig_fg')
crop_water_df <- crop_water_df %>% group_by(Year, FARS_crop) %>% summarise(irrig_fg = sum(irrig_fg))
graph_name = "Utility"
utility_water_plot <- ggplot(crop_water_df, aes(x=Year, y=irrig_fg, fill = FARS_crop)) +
  geom_area(position = "fill", colour = "black", size = .2, alpha = 1, linetype=ifelse(graph_name=="Commercial", "solid", "dashed")) +
  ylab("Fg kcal (% of total)") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  scale_fill_manual(values = crop_plot_cols, name = "kcal (%)") + 
  theme_bw() + theme(axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_blank(), plot.margin=right_margin_plot, panel.border = element_rect(color = "black", fill = NA, size = border_size))

# Get legend to save and remove from plots
legend_save <- get_legend(commercial_water_plot) %>% as_ggplot() + theme(plot.margin=legend_margin_plot)
commercial_water_plot <- commercial_water_plot + theme(legend.position = "none")

# Save plot out to list
Water_fig <- list(water_proj_plot, commercial_water_plot, utility_water_plot, legend_save)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pull figures into a single exportable fig
# Split figs for modifying widths, heights, margins
kcal_proj_plot <- Food_fig[[1]]
commercial_kcal_plot <- Food_fig[[2]]
utility_kcal_plot <- Food_fig[[3]]
energy_proj_plot <- Energy_fig[[1]]
commercial_energy_plot <- Energy_fig[[2]]
utility_energy_plot <- Energy_fig[[3]]
wu_proj_plot <- Water_fig[[1]]
commercial_irrig_plot <- Water_fig[[2]]
utility_irrig_plot <- Water_fig[[3]]
# Get legend
legend_save <- Food_fig[[4]]


# -------------------------- #
#       Compile figures      # 
# -------------------------- #
fig <- ggarrange(kcal_proj_plot, commercial_kcal_plot, utility_kcal_plot, NULL,
                 energy_proj_plot, commercial_energy_plot, utility_energy_plot, legend_save,
                 wu_proj_plot, commercial_irrig_plot, utility_irrig_plot, NULL, nrow = 3, ncol = 4, widths = c(2,2,2,1), align = c("hv")) 
fig <- annotate_figure(fig, bottom = text_grob("Year", size = 12, hjust = 0))
ggsave(path = 'Outputs/Figures', filename = paste("FEW_results_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".pdf", sep=""), width = 7, height = 5, units = "in", useDingbats=FALSE)
# Get runtime
end_time <- Sys.time()
few_tot_time = (end_time - start_time) %>% format(format = "%H:%M:%S")
print(paste(model_name, " Resource-Model Runtime: ", few_tot_time, ". See '", wd, "/Outputs' for results.", sep=""))