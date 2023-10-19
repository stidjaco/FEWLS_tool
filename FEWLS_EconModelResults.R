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
This file executes all model functions derived in *FEWLS_model.R* for economic implications
of lifespan agrisolar co-location for both commerical- and utility-scale installations. 
Exports dataframe containing all cumulative econoimc values for every array during every 
year of the lifespan (EconModel_[Comm|Util]ArrayResults_##yrLS_##MW.csv), a dataframe of 
grouped cumulative FEW resource reuslts (FEW_Resource_Impacts_##yrLS_##MW.csv), along with
a .png of the same dataframe, and an economic impact figure.
"
# Get economic model results.

# --------------------------- #
# Get dataset info and figure # 
# --------------------------- #
# Get start time
start_time <- Sys.time()

# Separate by threshold
in_df <- in_solar_df
Commercial_df <- in_df[which(in_df$Capacity<=capacity_threshold), ]
Utility_df <- in_df[which(in_df$Capacity>capacity_threshold), ]

# Get dataset info
getData_info(in_solar_df)
#getPlot_Array_Tech_Dist(in_solar_df) # Only needed on first run for dataset information


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                       DATAFRAMES                                          #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


####### Get input dataset info (Commercial) #######
commercial <- data.frame(Year = integer(), food_min = numeric(), food_max = numeric(), food_base = numeric(), energy_base = numeric(), energy_max = numeric(), energy_min = numeric(), water_base = numeric(), 
                         water_max = numeric(), water_min = numeric(), operation_base = numeric(), operation_max = numeric(), operation_min = numeric(), oandm_base = numeric(), oandm_max = numeric(), 
                         oandm_min = numeric(), install_base = numeric(), install_max = numeric(), install_min = numeric(), Tot_Budg_base = numeric(), Tot_Budg_min = numeric(), Tot_Budg_max = numeric(), 
                         Index = integer(), Crop = character(), Yr_inst = integer(), Capacity = numeric())
# For every commerical installation, run the economic model
for(i in c(1:nrow(Commercial_df))){
  # Get array
  in_df <- Commercial_df[i, ]
  # Source non user inputs
  source("FEWLS_NonUserInputs.R", local = parent.frame())
  # Run commercial econ model
  commercial_temp <- getEconomic_Commercial(in_df)
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
# Convert values using scale factor
commercial[, 2:ncol(commercial)] <- commercial[, 2:ncol(commercial)] / econ_scale



####### Get input dataset info (Utility) #######
# Ensure in_solar_df returns to in_df
in_df <- in_solar_df
# Build new dataframe for utility output
utility <- data.frame(Year = integer(), food_min = numeric(), food_max = numeric(), food_base = numeric(), water_base = numeric(), water_max = numeric(), water_min = numeric(), operation_base = numeric(),
                      operation_max = numeric(), operation_min = numeric(), lease_base = numeric(), lease_max = numeric(), lease_min = numeric(), Tot_Budg_base = numeric(), Tot_Budg_min = numeric(),
                      Tot_Budg_max = numeric(), Index = integer(), Crop = character(), Yr_inst = integer(), Capacity = numeric())
# For every utility scale installation, run the economic model 
for(i in c(1:nrow(Utility_df))){
  # Get array
  in_df <-  Utility_df[i, ]
  # Source non user inputs
  source("FEWLS_NonUserInputs.R", local = parent.frame())
  # Run Utility econ model
  utility_temp <- getEconomic_Utility(in_df)
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
# Convert values using scale factor
utility[, 2:ncol(utility)] <- utility[, 2:ncol(utility)] / econ_scale

# Export Commercial and Utility Dataframes
write.csv(commercial_save, paste("Outputs/Dataframes/EconModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""), row.names=FALSE)
write.csv(utility_save, paste("Outputs/Dataframes/EconModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""), row.names=FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Prepare export dataframes for manuscript 


# Get resource model results into clean dataframe 
# ------------------------ #
# Get economic projections # 
# ------------------------ #
# Run Econ functions
comm_df <- commercial
util_df <- utility

# Create dataframe to export
mil = 1e-3 # already in billion, get to million
rnd = 0 # decimals to round to
df <- data.frame(Resource = character(), Scale = character(), Base = numeric(), Best = numeric(), Worst = numeric())
df <- add_row(df, Resource = "Fg Food Revenue", Scale = "Comm", Base = round(tail(comm_df$food_base, 1)/mil, rnd), Best = round(tail(comm_df$food_min, 1)/mil, rnd), Worst = round(tail(comm_df$food_max, 1)/mil, rnd))
df <- add_row(df, Resource = "Fg Food Revenue", Scale = "Utility", Base = round(tail(util_df$food_base, 1)/mil, rnd), Best = round(tail(util_df$food_min, 1)/mil, rnd), Worst = round(tail(util_df$food_max, 1)/mil, rnd))
df <- add_row(df, Resource = "Electricity Value", Scale = "Comm", Base = round(tail(comm_df$energy_base, 1)/mil, rnd), Best = round(tail(comm_df$energy_max, 1)/mil, rnd), Worst = round(tail(comm_df$energy_min, 1)/mil, rnd))
df <- add_row(df, Resource = "Landlease Value", Scale = "Utility", Base = round(tail(util_df$lease_base, 1)/mil, rnd), Best = round(tail(util_df$lease_max, 1)/mil, rnd), Worst = round(tail(util_df$lease_min, 1)/mil, rnd))
df <- add_row(df, Resource = paste("\u0394", "Water Use Value", sep = ""), Scale = "Comm", Base = round(tail(comm_df$water_base, 1)/mil, rnd+1), Best = round(tail(comm_df$water_max, 1)/mil, rnd+1), Worst = round(tail(comm_df$water_min, 1)/mil, rnd+1))
df <- add_row(df, Resource = paste("\u0394", "Water Use Value", sep = ""), Scale = "Utility", Base = round(tail(util_df$water_base, 1)/mil, rnd+1), Best = round(tail(util_df$water_max, 1)/mil, rnd+1), Worst = round(tail(util_df$water_min, 1)/mil, rnd+1))
df <- add_row(df, Resource = "Fg Operation Cost", Scale = "Comm", Base = round(tail(comm_df$operation_base, 1)/mil, rnd), Best = round(tail(comm_df$operation_max, 1)/mil, rnd), Worst = round(tail(comm_df$operation_min, 1)/mil, rnd))
df <- add_row(df, Resource = "Fg Operation Cost", Scale = "Utility", Base = round(tail(util_df$operation_base, 1)/mil, rnd), Best = round(tail(util_df$operation_max, 1)/mil, rnd), Worst = round(tail(util_df$operation_min, 1)/mil, rnd))
df <- add_row(df, Resource = "O&M Cost", Scale = "Comm", Base = round(tail(comm_df$oandm_base, 1)/mil, rnd), Best = round(tail(comm_df$oandm_min, 1)/mil, rnd), Worst = round(tail(comm_df$oandm_max, 1)/mil, rnd))
df <- add_row(df, Resource = "Installation Cost", Scale = "Comm", Base = round(tail(comm_df$install_base, 1)/mil, rnd), Best = round(tail(comm_df$install_min, 1)/mil, rnd), Worst = round(tail(comm_df$install_max, 1)/mil, rnd))
df <- add_row(df, Resource = "Total Budget", Scale = "Comm", Base = round(tail(comm_df$Tot_Budg_base, 1)/mil, rnd), Best = round(tail(comm_df$Tot_Budg_max, 1)/mil, rnd), Worst = round(tail(comm_df$Tot_Budg_min, 1)/mil, rnd))
df <- add_row(df, Resource = "Total Budget", Scale = "Utility", Base = round(tail(util_df$Tot_Budg_base, 1)/mil, rnd), Best = round(tail(util_df$Tot_Budg_max, 1)/mil, rnd), Worst = round(tail(util_df$Tot_Budg_min, 1)/mil, rnd))
# Make negative budget items negative
df$Base <- ifelse(df$Resource %in% c("Fg Food Revenue", "O&M Cost", "Installation Cost"), df$Base * -1, df$Base)
df$Best <- ifelse(df$Resource %in% c("Fg Food Revenue", "O&M Cost", "Installation Cost"), df$Best * -1, df$Best)
df$Worst <- ifelse(df$Resource %in% c("Fg Food Revenue", "O&M Cost", "Installation Cost"), df$Worst * -1, df$Worst)

# Get per MW and per year values
m2_ha = 10000
df$Per_MW <- ifelse(df$Scale=="Comm", df$Base / sum(Commercial_df$Capacity), df$Base / sum(Utility_df$Capacity)) %>% round(2)
df$Per_ha <- ifelse(df$Scale=="Comm", df$Base / sum(Commercial_df$dir_a/m2_ha), df$Base / sum(Utility_df$dir_a/m2_ha)) %>% round(3)
df$Per_kW_Yr <- ifelse(df$Scale=="Comm", df$Base / sum(Commercial_df$Capacity) / system_lifespan, df$Base / sum(Utility_df$Capacity) / system_lifespan) %>% round(9) # NOTE: Here, results are in Unit/MW/yr, not converted until line 209

# Correct irrigated vs non-irrigated arrays
df$Per_MW <- ifelse(df$Resource %in% c("Fg Irrig Cost") & df$Scale=="Comm", df$Base / sum(Commercial_df$Capacity[which(Commercial_df$irrig==1)]), df$Per_MW) %>% round(2)
df$Per_MW <- ifelse(df$Resource %in% c("Fg Irrig Cost") & df$Scale=="Util", df$Base / sum(Utility_df$Capacity[which(Utility_df$irrig==1)]), df$Per_MW) %>% round(2)
df$Per_ha <- ifelse(df$Resource %in% c("Fg Irrig Cost") & df$Scale=="Comm", df$Base / sum(Commercial_df$dir_a[which(Commercial_df$irrig==1)]/m2_ha), df$Per_ha) %>% round(3)
df$Per_ha <- ifelse(df$Resource %in% c("Fg Irrig Cost") & df$Scale=="Util", df$Base / sum(Utility_df$dir_a[which(Utility_df$irrig==1)]/m2_ha), df$Per_ha) %>% round(3)
df$Per_kW_Yr <- ifelse(df$Resource %in% c("Fg Irrig Cost") & df$Scale=="Comm", df$Base / sum(Commercial_df$Capacity[which(Commercial_df$irrig==1)]) / system_lifespan, df$Per_kW_Yr) %>% round(9)
df$Per_kW_Yr <- ifelse(df$Resource %in% c("Fg Irrig Cost") & df$Scale=="Util", df$Base / sum(Utility_df$Capacity[which(Utility_df$irrig==1)]) / system_lifespan, df$Per_kW_Yr) %>% round(9)
df$Unit <- ifelse(df$Resource %in% c("Fg Irrig Cost"), "*2018 Mil USD", "2018 Mil USD")

# Simplify significant figures
# Set to three significant figures for base/best/worst
df$Base <- df$Base %>% signif(3) %>% format(scientific = F) %>% as.numeric()
df$Best <- df$Best %>% signif(3) %>% format(scientific = F) %>% as.numeric()
df$Worst <- df$Worst %>% signif(3) %>% format(scientific = F) %>% as.numeric()
df$Per_MW <- df$Per_MW %>% signif(3) %>% format(scientific = F) %>% as.numeric()
df$Per_ha <- df$Per_ha %>% signif(3) %>% format(scientific = F) %>% as.numeric()
df$Per_kW_Yr <- (df$Per_kW_Yr / kW_MW * 1e6) %>% round(2) %>% signif(3) %>% format(scientific = F) %>% as.numeric() # funky first portion is to get units into $/kW/yr in final df
names(df)[which(names(df)=="Per_kW_Yr")] <- "$/kW/Yr"
names(df)[which(names(df)=="Per_ha")] <- "/ha"
names(df)[which(names(df)=="Per_MW")] <- "/MW"
write.csv(df, paste("Outputs/Dataframes/Economic_Impacts_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""), row.names=FALSE)

# Clean pub datatable
library(gt)
df_export <- df %>% gt() %>% fmt_currency(columns = 3:5,currency = "USD",decimals = 0) %>% fmt_currency(columns = 6,currency = "USD",decimals = 2) %>% fmt_currency(columns = 7,currency = "USD",decimals = 3) %>% fmt_currency(columns = 8,currency = "USD",decimals = 2)
gtsave(df_export, paste("Outputs/Dataframes/Economic_Impacts_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".png", sep = ""))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                        FIGURES                                            #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COMMERCIAL ECON PLOTTING

# Ensure export quality, call in dataframes from above
commercial <- read.csv(paste("Outputs/Dataframes/EconModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
utility <- read.csv(paste("Outputs/Dataframes/EconModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))

# Group by year & Convert values using scale factor
commercial <- commercial %>% dplyr::select(-matches("Crop|Index|Yr_inst|Capacity")) %>% group_by(Year) %>% summarise_all(sum, na.rm=TRUE)
commercial[, 2:ncol(commercial)] <- commercial[, 2:ncol(commercial)] / econ_scale
utility <- utility %>% dplyr::select(-matches("Crop|Index|Yr_inst|Capacity")) %>% group_by(Year) %>% summarise_all(sum, na.rm=TRUE)
utility[, 2:ncol(utility)] <- utility[, 2:ncol(utility)] / econ_scale

# Return in_df to in_solar_df
in_df <- in_solar_df

# Get max and min y's 
max_y = (dplyr::select(commercial, -matches("Year|install|food|oandm")) %>% max() %>% signif(2) %>% plyr::round_any(0.2, f=ceiling)) %>% as.numeric()
min_y = (dplyr::select(commercial, matches("install|food|oandm")) %>% max() %>% signif(2) %>% plyr::round_any(0.2, f=ceiling) * -1) %>% as.numeric()
break_rate = 0.3

# Rectangles for phases of operational LCA
start_add <-c(range(in_df$Year)) %>% min() %>% as.numeric()
start_constant <- c(range(in_df$Year)) %>% max() %>% as.numeric()
start_removal <- start_add + system_lifespan - 1

# Get economic plot for commercial arrays
commercial_economic_plot <- ggplot(commercial, aes(x=Year)) + 
  annotate(geom = "rect", xmin = start_add, xmax = start_constant, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray75") +
  annotate(geom = "rect", xmin = start_constant, xmax = start_removal, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray50") +
  annotate(geom = "rect", xmin = start_removal, xmax = start_constant+system_lifespan-1, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray25") +
  geom_ribbon(aes(ymin = Tot_Budg_min, ymax = Tot_Budg_max), fill = 'green', alpha = 0.2) +
  geom_ribbon(aes(ymin = energy_min, ymax = energy_max), fill = energy_col, alpha = 0.2) +
  geom_ribbon(aes(ymin = install_max * -1, ymax = install_min * -1), fill = 'red', alpha = 0.2) +
  geom_ribbon(aes(ymin = oandm_max * -1, ymax = oandm_min * -1), fill = 'black', alpha = 0.2) +
  geom_ribbon(aes(ymin = food_max * -1, ymax = food_min * -1), fill = food_col, alpha = 0.2) +
  geom_ribbon(aes(ymin = water_min, ymax = water_max), fill = water_col, alpha = 0.2) +
  geom_ribbon(aes(ymin = operation_min, ymax = operation_max), fill = 'brown', alpha = 0.2) +
  geom_line(aes(y=Tot_Budg_base, color = "Total Budget"), size = 1) + 
  geom_line(aes(y=energy_base, color = "$ Energy Generated"), size = 1) + 
  geom_line(aes(y=water_base, color = "$ Forgone Irrigation"), size = 1) + 
  geom_line(aes(y=oandm_base * -1, color = "$ O&M"), size = 1) + 
  geom_line(aes(y=food_base * -1, color = "$ Forgone Food"), size = 1) + 
  geom_line(aes(y=install_base * -1, color = "$ Installation Cost"), size = 1) +
  geom_line(aes(y=operation_base, color = "$ Forgone Operation"), size = 1) +
  xlab("Year") + 
  ylab("Real Economic Effects (Billion 2018 USD)") +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  #scale_y_continuous(breaks=seq(min_y, max_y, 0.25), limits = c(-0.75, 1.75), expand = c(0, 0)) +
  scale_y_continuous(breaks=seq(signif(min_y, 2), signif(max_y, 2), break_rate), limits = c(min_y*1.025, max_y*1.025), expand = c(0, 0), labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("Total Budget" = "green", "$ Forgone Food" = food_col, "$ Energy Generated" = energy_col, "$ Forgone Irrigation" = water_col, "$ O&M" = "black", "$ Installation Cost" = "red", "$ Forgone Operation" = "brown")) +
  labs(color = "Budget Item") +
  geom_hline(yintercept = 0) +
  theme(plot.title = element_blank(), legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = border_size))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ UTILITY PLOTTING
# Get max and min y's 
max_y = (dplyr::select(utility, -matches("Year|install|food|oandm")) %>% max() %>% signif(2) %>% plyr::round_any(0.2, f=ceiling)) %>% as.numeric()
min_y = (dplyr::select(utility, matches("install|food|oandm")) %>% max() %>% signif(2) %>% plyr::round_any(0.2, f=ceiling) * -1) %>% as.numeric()
break_rate = 0.2

# Get economic plot for utility arrays
utility_economic_plot <- ggplot(utility, aes(x=Year)) + 
  annotate(geom = "rect", xmin = start_add, xmax = start_constant, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray75") +
  annotate(geom = "rect", xmin = start_constant, xmax = start_removal, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray50") +
  annotate(geom = "rect", xmin = start_removal, xmax = start_constant+system_lifespan-1, ymin = -Inf, ymax = Inf, alpha = phase_alpha, fill = "gray25") +
  geom_ribbon(aes(ymin = Tot_Budg_min, ymax = Tot_Budg_max), fill = 'green', alpha = 0.2) +
  geom_ribbon(aes(ymin = lease_min, ymax = lease_max), fill = energy_col, alpha = 0.2) +
  geom_ribbon(aes(ymin = food_max * -1, ymax = food_min * -1), fill = food_col, alpha = 0.2) +
  geom_ribbon(aes(ymin = water_min, ymax = water_max), fill = water_col, alpha = 0.2) +
  geom_ribbon(aes(ymin = operation_min, ymax = operation_max), fill = 'brown', alpha = 0.2) +
  geom_line(aes(y=Tot_Budg_base, color = "Total Budget"), size = 1) + 
  geom_line(aes(y=lease_base, color = "$ Lease Generated"), size = 1) + 
  geom_line(aes(y=water_base, color = "$ Forgone Irrigation"), size = 1) + 
  geom_line(aes(y=food_base * -1, color = "$ Forgone Food"), size = 1) + 
  geom_line(aes(y=operation_base, color = "$ Forgone Operation"), size = 1) +
  xlab("Year") + 
  ylab("Net Dollar Value (2018 USD)") +
  scale_x_continuous(breaks=seq(start_add, start_constant+system_lifespan, 8), expand = c(0, 0)) +
  #scale_y_continuous(breaks=seq(min_y, max_y, 0.25), limits = c(-0.75, 1.75), expand = c(0, 0)) +
  scale_y_continuous(breaks=seq(signif(min_y, 2), signif(max_y, 2), break_rate), limits = c(min_y*1.025, max_y*1.025), expand = c(0, 0), labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("Total Budget" = "green", "$ Forgone Food" = food_col, "$ Lease Generated" = energy_col, "$ Forgone Irrigation" = water_col, "$ Forgone Operation" = "brown")) +
  labs(color = "Budget Item") +
  geom_hline(yintercept = 0) +
  theme(plot.title = element_blank(), legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = border_size))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Generate Save Plot
# Arrange & Save
fig <- ggarrange(commercial_economic_plot, 
                 utility_economic_plot+theme(axis.title.y = element_blank()), nrow = 1, ncol = 2, widths = c(9,10), align = c("hv")) 
fig <- annotate_figure(fig, bottom = text_grob("Year", size = 12, hjust = 0))
ggsave(path = 'Outputs/Figures', filename = paste("ECON_results_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".pdf", sep=""), width = 7, height = 3.25, units = "in", useDingbats=FALSE)

# Get runtime
end_time <- Sys.time()
econ_tot_time = (end_time - start_time) %>% format(format = "%H:%M:%S")
print(paste(model_name, " Econ-Model Runtime: ", econ_tot_time, ". See '", wd, "/Outputs' for results.", sep=""))