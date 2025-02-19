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
This is a bonus script purely for generating a bar chart-style figure for the
Agrisolar manuscript of the primary FEWLS results. 
"



# Call in result FEW dataframes
commFEW <- read.csv(paste("Outputs/Dataframes/ResourceModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
utilFEW <- read.csv(paste("Outputs/Dataframes/ResourceModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))

# Call in result econ dataframes
commEcon <- read.csv(paste("Outputs/Dataframes/EconModel_CommArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))
utilEcon <- read.csv(paste("Outputs/Dataframes/EconModel_UtilArrayResults_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".csv", sep = ""))

# Filter each dataframe for Year == 2042
utilFEW_2042 <- utilFEW[utilFEW$Year == 2042, ]
commFEW_2042 <- commFEW[commFEW$Year == 2042, ]
utilEcon_2042 <- utilEcon[utilEcon$Year == 2042, ]
commEcon_2042 <- commEcon[commEcon$Year == 2042, ]

# Set variables
geomBarWidth = 1
geomBarBorderWidth = 0.50
errorBarWidth = 0.25
errorBarThickness = 1
hlineLineWidth = 0.5
hlineLineType = "dashed"
hlineColor = 'darkgrey'

# Updated Color Scheme
utilColor <- "#377eb8"  # Muted Blue for Utility
commColor <- "#e41a1c"  # Muted Red for Commercial

# Plot Variables
hjust_label = 1.75
vjust_label = -0.75
border_size = 0.75
border_size = 1

# Adjustment factors
econAdj = 1e3 # thousand of dollars
kcalAdj = 1e6 # Million of kcal
waterAdj = 1e3 # thousand m3
energyAdj = 1e0 # GWh
irrEngyAdj = 1e-3 # MWh
oandmWaterAdj = waterAdj # thousand m3

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Vertical FEW Barchart

# Create labels for each resource section
FEWLabs <- c("Food Production", "Elect. Production", "Irrig. Elect. Use", "Irrig. Water Use", "O&M Water Use")

# Create data for diverging bar chart with weighted avg, min, and max values for Utility and Commercial scale
utilFEWDiverging <- data.frame(
  Metric = FEWLabs,
  Avg = c(-sum(utilFEW_2042$kcal) / sum(utilFEW_2042$Capacity) / kcalAdj / system_lifespan, 
          sum(utilFEW_2042$energy_base) / sum(utilFEW_2042$Capacity) / energyAdj / system_lifespan, 
          -sum(utilFEW_2042$irrEngy_fg) / sum(utilFEW_2042$Capacity) / irrEngyAdj / system_lifespan,
          -sum(utilFEW_2042$irrig_fg) / sum(utilFEW_2042$Capacity) / waterAdj / system_lifespan,
          sum(utilFEW_2042$oandm_base) / sum(utilFEW_2042$Capacity) / oandmWaterAdj / system_lifespan), 
  Min = c(-sum(utilFEW_2042$kcal_min) / sum(utilFEW_2042$Capacity) / kcalAdj / system_lifespan, 
          sum(utilFEW_2042$energy_min) / sum(utilFEW_2042$Capacity) / energyAdj / system_lifespan, 
          -sum(utilFEW_2042$irrEngy_fg_dry) / sum(utilFEW_2042$Capacity) / irrEngyAdj / system_lifespan,
          -sum(utilFEW_2042$irrig_fg_dry) / sum(utilFEW_2042$Capacity) / waterAdj / system_lifespan,
          sum(utilFEW_2042$oandm_min) / sum(utilFEW_2042$Capacity) / oandmWaterAdj / system_lifespan),
  Max = c(-sum(utilFEW_2042$kcal_max) / sum(utilFEW_2042$Capacity) / kcalAdj / system_lifespan, 
          sum(utilFEW_2042$energy_max) / sum(utilFEW_2042$Capacity) / energyAdj / system_lifespan, 
          -sum(utilFEW_2042$irrEngy_fg_wet) / sum(utilFEW_2042$Capacity) / irrEngyAdj / system_lifespan,
          -sum(utilFEW_2042$irrig_fg_wet) / sum(utilFEW_2042$Capacity) / waterAdj / system_lifespan,
          sum(utilFEW_2042$oandm_max) / sum(utilFEW_2042$Capacity) / oandmWaterAdj / system_lifespan))

commFEWDiverging <- data.frame(
  Metric = FEWLabs,
  Avg = c(-sum(commFEW_2042$kcal) / sum(commFEW_2042$Capacity) / kcalAdj / system_lifespan, 
          sum(commFEW_2042$energy_base) / sum(commFEW_2042$Capacity) / energyAdj / system_lifespan, 
          -sum(commFEW_2042$irrEngy_fg) / sum(commFEW_2042$Capacity) / irrEngyAdj / system_lifespan,
          -sum(commFEW_2042$irrig_fg) / sum(commFEW_2042$Capacity) / waterAdj / system_lifespan, 
          sum(commFEW_2042$oandm_base) / sum(commFEW_2042$Capacity) / oandmWaterAdj / system_lifespan),
  Min = c(-sum(commFEW_2042$kcal_min) / sum(commFEW_2042$Capacity) / kcalAdj / system_lifespan, 
          sum(commFEW_2042$energy_min) / sum(commFEW_2042$Capacity) / energyAdj / system_lifespan, 
          -sum(commFEW_2042$irrEngy_fg_dry) / sum(commFEW_2042$Capacity) / irrEngyAdj / system_lifespan,
          -sum(commFEW_2042$irrig_fg_dry) / sum(commFEW_2042$Capacity) / waterAdj / system_lifespan,
          sum(commFEW_2042$oandm_min) / sum(commFEW_2042$Capacity) / oandmWaterAdj / system_lifespan),
  Max = c(-sum(commFEW_2042$kcal_max) / sum(commFEW_2042$Capacity) / kcalAdj / system_lifespan, 
          sum(commFEW_2042$energy_max) / sum(commFEW_2042$Capacity) / energyAdj / system_lifespan, 
          -sum(commFEW_2042$irrEngy_fg_wet) / sum(commFEW_2042$Capacity) / irrEngyAdj / system_lifespan,
          -sum(commFEW_2042$irrig_fg_wet) / sum(commFEW_2042$Capacity) / waterAdj / system_lifespan,
          sum(commFEW_2042$oandm_max) / sum(commFEW_2042$Capacity) / oandmWaterAdj / system_lifespan))

# Set scale
utilFEWDiverging$Scale <- "Utility"
commFEWDiverging$Scale <- "Commercial"

# Combine dataframes
data_FEW_diverging <- rbind(utilFEWDiverging, commFEWDiverging)

# Set the factor levels for Metric to ensure the correct order (reverse actual to account for flip)
data_FEW_diverging$Metric <- factor(data_FEW_diverging$Metric, levels = rev(FEWLabs))

# Create dodged position for bars and error bars
dodge <- position_dodge(width = geomBarWidth)

# Create the diverging bar chart with error bars
FEWPlot <- ggplot(data_FEW_diverging, aes(x = Metric, y = Avg, fill = Scale)) +
  geom_hline(yintercept = 0, color = hlineColor, size = hlineLineWidth, lty = hlineLineType) +  # Zero line
  geom_bar(stat = "identity", position = dodge, color = "black", size = geomBarBorderWidth, width = geomBarWidth) +
  geom_errorbar(aes(ymin = Min, ymax = Max), position = dodge, width = errorBarWidth, size = errorBarThickness) +
  scale_fill_manual(values = c("Utility" = utilColor, "Commercial" = commColor)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Resource") + 
  ylab("FEW Impact (MW-1 yr-1)") +
  theme_pubr() + 
  theme(legend.position = c(0.85, 0.10), axis.title = element_text(face = "bold"), axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = border_size), panel.background = element_blank(), 
        panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
        panel.grid.minor.x = element_blank(),  # Remove minor x-axis grid lines
        panel.grid.major.y = element_line(color = "grey90"),  # Keep y-axis grid lines
        panel.grid.minor.y = element_blank()) + # Remove minor y-axis grid lines (optional)) +
  #coord_flip(ylim = c(-40, 40)) 
  coord_polar()
FEWPlot

# Export
#ggsave(paste("Outputs/Figures/FEWLS_AgrisolarFEWimpacts_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".pdf", sep = ""), FEWPlot + theme(axis.text.y = element_blank()), width = 7, height = 3, dpi = 300)
#ggsave(paste("Outputs/Figures/FEWLS_AgrisolarFEWimpactsAxis_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".pdf", sep = ""), FEWPlot, width = 7, height = 3, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Vertical Econ Barchart

# Create labels for each economic resource section
econLabs <- c("Total Budget", "Food Revenue", "NEM Offset+Return", "Water Use", "O&M Costs", "Operation Costs", "Lease Payments", "Install Costs")

# Create data for diverging bar chart with weighted avg, min, and max values for Utility, with zeros for non-relevant variables
utilEconDiverging <- data.frame(
  Metric = econLabs,
  Avg = c(sum(utilEcon$Tot_Budg_base) / sum(utilEcon$Capacity) / econAdj / system_lifespan,
          -sum(utilEcon$food_base) / sum(utilEcon$Capacity) / econAdj / system_lifespan,  # Food is a loss for utility
          0,  # Energy is not relevant for utility
          sum(utilEcon$water_base) / sum(utilEcon$Capacity) / econAdj / system_lifespan,  # Water is a loss for utility
          0,  # O&M Costs are not relevant for utility
          sum(utilEcon$operation_base) / sum(utilEcon$Capacity) / econAdj / system_lifespan,
          sum(utilEcon$lease_base) / sum(utilEcon$Capacity) / econAdj / system_lifespan, 
          0),  # Install Costs are not relevant for utility
  Min = c(sum(utilEcon$Tot_Budg_min) / sum(utilEcon$Capacity) / econAdj / system_lifespan,
          -sum(utilEcon$food_max) / sum(utilEcon$Capacity) / econAdj / system_lifespan,  # Food min is a loss
          0,  # Energy is not relevant for utility
          sum(utilEcon$water_min) / sum(utilEcon$Capacity) / econAdj / system_lifespan,  # Water min is a loss
          0,  # O&M Costs are not relevant for utility
          sum(utilEcon$operation_min) / sum(utilEcon$Capacity) / econAdj / system_lifespan,
          sum(utilEcon$lease_min) / sum(utilEcon$Capacity) / econAdj / system_lifespan, 
          0),  # Install Costs are not relevant for utility
  Max = c(sum(utilEcon$Tot_Budg_max) / sum(utilEcon$Capacity) / econAdj / system_lifespan,
          -sum(utilEcon$food_min) / sum(utilEcon$Capacity) / econAdj / system_lifespan,  # Food max is a loss
          0,  # Energy is not relevant for utility
          sum(utilEcon$water_max) / sum(utilEcon$Capacity) / econAdj / system_lifespan,  # Water max is a loss
          0,  # O&M Costs are not relevant for utility
          sum(utilEcon$operation_max) / sum(utilEcon$Capacity) / econAdj / system_lifespan,
          sum(utilEcon$lease_max) / sum(utilEcon$Capacity) / econAdj / system_lifespan, 
          0))  # Install Costs are not relevant for utility

# Create data for diverging bar chart with weighted avg, min, and max values for Commercial, with zeros for non-relevant variables
commEconDiverging <- data.frame(
  Metric = econLabs,
  Avg = c(sum(commEcon$Tot_Budg_base) / sum(commEcon$Capacity) / econAdj / system_lifespan,
          -sum(commEcon$food_base) / sum(commEcon$Capacity) / econAdj / system_lifespan,  # Food is a loss for commercial
          sum(commEcon$energy_base) / sum(commEcon$Capacity) / econAdj / system_lifespan, 
          sum(commEcon$water_base) / sum(commEcon$Capacity) / econAdj / system_lifespan,  # Water is a loss for commercial
          -sum(commEcon$oandm_base) / sum(commEcon$Capacity) / econAdj / system_lifespan,  # O&M is a loss for commercial
          sum(commEcon$operation_base) / sum(commEcon$Capacity) / econAdj / system_lifespan,
          0,  # Lease payments are not relevant for commercial
          -sum(commEcon$install_base) / sum(commEcon$Capacity) / econAdj / system_lifespan),  # Install is a loss
  Min = c(sum(commEcon$Tot_Budg_min) / sum(commEcon$Capacity) / econAdj / system_lifespan,
          -sum(commEcon$food_max) / sum(commEcon$Capacity) / econAdj / system_lifespan,  # Food min is a loss
          sum(commEcon$energy_min) / sum(commEcon$Capacity) / econAdj / system_lifespan, 
          sum(commEcon$water_min) / sum(commEcon$Capacity) / econAdj / system_lifespan,  # Water min is a loss
          -sum(commEcon$oandm_max) / sum(commEcon$Capacity) / econAdj / system_lifespan,  # O&M min is a loss
          sum(commEcon$operation_min) / sum(commEcon$Capacity) / econAdj / system_lifespan,
          0,  # Lease payments are not relevant for commercial
          -sum(commEcon$install_min) / sum(commEcon$Capacity) / econAdj / system_lifespan),  # Install min is a loss
  Max = c(sum(commEcon$Tot_Budg_max) / sum(commEcon$Capacity) / econAdj / system_lifespan,
          -sum(commEcon$food_min) / sum(commEcon$Capacity) / econAdj / system_lifespan,  # Food max is a loss
          sum(commEcon$energy_max) / sum(commEcon$Capacity) / econAdj / system_lifespan, 
          sum(commEcon$water_max) / sum(commEcon$Capacity) / econAdj / system_lifespan,  # Water max is a loss
          -sum(commEcon$oandm_min) / sum(commEcon$Capacity) / econAdj / system_lifespan,  # O&M max is a loss
          sum(commEcon$operation_max) / sum(commEcon$Capacity) / econAdj / system_lifespan,
          0,  # Lease payments are not relevant for commercial
          -sum(commEcon$install_max) / sum(commEcon$Capacity) / econAdj / system_lifespan))  # Install max is a loss


# Set scale
utilEconDiverging$Scale <- "Utility"
commEconDiverging$Scale <- "Commercial"

# Combine dataframes
data_econ_diverging <- rbind(utilEconDiverging, commEconDiverging)

# Set the factor levels for Metric to ensure the correct order (reverse actual to account for flip)
data_econ_diverging$Metric <- factor(data_econ_diverging$Metric, levels = rev(econLabs))

# Create dodged position for bars and error bars
dodge <- position_dodge(width = geomBarWidth)

# Create the diverging bar chart with error bars for economic impact
econPlot <- ggplot(data_econ_diverging, aes(x = Metric, y = Avg, fill = Scale)) +
  geom_hline(yintercept = 0, color = hlineColor, size = hlineLineWidth, lty = hlineLineType) +  # Zero line
  geom_bar(stat = "identity", position = dodge, color = "black", size = geomBarBorderWidth, width = geomBarWidth) +
  geom_errorbar(aes(ymin = Min, ymax = Max), position = dodge, width = errorBarWidth, size = errorBarThickness) + 
  scale_fill_manual(values = c("Utility" = utilColor, "Commercial" = commColor)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Resource") + 
  ylab("Economic Impact (Thousand USD MW-1 yr-1)") +
  theme_pubr() + 
  theme(legend.position = c(0.80, 0.20), axis.title = element_text(face = "bold"), axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = border_size), panel.background = element_blank(), 
        panel.grid.major.x = element_blank(),  # Remove major x-axis grid lines
        panel.grid.minor.x = element_blank(),  # Remove minor x-axis grid lines
        panel.grid.major.y = element_line(color = "grey90"),  # Keep y-axis grid lines
        panel.grid.minor.y = element_blank()) + # Remove minor y-axis grid lines (optional)) +
  coord_flip(ylim = c(-160, 160))
econPlot

# Export 
#ggsave(paste("Outputs/Figures/FEWLS_AgrisolarEconimpacts_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".pdf", sep = ""), econPlot + theme(axis.text.y = element_blank()), width = 7, height = 4.8, dpi = 300)
#ggsave(paste("Outputs/Figures/FEWLS_AgrisolarEconimpactsAxis_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".pdf", sep = ""), econPlot, width = 7, height = 4.8, dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Plot Together

library(patchwork)

# Save plots 
plot1 <- FEWPlot + theme(legend.position = 'none')
plot2 <- econPlot 

# Customizing layout proportions
layout <- (plot1 / plot2)  +  plot_layout(heights = c(5, 8), widths = c(1,1))# + plot_annotation(tag_levels = 'A')

# Print the adjusted layout
layout
ggsave(paste("Outputs/Figures/FEWLS_AgrisolarImpactsFull_", system_lifespan, "yrLS_", capacity_threshold, "MW_", model_name, ".pdf", sep = ""), layout, width = 7, height = 6, dpi = 300)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Radar Plot

# Split by crop types
grp_1_crops <- c("grapes", "orchards", "vegetables")
grp_2_crops <- c("cotton", "crops, other", "hay/pasture", "grain")

getRadarCropPlot = function(cropGroup){
  
  # Summarize the contributions of each crop type to food, energy, and water
  radar_data <- commFEW_2042 %>%
    group_by(FARS_crop) %>%
    summarise(
      Food = sum(food_base) / sum(Area) / system_lifespan,
      Energy = sum(energy_base) / sum(Area) / system_lifespan,
      Water = sum(water_base) / sum(Area) / system_lifespan
    )
  
  # Normalize the values to a scale of 0 to 1 for the radar plot
  radar_data_normalized <- radar_data %>%
    mutate(across(Food:Water, ~ (. - min(.)) / (max(.) - min(.))))
  
  # Set row names as the crop type for plotting
  radar_data_normalized <- as.data.frame(radar_data_normalized)
  rownames(radar_data_normalized) <- radar_data$FARS_crop
  
  # Select for Group 1 crops
  radar_data_normalized <- radar_data_normalized[rownames(radar_data_normalized) %in% cropGroup, ]
  
  # Add max and min rows to create the radar plot's boundaries
  radar_data_for_plot <- rbind(
    max = rep(1, ncol(radar_data_normalized)),
    min = rep(0, ncol(radar_data_normalized)),
    radar_data_normalized
  )
  
  # Drop FARS_crop column
  radar_data_for_plot <- radar_data_for_plot[, -1]
  
  # Get the custom colors for each crop type from `crop_plot_cols`
  crop_colors <- crop_plot_cols[rownames(radar_data_normalized)]
  
  # Create the radar plot using the custom colors
  radarchart(grp1_data_for_plot, axistype = 1, pcol = grp1_crop_colors, 
             plwd = 2, plty = 1, 
             title = "Crop Type Resource Contribution Radar Plot (Group 1)",
             vlabels = colnames(grp1_data_normalized))
  
  
  # Create the radar plot using the custom colors
  radarchart(radar_data_for_plot, axistype = 1, pcol = crop_colors, 
             plwd = 2, plty = 1, 
             title = "Crop Type Resource Contribution Radar Plot",
             vlabels = colnames(radar_data_normalized))
  
  # Add a legend using the custom colors for each crop type
  legend(x = 1.5, y = 1, legend = rownames(radar_data_normalized), 
         col = crop_colors, lty = 1, lwd = 2)
}

# Apply function
getRadarCropPlot(grp_1_crops)
getRadarCropPlot(grp_2_crops)