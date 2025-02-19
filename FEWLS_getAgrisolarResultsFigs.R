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

library(patchwork)

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

# Attach area column to each dataframe from in_solar_df by Index
utilFEW_2042$Area <- in_solar_df$dir_a[match(utilFEW_2042$Index, in_solar_df$Index)]
commFEW_2042$Area <- in_solar_df$dir_a[match(commFEW_2042$Index, in_solar_df$Index)]
utilEcon_2042$Area <- in_solar_df$dir_a[match(utilEcon_2042$Index, in_solar_df$Index)]
commEcon_2042$Area <- in_solar_df$dir_a[match(commEcon_2042$Index, in_solar_df$Index)]

# Convert m2 to hectares
utilFEW_2042$Area <- utilFEW_2042$Area / 10000
commFEW_2042$Area <- commFEW_2042$Area / 10000
utilEcon_2042$Area <- utilEcon_2042$Area / 10000
commEcon_2042$Area <- commEcon_2042$Area / 10000

# Get area for two area bias estimation methods
# Area bias sensititivty analysis (Ong et al., 2013)
ongBias <- in_solar_df
ongBias$dir_a <- ifelse(ongBias$Capacity<20, ongBias$dir_a/0.71, ongBias$dir_a/0.91)
# Area bias sensititivty analysis (Fire Code Buffer)
fireCode <- in_solar_df
st_crs(fireCode) <- 4326
fireCode <- st_transform(fireCode, crs = 7801)
fireCode <- st_buffer(fireCode, dist = 7)
fireCode$dir_a <- st_area(fireCode) %>% as.numeric()
fireCode <- st_transform(fireCode, crs = 4326)
# Convert to hectares
ongBias$dir_a <- ongBias$dir_a / 10000
fireCode$dir_a <- fireCode$dir_a / 10000
# Save commercial and utility scale values for databatable figure
ongBiasComm <- ongBias$dir_a[ongBias$Index %in% commFEW_2042$Index] %>% sum()
ongBiasUtil <- ongBias$dir_a[ongBias$Index %in% utilFEW_2042$Index] %>% sum()
fireCodeComm <- fireCode$dir_a[fireCode$Index %in% commFEW_2042$Index] %>% sum()
fireCodeUtil <- fireCode$dir_a[fireCode$Index %in% utilFEW_2042$Index] %>% sum()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get FEW Graphical Table Fig

# Prep the dataframes. We need to get units in correct format, and ensure that FOOD/ENERGY 
# First, get food, energy, and water columns accounting for all input variables.
# For utility-scale
utilFEW_2042$food_base <- utilFEW_2042$kcal
utilFEW_2042$food_min <- utilFEW_2042$kcal_min
utilFEW_2042$food_max <- utilFEW_2042$kcal_max
utilFEW_2042$energy_base <- utilFEW_2042$energy_base + utilFEW_2042$irrEngy_fg
utilFEW_2042$energy_min <- utilFEW_2042$energy_min + utilFEW_2042$irrEngy_fg_wet
utilFEW_2042$energy_max <- utilFEW_2042$energy_max + utilFEW_2042$irrEngy_fg_dry
utilFEW_2042$water_base <- utilFEW_2042$irrig_fg - utilFEW_2042$oandm_base
utilFEW_2042$water_min <- utilFEW_2042$irrig_fg_wet - utilFEW_2042$oandm_max
utilFEW_2042$water_max <- utilFEW_2042$irrig_fg_dry - utilFEW_2042$oandm_min
# For commercial-scale
commFEW_2042$food_base <- commFEW_2042$kcal
commFEW_2042$food_min <- commFEW_2042$kcal_min
commFEW_2042$food_max <- commFEW_2042$kcal_max
commFEW_2042$energy_base <- commFEW_2042$energy_base + commFEW_2042$irrEngy_fg
commFEW_2042$energy_min <- commFEW_2042$energy_min + commFEW_2042$irrEngy_fg_wet
commFEW_2042$energy_max <- commFEW_2042$energy_max + commFEW_2042$irrEngy_fg_dry
commFEW_2042$water_base <- commFEW_2042$irrig_fg - commFEW_2042$oandm_base
commFEW_2042$water_min <- commFEW_2042$irrig_fg_wet - commFEW_2042$oandm_max
commFEW_2042$water_max <- commFEW_2042$irrig_fg_dry - commFEW_2042$oandm_min

# Get total budget from econ dataframes
utilFEW_2042$Tot_Budg_base <- utilEcon_2042$Tot_Budg_base[match(utilFEW_2042$Index, utilEcon_2042$Index)]
utilFEW_2042$Tot_Budg_min <- utilEcon_2042$Tot_Budg_min[match(utilFEW_2042$Index, utilEcon_2042$Index)]
utilFEW_2042$Tot_Budg_max <- utilEcon_2042$Tot_Budg_max[match(utilFEW_2042$Index, utilEcon_2042$Index)]
commFEW_2042$Tot_Budg_base <- commEcon_2042$Tot_Budg_base[match(commFEW_2042$Index, commEcon_2042$Index)]
commFEW_2042$Tot_Budg_min <- commEcon_2042$Tot_Budg_min[match(commFEW_2042$Index, commEcon_2042$Index)]
commFEW_2042$Tot_Budg_max <- commEcon_2042$Tot_Budg_max[match(commFEW_2042$Index, commEcon_2042$Index)]

# Plot colors
food_col = "wheat4"
energy_col = "orange"
water_col = "blue"
phase_alpha = 0.2
crop_plot_cols = c("grapes" = "purple3", 
                   "cotton" = "chocolate1", 
                   "crops, other" = "lightskyblue1", 
                   "grain" = "peru", 
                   "hay/pasture" = "khaki2",
                   "orchards" = "red2", 
                   "vegetables" = "yellowgreen")

# Set crop groups
crp_grp1 = "wheat|sorghum|grain|rice"
crp_grp2 = "vegetable|tomatoes|potatoes|beans|corn"
crp_grp3 = "hay|pastureland"
crp_nme1 = "grain"
crp_nme2 = "vegetables"
crp_nme3 = "hay/pasture"

# Crop order 
crpOrder = c("grain", "hay/pasture", "orchards", "grapes", "vegetables", "cotton", "crops, other")

# Get distribution by crop
getResource_CropDistFig = function(df){
  # Get df  
  crop_df <- df
  nass = Nass_Classifications
  nass$Crop <- tolower(nass$Crop)
  crop_df$FARS_crop <- nass$Irrig_Crop[match(tolower(crop_df$Crop), nass$Crop)] %>% tolower()
  # Group further to simplify
  crop_df$FARS_crop <- ifelse(grepl(crp_grp1,crop_df$FARS_crop), crp_nme1, crop_df$FARS_crop)
  crop_df$FARS_crop <- ifelse(grepl(crp_grp2,crop_df$FARS_crop), crp_nme2, crop_df$FARS_crop)
  crop_df$FARS_crop <- ifelse(grepl(crp_grp3,crop_df$FARS_crop), crp_nme3, crop_df$FARS_crop)
  # Rename 'berry totals' to grapes
  crop_df$FARS_crop <- ifelse(crop_df$FARS_crop == "berry totals", "grapes", crop_df$FARS_crop)
  # Set the factor levels for FARS_crop to desired order
  crop_df$FARS_crop  <- factor(crop_df$FARS_crop , levels = rev(crpOrder))
  return(crop_df)
}

# Apply function
utilFEW_2042 <- getResource_CropDistFig(utilFEW_2042)
commFEW_2042 <- getResource_CropDistFig(commFEW_2042)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

df = utilFEW_2042
resource_col = "Tot_Budg_base"
title = "Total Budget"
Scale = "Utility"

# Function to create horizontal stacked bar chart with independent sorting for Utility and Commercial
plot_stacked_percentage_scale <- function(df, resource_col, title, Scale) {
  
  # Summarize and calculate percentage for the given Scale (Utility or Commercial)
  df_summary <- df %>%
    group_by(FARS_crop) %>%
    summarise(Total = sum(get(resource_col))) %>%
    mutate(Percent = Total / sum(Total) * 100, Scale = Scale) %>%
    arrange(Percent) %>%
    mutate(FARS_crop = factor(FARS_crop, levels = FARS_crop))
  
  # Create the plot
  ggplot(df_summary, aes(x = Scale, y = Percent, fill = FARS_crop)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = crop_plot_cols) +
    labs(fill = "Crop", title = title) +
    theme_pubr() +
    theme(legend.position = "none", axis.title = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          panel.grid = element_blank()) + # Customize the theme
    coord_flip() # Flip the bars horizontally
}

# Function to create and combine two plots (one for Utility, one for Commercial)
plot_stacked_percentage <- function(util_df, comm_df, resource_col, title) {
  
  # Create separate plots for Utility and Commercial
  util_plot <- plot_stacked_percentage_scale(util_df, resource_col, title, "Utility")
  comm_plot <- plot_stacked_percentage_scale(comm_df, resource_col, title, "Commercial")
  
  # Combine the plots using patchwork
  combined_plot <- util_plot | comm_plot #+ plot_layout(heights = c(1, 1))
  
  return(combined_plot)
}

# Generate stacked percentage bar charts
land_plot <- plot_stacked_percentage(utilFEW_2042, commFEW_2042, "Area", "Land Use")
food_plot <- plot_stacked_percentage(utilFEW_2042, commFEW_2042, "food_base", "Food Contribution")
energy_plot <- plot_stacked_percentage(utilFEW_2042, commFEW_2042, "energy_base", "Energy Contribution")
water_plot <- plot_stacked_percentage(utilFEW_2042, commFEW_2042, "water_base", "Water Contribution")
econ_plot <- plot_stacked_percentage(utilFEW_2042, commFEW_2042, "Tot_Budg_base", "Total Budget")

# Combine the plots into a grid layout using grid.arrange
FEWplot <- land_plot / food_plot / energy_plot / water_plot / econ_plot #+ plot_layout(heights = c(1, 1, 1, 1, 1))
FEWplot

# Export
ggsave("Outputs/Figures/FEWLS_AgrisolarFEWResults.pdf", FEWplot, width = 10, height = 10, units = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set adjustment factors for reporting in units
kcalTotAdj = 1e12 # Trillion of kcal
kcalFootAdj = 1e6 # million kcal
waterTotAdj = 1e6 # million m3
waterFootAdj = 1e3 # thousand m3
energyTotAdj = 1e3 # TWh, start in GWh
energyFootAdj = 1e0 # GWh, start in GWh
econTotAdj = 1e6 # million USD
economicFootAdj = 1e3 # thousand USD

# Set footprint correction, we divide the sum of each variable by the area*years to get unit per ha per year
utilFootprintCorr = sum(utilFEW_2042$Area) * system_lifespan
commFootprintCorr = sum(commFEW_2042$Area) * system_lifespan

# Create the resources dataframe with placeholders
resources <- data.frame(
  Resource = c("Land", "Food", "Energy", "Water", "Econ"),
  Utility_Total_Impact = NA,
  Commercial_Total_Impact = NA,
  Utility_Footprint_per_ha_per_year = NA,
  Commercial_Footprint_per_ha_per_year = NA)

# Calculate Total Impact for Utility and Commercial (numeric values, no strings)
resources$Utility_Total_Impact <- c(
  sum(utilFEW_2042$Area),
  sum(utilFEW_2042$food_base) / kcalTotAdj,
  sum(utilFEW_2042$energy_base) / energyTotAdj,
  sum(utilFEW_2042$water_base) / waterTotAdj, 
  sum(utilEcon_2042$Tot_Budg_base) / econTotAdj)
resources$Commercial_Total_Impact <- c(
  sum(commFEW_2042$Area),
  sum(commFEW_2042$food_base) / kcalTotAdj,
  sum(commFEW_2042$energy_base) / energyTotAdj,
  sum(commFEW_2042$water_base) / waterTotAdj, 
  sum(commEcon_2042$Tot_Budg_base) / econTotAdj)

# Add min and max values for impact (numeric)
resources$Utility_Min_Max <- c(
  paste0("(", "NA", ", ", "NA", ")"),
  paste0("(", sum(utilFEW_2042$food_min)/kcalTotAdj, ", ", sum(utilFEW_2042$food_max)/kcalTotAdj, ")"),
  paste0("(", sum(utilFEW_2042$energy_min)/energyTotAdj, ", ", sum(utilFEW_2042$energy_max)/energyTotAdj, ")"),
  paste0("(", sum(utilFEW_2042$water_min)/waterTotAdj, ", ", sum(utilFEW_2042$water_max)/waterTotAdj, ")"), 
  paste0("(", sum(utilEcon_2042$Tot_Budg_min)/econTotAdj, ", ", sum(utilEcon_2042$Tot_Budg_max)/econTotAdj, ")"))
resources$Commercial_Min_Max <- c(
  paste0("(", "NA", ", ", "NA", ")"),
  paste0("(", sum(commFEW_2042$food_min)/kcalTotAdj, ", ", sum(commFEW_2042$food_max)/kcalTotAdj, ")"),
  paste0("(", sum(commFEW_2042$energy_min)/energyTotAdj, ", ", sum(commFEW_2042$energy_max)/energyTotAdj, ")"),
  paste0("(", sum(commFEW_2042$water_min)/waterTotAdj, ", ", sum(commFEW_2042$water_max)/waterTotAdj, ")"), 
  paste0("(", sum(commEcon_2042$Tot_Budg_min)/econTotAdj, ", ", sum(commEcon_2042$Tot_Budg_max)/econTotAdj, ")"))

# Calculate Footprint per ha per year for Utility and Commercial (numeric)
resources$Utility_Footprint_per_ha_per_year <- c(
  "NA",
  round(sum(utilFEW_2042$food_base) / utilFootprintCorr / kcalFootAdj, 8),
  round(sum(utilFEW_2042$energy_base) / utilFootprintCorr / energyFootAdj, 8),
  round(sum(utilFEW_2042$water_base) / utilFootprintCorr / waterFootAdj, 8), 
  round(sum(utilEcon_2042$Tot_Budg_base) / utilFootprintCorr / economicFootAdj, 8))

resources$Commercial_Footprint_per_ha_per_year <- c(
  "NA",
  round(sum(commFEW_2042$food_base) / commFootprintCorr / kcalFootAdj, 8),
  round(sum(commFEW_2042$energy_base) / commFootprintCorr / energyFootAdj, 8),
  round(sum(commFEW_2042$water_base) / commFootprintCorr / waterFootAdj, 8), 
  round(sum(commEcon_2042$Tot_Budg_base) / commFootprintCorr / economicFootAdj, 8))

"
# Calculate Footprint per ha per year for Utility and Commercial (numeric)
resources$Utility_Footprint_per_ha_per_year <- c(
  round(sum(utilFEW_2042$Capacity) / utilFootprintCorr / 1, 8),
  round(sum(utilFEW_2042$food_base) / utilFootprintCorr / kcalFootAdj, 8),
  round(sum(utilFEW_2042$energy_base) / utilFootprintCorr / energyFootAdj, 8),
  round(sum(utilFEW_2042$water_base) / utilFootprintCorr / waterFootAdj, 8))

resources$Commercial_Footprint_per_ha_per_year <- c(
  round(sum(commFEW_2042$Capacity) / commFootprintCorr / 1, 8),
  round(sum(commFEW_2042$food_base) / commFootprintCorr / kcalFootAdj, 8),
  round(sum(commFEW_2042$energy_base) / commFootprintCorr / energyFootAdj, 8),
  round(sum(commFEW_2042$water_base) / commFootprintCorr / waterFootAdj, 8))
"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get Economic Waterfall Plot

# Plot colors
food_col = "red" # "wheat4"
energy_col = "blue" # "orange"
water_col = "blue"
operation_col = "blue" #"brown"
install_col = "red" # "yellow"
oandm_col = "red"
total_col = "darkgreen"

# Set plot variables 
geomBarWidth = 0.75
geomBarBorderWidth = 0.50
errorBarWidth = 0.125
errorBarThickness = 0.65
hlineLineWidth = 0.5
hlineLineType = "solid"
hlineColor = 'darkgrey'
conSegLineCol = 'black'
conSegLineSize = 0.5
conSegLineType = 'solid'

# Set new econ adjust
econAdj = 1e3 # Thousand USD

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Utility scale

# Set variable order
utilVarOrder = c("Lease", "Operation", "Water", "Food", "Profit")

# Create a new dataframe for the waterfall plot (based on utilEcon_2042), creating a footprint (ha-1 yr-1)
utilEcon_waterfall <- data.frame(
  Metric = utilVarOrder,
  Base_Value = c(
                 sum(utilEcon_2042$lease_base) / utilFootprintCorr / econAdj, 
                 sum(utilEcon_2042$operation_base) / utilFootprintCorr / econAdj, 
                 sum(utilEcon_2042$water_base) / utilFootprintCorr / econAdj, 
                 -sum(utilEcon_2042$food_base) / utilFootprintCorr / econAdj, 
                 sum(utilEcon_2042$Tot_Budg_base) / utilFootprintCorr / econAdj),
  Min_Value = c(
                sum(utilEcon_2042$lease_min) / utilFootprintCorr / econAdj, 
                sum(utilEcon_2042$operation_min) / utilFootprintCorr / econAdj, 
                sum(utilEcon_2042$water_min) / utilFootprintCorr / econAdj, 
                -sum(utilEcon_2042$food_max) / utilFootprintCorr / econAdj, # Note: max and min are switched for negative variables
                sum(utilEcon_2042$Tot_Budg_min) / utilFootprintCorr / econAdj), 
  Max_Value = c(
                sum(utilEcon_2042$lease_max) / utilFootprintCorr / econAdj, 
                sum(utilEcon_2042$operation_max) / utilFootprintCorr / econAdj,
                sum(utilEcon_2042$water_max) / utilFootprintCorr / econAdj,
                -sum(utilEcon_2042$food_min) / utilFootprintCorr / econAdj,
                sum(utilEcon_2042$Tot_Budg_max) / utilFootprintCorr / econAdj))

# Round the values to zero decimal places
utilEcon_waterfall <- utilEcon_waterfall %>%
  mutate(Base_Value = round(Base_Value, 3),
         Min_Value = round(Min_Value, 3),
         Max_Value = round(Max_Value, 3))

# Order the dataframe rows by Metric = c("Food", "Lease", "Operation", "Water", "Total Budget")
utilEcon_waterfall <- utilEcon_waterfall %>% mutate(Metric = factor(Metric, levels = utilVarOrder))

# Set empty Start and End columns
utilEcon_waterfall$Start <- 0
utilEcon_waterfall$End <- 0

# For each varaible, get start and end
for(var in utilVarOrder){
  # Get the index of the variable
  idx <- which(utilEcon_waterfall$Metric == var)
  # Set the start and end values
  utilEcon_waterfall$Start[idx] <- sum(utilEcon_waterfall$Base_Value[1:idx-1])
  utilEcon_waterfall$End[idx] <- sum(utilEcon_waterfall$Base_Value[1:idx])
  # If var is Total Budget, set the bar start to zero and bar end to base value
  if(var == "Profit"){
    utilEcon_waterfall$Start[idx] <- 0
    utilEcon_waterfall$End[idx] <- utilEcon_waterfall$Base_Value[idx]
  }
  # Set the min and max values
  utilEcon_waterfall$Min_Value[idx] <- utilEcon_waterfall$End[idx] - abs(utilEcon_waterfall$Base_Value[idx] - utilEcon_waterfall$Min_Value[idx])
  utilEcon_waterfall$Max_Value[idx] <- utilEcon_waterfall$End[idx] + abs(utilEcon_waterfall$Base_Value[idx] - utilEcon_waterfall$Max_Value[idx])
}

# Create the waterfall plot
util_waterfall_plot <- ggplot(utilEcon_waterfall, aes(x = Metric)) +
  geom_hline(yintercept = 0, color = hlineColor, size = hlineLineWidth, lty = hlineLineType) +  # Zero line
  geom_rect(aes(xmin = as.numeric(Metric) - (geomBarWidth/2), xmax = as.numeric(Metric) + (geomBarWidth/2), 
                ymin = Start, ymax = End, fill = Metric), color = "black", size = geomBarBorderWidth) +
  geom_errorbar(aes(x = Metric, ymin = Min_Value, ymax = Max_Value), width = errorBarWidth, size = errorBarThickness) +
  geom_segment(data = utilEcon_waterfall[utilEcon_waterfall$Metric != "Profit", ], aes(x = as.numeric(Metric) + (geomBarWidth/2), xend = as.numeric(Metric) + 1 - (geomBarWidth/2), y = End, yend = End), 
               color = conSegLineCol, size = conSegLineSize, lty = conSegLineType) +
  # Add text to display the baseline value at the midpoint of each bar
  geom_text(aes(x = Metric, y = (Start + End) / 2, label = round(Base_Value, 2)), color = "white", size = 7/.pt, vjust = 0.5) +
  scale_fill_manual(values = c(energy_col, operation_col, water_col, food_col, total_col)) +
  scale_y_continuous(breaks = seq(-500, 500, 50), limits = c(0, 225), expand = c(0,0)) +
  #scale_y_continuous(breaks = seq(-100, 100, 2), limits = c(-0.75, 7.1), expand = c(0,0)) +
  xlab("Cash Flow Source") + 
  ylab("Economic Footprint (thousand USD ha-1 yr-1)") +
  theme_pubr() + 
  theme(legend.position = "none", axis.title = element_text(face = "bold"), panel.border = element_rect(color = "black", fill = NA, size = border_size), panel.background = element_blank()) # c(0.15, 0.7)
util_waterfall_plot

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Commercial scale

# Set variable order for Commercial scale
commVarOrder = c("NEM", "Operation", "Water", "Installation", "O&M", "Food", "Profit")

# Create a new dataframe for the waterfall plot (based on commEcon_2042), creating a footprint (ha-1 yr-1)
commEcon_waterfall <- data.frame(
  Metric = commVarOrder,
  Base_Value = c(
                 sum(commEcon_2042$energy_base) / commFootprintCorr / econAdj, 
                 sum(commEcon_2042$operation_base) / commFootprintCorr / econAdj, 
                 sum(commEcon_2042$water_base) / commFootprintCorr / econAdj, 
                 -sum(commEcon_2042$install_base) / commFootprintCorr / econAdj,
                 -sum(commEcon_2042$oandm_base) / commFootprintCorr / econAdj,
                 -sum(commEcon_2042$food_base) / commFootprintCorr / econAdj, 
                 sum(commEcon_2042$Tot_Budg_base) / commFootprintCorr / econAdj),
  Min_Value = c(
                sum(commEcon_2042$energy_min) / commFootprintCorr / econAdj, 
                sum(commEcon_2042$operation_min) / commFootprintCorr / econAdj, 
                sum(commEcon_2042$water_min) / commFootprintCorr / econAdj, 
                -sum(commEcon_2042$install_max) / commFootprintCorr / econAdj, # Note: max and min are switched for negative variables
                -sum(commEcon_2042$oandm_max) / commFootprintCorr / econAdj,
                -sum(commEcon_2042$food_max) / commFootprintCorr / econAdj, 
                sum(commEcon_2042$Tot_Budg_min) / commFootprintCorr / econAdj), 
  Max_Value = c(
                sum(commEcon_2042$energy_max) / commFootprintCorr / econAdj, 
                sum(commEcon_2042$operation_max) / commFootprintCorr / econAdj,
                sum(commEcon_2042$water_max) / commFootprintCorr / econAdj,
                -sum(commEcon_2042$install_min) / commFootprintCorr / econAdj,
                -sum(commEcon_2042$oandm_min) / commFootprintCorr / econAdj,
                -sum(commEcon_2042$food_min) / commFootprintCorr / econAdj,
                sum(commEcon_2042$Tot_Budg_max) / commFootprintCorr / econAdj))

# Round the values to zero decimal places
commEcon_waterfall <- commEcon_waterfall %>%
  mutate(Base_Value = round(Base_Value, 3),
         Min_Value = round(Min_Value, 3),
         Max_Value = round(Max_Value, 3))

# Order the dataframe rows by Metric = c("Food Loss", "Energy Return", "Operation Offset", "Water Saved", "Total Budget")
commEcon_waterfall <- commEcon_waterfall %>% mutate(Metric = factor(Metric, levels = commVarOrder))

# Set empty Start and End columns
commEcon_waterfall$Start <- 0
commEcon_waterfall$End <- 0

# For each variable, get start and end
for(var in commVarOrder){
  # Get the index of the variable
  idx <- which(commEcon_waterfall$Metric == var)
  # Set the start and end values
  commEcon_waterfall$Start[idx] <- sum(commEcon_waterfall$Base_Value[1:idx-1])
  commEcon_waterfall$End[idx] <- sum(commEcon_waterfall$Base_Value[1:idx])
  # If var is Total Budget, set the bar start to zero and bar end to base value
  if(var == "Profit"){
    commEcon_waterfall$Start[idx] <- 0
    commEcon_waterfall$End[idx] <- commEcon_waterfall$Base_Value[idx]
  }
  # Set the min and max values
  commEcon_waterfall$Min_Value[idx] <- commEcon_waterfall$End[idx] - abs(commEcon_waterfall$Base_Value[idx] - commEcon_waterfall$Min_Value[idx])
  commEcon_waterfall$Max_Value[idx] <- commEcon_waterfall$End[idx] + abs(commEcon_waterfall$Base_Value[idx] - commEcon_waterfall$Max_Value[idx])
}

# Create the waterfall plot for Commercial scale
comm_waterfall_plot <- ggplot(commEcon_waterfall, aes(x = Metric)) +
  geom_hline(yintercept = 0, color = hlineColor, size = hlineLineWidth, lty = hlineLineType) +  # Zero line
  geom_rect(aes(xmin = as.numeric(Metric) - (geomBarWidth/2), xmax = as.numeric(Metric) + (geomBarWidth/2), 
                ymin = Start, ymax = End, fill = Metric), color = "black", size = geomBarBorderWidth) +
  geom_errorbar(aes(x = Metric, ymin = Min_Value, ymax = Max_Value), width = errorBarWidth, size = errorBarThickness) +
  geom_segment(data = commEcon_waterfall[commEcon_waterfall$Metric != "Profit", ], aes(x = as.numeric(Metric) + (geomBarWidth/2), xend = as.numeric(Metric) + 1 - (geomBarWidth/2), y = End, yend = End), 
               color = conSegLineCol, size = conSegLineSize, lty = conSegLineType) +
  # Add text to display the baseline value at the midpoint of each bar
  geom_text(aes(x = Metric, y = (Start + End) / 2, label = round(Base_Value, 2)), color = "white", size = 7/.pt, vjust = 0.5) +
  scale_fill_manual(values = c(energy_col, operation_col, water_col, install_col, oandm_col, food_col, total_col)) +
  scale_y_continuous(breaks = seq(-500, 500, 50), limits = c(0, 225), expand = c(0,0)) +
  xlab("Cash Flow Source") + 
  ylab("Economic Footprint (thousand USD ha-1 yr-1)") +
  theme_pubr() + 
  theme(legend.position = "none", axis.title = element_text(face = "bold"), panel.border = element_rect(color = "black", fill = NA, size = border_size), panel.background = element_blank()) 
#comm_waterfall_plot

# ~~~~~~~~~~~~~~~~~~~~~~~~~

# Rotate the x-axis text for both plots 
util_waterfall_plot <- util_waterfall_plot + theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
comm_waterfall_plot <- comm_waterfall_plot + theme(axis.text.x = element_text(angle = 30, hjust = 1)) 

# Combine the two waterfall plots and add a common x-axis label using plot_layout
combined_plot <- (comm_waterfall_plot + util_waterfall_plot) + plot_layout(widths = c(7,5)) & theme(axis.title.x = element_text(face = "bold", size = 12)) + plot_annotation(tag_levels = 'A')

# Display the combined plot
combined_plot
# Save the plot
ggsave("Outputs/Figures/FEWLS_AgrisolarEconResults.pdf", combined_plot, width = 7, height = 3.75, units = "in", dpi = 300)

