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
of CDL crop operational costs per m2 from UC Davis Crop Crop Cost Studies. Only needs to 
be run if cleanFEWLS() is run, or if completed reports from 
*https://coststudies.ucdavis.edu/en/current/* are updated within 
*Data/Downloaded/UCDavis_CropExt*.

Note: Here, we remove the cost to irrigate, because we calculate that based on electricity 
consumption (the reports do as well). However, we retain irrigation labor, because that is
maintained regarldless of the source of electricity. 
"
## ----------------------------------------------------------- ##
##                                                             ##
##    UC Davis Extension Crop Economic Values of Operation     ##
##                                                             ##
## ----------------------------------------------------------- ##

#install.packages("pdftools")
library(pdftools)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get Cost per acre production and year of study

# Get pdf locations
UCDavis_CropExt_path <- "Data/Downloaded/UCDavis_CropExt/"
cropbudg_locs <- list.files(path = UCDavis_CropExt_path, pattern = ".pdf", full.names = FALSE)

# For every file name in UC Davis Crop Extension Repository, get total operational production costs without irrigation application costs
cost_acre_df <- data.frame(Crop = character(), cost_acre = numeric(), Year = numeric())
for(cropbudg_loc in cropbudg_locs){
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get total cost per acre page number and year of report
  
  # Get crop pdf
  crop <- paste(UCDavis_CropExt_path, cropbudg_loc, sep = "") %>% pdf_text()
  # Get table of contents and split
  tbl <- crop[2:3] %>% strsplit("\n") %>% unlist()
  # Get location of table with "Costs per Acre to Produce" this crop
  cpa_tbl_loc <- tbl[grepl("cost per acre to produce|costs per acre to produce|costs and returns per 1,000lf to establish", tolower(tbl))][1]
  #cpa_tbl_loc <- tbl[grepl("Cost Per Acre to Produce|Costs per acre to Produce|Costs per Acre to Produce|COSTS PER ACRE to PRODUCE|COSTS PER ACRE TO PRODUCE|COSTS AND RETURNS PER 1,000LF TO ESTABLISH", tbl)][1]
  pg_num <- regmatches(cpa_tbl_loc, gregexpr("[[:digit:]]+", cpa_tbl_loc)) %>% unlist() %>% as.numeric() %>% tail(1)# %>% max()
  cpa_tbl_loc_CANE <- tbl[grepl("Spring and Fall Crop Production", tbl)]
  pg_numCANE <- regmatches(cpa_tbl_loc_CANE, gregexpr("[[:digit:]]+", cpa_tbl_loc_CANE)) %>% unlist() %>% as.numeric() %>% tail(1)
  pg_num <- ifelse(length(pg_num)==0, pg_numCANE, pg_num)
  
  # Get report year
  yr_lineHead <- head(crop) %>% strsplit("\n") %>% unlist() %>% head(10) # Number of lines on title page
  yr_lineTail <- tail(crop) %>% strsplit("\n") %>% unlist() %>% tail(1) # last line of pdf is report reference
  yearsHead <- regmatches(yr_lineHead, gregexpr("[[:digit:]]+", yr_lineHead)) %>% unlist() %>% as.numeric() # Gets all numerical strings in larger list of strings - Check title
  yearsTail <- regmatches(yr_lineTail, gregexpr("[[:digit:]]+", yr_lineTail)) %>% unlist() %>% as.numeric() # Gets all numerical strings in larger list of strings - Check bottom page reference
  year <- ifelse(length(yearsHead)==0, yearsTail[which(yearsTail>=1980 & yearsTail<=2100)], yearsHead[which(yearsHead>=1980 & yearsHead<=2100)]) # If title does not contain year, check last page citations
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get Total cost per acre and subtract irrigation costs (we independently estimate these)
  
  # Get Cost per acre to produce page
  tbl <- crop[pg_num:(pg_num+1)] %>% strsplit("\n") %>% unlist() # pg_num plus 1 gets rid of issue with cont'd tables
  
  # Get irrig cost
  irr_tbl_loc <- tbl[grepl("Irrigation|Irrigate", tbl)]
  irr_tbl_loc <- irr_tbl_loc[!grepl("Labor", irr_tbl_loc)] # Keep labor, we do not account for this in our estimation
  df <- data.frame(category = numeric(), cost = numeric())
  if(length(irr_tbl_loc)!=0){
    for(i in c(1:length(irr_tbl_loc))){
      cost <- irr_tbl_loc[i] %>% str_sub(-5)
      cost <- gsub(",", "", cost)
      cost <- regmatches(cost, gregexpr("[[:digit:]]+", cost)) %>% unlist() %>% as.numeric() %>% max()
      costdf <- data.frame(category = i, cost = cost)
      df <- rbind(df, costdf)
    }
  }
  # Remove NA, Inf, and -Inf from df
  df <- df[which(!df$cost %in% c(NA, Inf, -Inf)), ]
  
  # Get get total per acre cost and subtract irrgation costs (not labor)
  tot_tbl_loc <- tbl[grepl("TOTALCOSTS/ACRE|TOTAL COSTS/ACRE|TOTAL OPERATING COSTS/ACRE|TOTAL COST/1,000LF", tbl)] # Get total cost from table
  cost <- tot_tbl_loc %>% str_sub(-6) # 6 characters allows for ten-thousand's of dollars
  cost <- gsub(",", "", cost)
  cost <- regmatches(cost, gregexpr("[[:digit:]]+", cost)) %>% unlist() %>% as.numeric() %>% max() # Get end value in table which is number, and get maximum total cost (some are preliminary)
  cost <- ifelse(length(tot_tbl_loc[grepl('TOTAL COST/1,000LF', tot_tbl_loc)])>0, cost/(10000/43560), cost) # If caneberry (measured in 1000 linear feet), convert to acre
  crop_cost_acre <- cost - sum(df$cost, na.rm = TRUE)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save to df with crop name and get next crop pdf
  
  # Get crop name
  crop_name <- sub("_.*", "", cropbudg_loc)
  crop_df <- data.frame(Crop = crop_name, cost_acre = crop_cost_acre, Year = year)
  cost_acre_df <- rbind(cost_acre_df, crop_df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Adjust for inflation and convert to $/m2

# Multiply cost_acre_df by inflations adjustment and group together -- PPI because farmers are producers and this data is "price recieved"
cost_acre_df$cost_acre <- cost_acre_df$cost_acre * inflation_ratePPI$infl_rate_mult[match(cost_acre_df$Year, inflation_ratePPI$Year)]

# Group by crop and convert to meter
#m2_acre = 4046.86 # m2/acre
cost_m2_df <- cost_acre_df %>% group_by(Crop) %>% summarize(cost_m2 = round(mean(cost_acre, na.rm=TRUE)/m2_acre, 4))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Compare and merge with NASS CDL crop classes

# Get nass crop classes
nass_temp <- Nass_Classifications 

# Get cost/m2 for California crops
df <- nass_temp[, c('Value', 'Crop', 'Irrig_Crop')]
df$Crop <- tolower(df$Crop)
df$Irrig_Crop <- tolower(df$Irrig_Crop)
df$cost_m2 <- 0
# Check if CDL crop compares
for(i in c(1:nrow(cost_m2_df))){
  # Get forgone crop type and year from survey data
  crop <- cost_m2_df[i, ]$Crop %>% as.character()
  cost <- cost_m2_df[i, ]$cost_m2 %>% as.numeric()
  # Get matching crop type from dataset
  if(length(df[which(grepl(crop, df$Crop)), ]$cost_m2)!=0){
    df[which(grepl(crop, df$Crop)), ]$cost_m2 <- cost}
}
# Check if FARS crop compares
for(i in c(1:nrow(cost_m2_df))){
  # Get forgone crop type and year from survey data
  crop <- cost_m2_df[i, ]$Crop %>% as.character()
  cost <- cost_m2_df[i, ]$cost_m2 %>% as.numeric()
  
  # Get matching crop type from dataset
  df[which(grepl(crop, df$Irrig_Crop)), ]$cost_m2 <- ifelse(df[which(grepl(crop, df$Irrig_Crop)), ]$cost_m2==0 & length(df[which(grepl(crop, df$Irrig_Crop)), ]$cost_m2)!=0, cost, df[which(grepl(crop, df$Irrig_Crop)), ]$cost_m2)
}

# Get FARS crop classes which contain a value, and elicit that value for all FARS crops without value in the same crop classes
df_grouped <- df %>% group_by(Irrig_Crop) %>% summarize(cost_m2 = max(cost_m2, na.rm = TRUE))
df_temp <- df[which(df$cost_m2==0), ]
df <- df[which(!df$cost_m2==0), ]
for(i in c(1:nrow(df_grouped))){
  # Get forgone crop type and year from survey data
  crop <- df_grouped[i, ]$Irrig_Crop %>% as.character()
  cost <- df_grouped[i, ]$cost_m2 %>% as.numeric()
  # Get matching crop type from dataset
  if(length(df_temp[which(grepl(crop, df_temp$Irrig_Crop)), ]$cost_m2)!=0){
    df_temp[which(grepl(crop, df_temp$Irrig_Crop)), ]$cost_m2 <- cost}
}
# Recombine
df <- rbind(df, df_temp)

# Remove previous erroneous double crops (artifacts from method described above) and add correct double crop sums
df <- df[which(!grepl("dbl", df$Crop)), ]
df$Irrig_Crop = NULL
df$Value = NULL
dbl_crop_func <- function(crop_1, crop_2){
  # Add double crop row by sum of irrig for both crops for both years of the study
  df <- add_row(df, 
                Crop = paste("dbl crop ", crop_1, "/", crop_2, sep = ""), 
                cost_m2 = as.numeric(df[which(df$Crop==crop_1), ]['cost_m2']) + as.numeric(df[which(df$Crop==crop_2), ]['cost_m2']))
  return(df) }
# Run doule crop function
df <- dbl_crop_func("barley", "corn")
df <- dbl_crop_func("barley", "sorghum")
df <- dbl_crop_func("barley", "soybeans")
df <- dbl_crop_func("corn", "soybeans")
df <- dbl_crop_func("durum wheat", "sorghum")
df <- dbl_crop_func("lettuce", "barley")
df <- dbl_crop_func("lettuce", "cantaloupes")
df <- dbl_crop_func("lettuce", "cotton")
df <- dbl_crop_func("lettuce", "durum wheat")
df <- dbl_crop_func("oats", "corn")
df <- dbl_crop_func("soybeans", "cotton")
df <- dbl_crop_func("soybeans", "oats")
df <- dbl_crop_func("winter wheat", "corn")
df <- dbl_crop_func("winter wheat", "cotton")
df <- dbl_crop_func("winter wheat", "sorghum")
df <- dbl_crop_func("winter wheat", "soybeans")

# Fix names of appreviated crops
df$Crop <- ifelse(df$Crop=="dbl crop durum wheat/sorghum", "dbl crop durum wht/sorghum", df$Crop)
df$Crop <- ifelse(df$Crop=="dbl crop lettuce/durum wheat", "dbl crop lettuce/durum wht", df$Crop)
df$Crop <- ifelse(df$Crop=="dbl crop lettuce/cantaloupes", "dbl crop lettuce/cantaloupe", df$Crop)
df$Crop <- ifelse(df$Crop=="dbl crop winter wheat/corn", "dbl crop winwht/corn", df$Crop)
df$Crop <- ifelse(df$Crop=="dbl crop winter wheat/cotton", "dbl crop winwht/cotton", df$Crop)
df$Crop <- ifelse(df$Crop=="dbl crop winter wheat/sorghum", "dbl crop winwht/sorghum", df$Crop)
df$Crop <- ifelse(df$Crop=="dbl crop winter wheat/soybeans", "dbl crop winwht/soybeans", df$Crop)

# Make assumptions about remaining missing crops
df$cost_m2[which(df$cost_m2==0)] <- median(cost_m2_df$cost_m2)
nass_temp$Crop <- tolower(nass_temp$Crop)
nass_temp$cost_m2 <- df$cost_m2[match(nass_temp$Crop, df$Crop)]
nass_temp$cost_m2 <- ifelse(is.na(nass_temp$Irrig_Crop), 0, nass_temp$cost_m2)
df <- nass_temp

# Export 
df <- df[, c("Crop", "cost_m2")]
df$Crop <- tolower(df$Crop)
write.csv(df, "Data/Derived/operational_cost_m2.csv", row.names=FALSE)