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
of CDL crop irrigation depth (mm) from the 2013 Farm and Ranch Survey and the 2018 Irrigation
and Water Management Survey. Only needs to be run if cleanFEWLS() is run, or if surveys are 
updated. 
"
## ----------------------------------------------------------- ##
##                                                             ##
##             Build Irrigation Depth Dataframe                ##
##                                                             ##
## ----------------------------------------------------------- ##

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ------------------------ #
# Get Water use by Crop df # 
# ------------------------ #
# Change column classes
irrig <- irrig_perCrop_state[which(tolower(irrig_perCrop_state$State)==tolower(state)), ]
irrig$Value <- as.numeric(gsub(",","",irrig$Value))
# Select only weighted irrigation tech total
irrig <- irrig[which(irrig$Domain=="TOTAL"), ]
irrig <- irrig[!is.na(irrig$Value), ]
colnames(irrig)[colnames(irrig)=="Commodity"] <- "Crop"
# Get unit column
irrig$Units <- sub(".*IN", "", irrig$Data.Item)
irrig$Units <- sub("/.*", "", irrig$Units)
# Correct Horticulture totals that are partially in gallons
irrig$Value <- ifelse(irrig$Units %in% c("GALLONS", " GALLONS", "GALLONS ", " GALLONS "), irrig$Value/gal_acft, irrig$Value) # horticulture values still vast overestimates, remove gallon reproted values
irrig <- irrig[!irrig$Units %in% c("GALLONS", " GALLONS", "GALLONS ", " GALLONS "), ]
# Group by commodity title
irrig <- irrig %>%
  group_by(Crop, Year) %>%
  mutate() %>%
  summarise(Value = mean(Value), Units = max(Units))
# Irrigation from all CONUS. This is to fill in any missing crops in the California data by using US total
irrig_CONUS <- irrig_perCrop_conus[which(irrig_perCrop_conus$State=="US TOTAL"), ]
# Change column classes, Get unit column, Group by commodity title, Add missing crops from the rest of CONUS
irrig_CONUS$Value <- as.numeric(gsub(",","",irrig_CONUS$Value))
# Select only weighted irrigation tech total
irrig_CONUS <- irrig_CONUS[which(irrig_CONUS$Domain=="TOTAL"), ]
irrig_CONUS <- irrig_CONUS[!is.na(irrig_CONUS$Value), ]
colnames(irrig_CONUS)[colnames(irrig_CONUS)=="Commodity"] <- "Crop"
irrig_CONUS$Units <- sub(".*IN", "", irrig_CONUS$Data.Item)
irrig_CONUS$Units <- sub("/.*", "", irrig_CONUS$Units)
# Correct Horticulture totals that are partially in gallons
irrig_CONUS$Value <- ifelse(irrig_CONUS$Units %in% c("GALLONS", " GALLONS", "GALLONS ", " GALLONS "), irrig_CONUS$Value/gal_acft, irrig_CONUS$Value) # horticulture values still vast overestimates, remove gallon reproted values
irrig_CONUS <- irrig_CONUS[!irrig_CONUS$Units %in% c("GALLONS", " GALLONS", "GALLONS ", " GALLONS "), ]
# Group by commodity title and year
irrig_CONUS <- irrig_CONUS %>%
  group_by(Crop, Year) %>%
  mutate() %>%
  summarise(Value = mean(Value), Units = max(Units))
# Combine CONUS irrigation dataframe to get missing crop irrigation values in California
irrig <- merge(irrig, irrig_CONUS, by = c("Crop", "Year"), all = T)
irrig$Value.x <- ifelse(is.na(irrig$Value.x), irrig$Value.y, irrig$Value.x)
irrig$Units.x = NULL
irrig$Value.y = NULL
colnames(irrig)[colnames(irrig)=="Units.y"] <- "Units"
colnames(irrig)[colnames(irrig)=="Value.x"] <- "Value"
# Change factor of crop to all lower case to merge data frames by crop type. Remove horticulture totals because of lack of consistent data
irrig$Crop <- tolower(irrig$Crop)
rm(irrig_CONUS)
# Horticulture totals were only measured in 2018, for now assume they deviated from 2018 based on similar crop to irrigate, wheat (same irrig depth in 2018) -- Use wheat as proxy if missing hort-total
if(nrow(irrig[which(irrig$Crop=="horticulture totals" & irrig$Year==2013), ])==0){
  irrig <- add_row(irrig, Crop = "horticulture totals", Year = 2013, Value = irrig$Value[which(irrig$Crop=="wheat" & irrig$Year==2013)], Units = "ACRE FEET") }
if(nrow(irrig[which(irrig$Crop=="horticulture totals" & irrig$Year==2018), ])==0){
  irrig <- add_row(irrig, Crop = "horticulture totals", Year = 2018, Value = irrig$Value[which(irrig$Crop=="wheat" & irrig$Year==2018)], Units = "ACRE FEET") }
# Add double crop sums
dbl_crop_func <- function(crop_1, crop_2){
  # Two Farm and Ranch Survey Years
  fars1 <- 2013
  fars2 <- 2018
  # Add double crop row by sum of irrig for both crops for both years of the study
  irrig <- add_row(irrig, Year = fars1, Crop = paste(crop_1, crop_2, sep = "_"), 
                   Value = (as.numeric(irrig[which(irrig$Crop==crop_1 & irrig$Year==fars1), ]['Value']) + 
                              as.numeric(irrig[which(irrig$Crop==crop_2 & irrig$Year==fars1), ]['Value'])), 
                   Units = "ACRE FEET")
  irrig <- add_row(irrig, Year = fars2, Crop = paste(crop_1, crop_2, sep = "_"), 
                   Value = (as.numeric(irrig[which(irrig$Crop==crop_1 & irrig$Year==fars2), ]['Value']) + 
                              as.numeric(irrig[which(irrig$Crop==crop_2 & irrig$Year==fars2), ]['Value'])), 
                   Units = "ACRE FEET")
  return(irrig)
}
# Run doule crop function
irrig <- dbl_crop_func("wheat", "soybeans")
irrig <- dbl_crop_func("wheat", "corn")
irrig <- dbl_crop_func("wheat", "sorghum")
irrig <- dbl_crop_func("wheat", "cotton")
irrig <- dbl_crop_func("small grains", "corn")
irrig <- dbl_crop_func("vegetable totals", "wheat")
irrig <- dbl_crop_func("vegetable totals", "crops, other")
irrig <- dbl_crop_func("vegetable totals", "cotton")
irrig <- dbl_crop_func("vegetable totals", "hay & haylage")
irrig <- dbl_crop_func("hay & haylage", "sorghum")
irrig <- dbl_crop_func("hay & haylage", "corn")
irrig <- dbl_crop_func("soybeans", "cotton")
irrig <- dbl_crop_func("soybeans", "small grains")
irrig <- dbl_crop_func("hay & haylage", "soybeans")
irrig <- dbl_crop_func("corn", "soybeans")
irrig$Crop <- as.factor(irrig$Crop)
irrig$Crop <- tolower(irrig$Crop)
# Change irrig value to mm
irrig$Value <- irrig$Value*mm_acft_acre
irrig$Units = NULL

# Export df
write.csv(irrig, "Data/Derived/mm_per_yearIrrig.csv", row.names=FALSE)