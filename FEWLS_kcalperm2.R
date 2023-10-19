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
of spatial caloric density (kcal/m2) from NASS Agricultual Yield Surveys and USDA FoodData.
Only needs to be run if cleanFEWLS() is run, or if surveys are updated. 
"
## ----------------------------------------------------------- ##
##                                                             ##
##           Build kcal/m2 dataframe for a region              ##
##                                                             ##
## ----------------------------------------------------------- ##

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get crop names from NASS CDL
nass_temp <- Nass_Classifications
crop_names <- nass_temp[which(!nass_temp$Value %in% cdl_non_ag), ]
crop_names <- crop_names$Crop %>% tolower()
crop_names <- crop_names[!grepl("dbl", crop_names)]
# Create single crop name string separated by or operands "|" 
string = as.character("")
for(crop in crop_names){
  string = paste(string, "|", crop, sep = "")  }
string <- substring(string, 2) %>% tolower() # removes | at start of string
# Create extended single crop name string separated by or operands "|" -- this is for the second iteration of the dataset
crop_names_ext <- scan(text = crop_names, what = " ") 
string_ext = as.character("")
for(crop in crop_names_ext){
  string_ext = paste(string_ext, "|", crop, sep = "")  }
string_ext <- substring(string_ext, 2) %>% tolower() # removes | at start of string

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USDA FoodData

#  id	              name	              unit_name nutrient_nbr	rank
# 2047	Energy (Atwater General Factors)	KCAL	      957	      280
# 2048	Energy (Atwater Specific Factors)	KCAL	      958	      290
# 1008	Energy	                          KCAL	      208	      300

# Call in FNDDS nutreint value database # -- kcal values are per 100 grams of food
food_df <- read.csv(FoodData_food_loc, stringsAsFactors = FALSE)
food_nutrient_df <- read.csv(FoodData_nutrient_loc, stringsAsFactors = FALSE) 

# Subset for kcal
id_list <- c(2047, 2048, 1008)
food_nutrient_df <- food_nutrient_df[which(food_nutrient_df$nutrient_id %in% id_list), ]

# Match food desciptions by fdc_id
food_nutrient_df$crop <- food_df$description[match(food_nutrient_df$fdc_id, food_df$fdc_id)] %>% tolower()
save <- food_nutrient_df

# Create non strong -- all altered crop values in some way
non_string = "cooked|baked|roasted|salted|grilled|flavored|hash|patty|strained|baby|added|fast food|blend|snacks|puffed|boil|toasted|frosted|infused|extra|boiling|chili lime|sweetened|juice|girl scouts|frozen|covered|soup|bread|fried|with|creamed|stuffed|coated|beef|dumpling|meal|chips|cake|paper|syrup|sauce|milk|pudding|sandwhich|omlet|fritter|beverage|canned|marachino|pop|greens|mixed salad|roll|pastry|pork|acorn|chocolate|cocoa|freeze|frosted|crunch|grapeseed|quaker|grits|puffs|flakes|cheerios|kellogg's|chick"

# Get foods which match cdl crop type strings and remove undesired crops. If we miss some, we will go back though after
food_nutrient_df <- food_nutrient_df[grepl(string, food_nutrient_df$crop), ]
food_nutrient_df <- food_nutrient_df[grepl("raw|nfs|oil|nuts|unsalted|winter|cereal|vegetable|fruit", food_nutrient_df$crop), ]
food_nutrient_df <- food_nutrient_df[!grepl(non_string, food_nutrient_df$crop), ]
food_nutrient_df$crop_name <- str_extract_all(food_nutrient_df$crop, string, simplify = TRUE)

# Get only single match rows
food_nutrient_df$rm <- NA
for(i in c(1:nrow(food_nutrient_df))){
  food_nutrient_df[i, ]$rm <- ifelse(food_nutrient_df[i, ]$crop_name[2]=="", 0, 1) }
food_nutrient_df <- food_nutrient_df[which(food_nutrient_df$rm==0), ]

# Get primary crop name
food_nutrient_df$crop_name_1 = NA
for(i in c(1:nrow(food_nutrient_df))){
  food_nutrient_df[i, ]$crop_name_1 <- food_nutrient_df[i, ]$crop_name[1] }
food_nutrient_df$crop_name = NULL

# Get oil crops
oil_crops <- food_nutrient_df[grepl("oil", food_nutrient_df$crop), ]
non_oil_crops <- food_nutrient_df[!grepl("oil", food_nutrient_df$crop), ]

# Group, average, and bring back together
seedoil_cdl_crops_temp <- seedoil_cdl_crops %>% tolower()
oil_crops <- oil_crops[which(oil_crops$crop_name_1 %in% seedoil_cdl_crops_temp), ] %>% group_by(crop_name_1) %>% summarise(kcal = mean(amount))
non_oil_crops <- non_oil_crops %>% group_by(crop_name_1) %>% summarise(kcal = mean(amount))
food_nutrient_df <- rbind(non_oil_crops, oil_crops)
colnames(food_nutrient_df)[colnames(food_nutrient_df)=="crop_name_1"] <- "crop"
df <- food_nutrient_df

# Group to reduce any crops that show up in both oil and non oil crops
df <- df %>% group_by(crop) %>% summarize(kcal = mean(kcal))

# Get still missing crops to focus on -- exclude double cropped crops, we will include those later
missing_crops <- setdiff(crop_names, df$crop)
missing_crops <- missing_crops[!grepl("dbl", missing_crops)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get extra missing crops which were initially filtered out

# Make new dataframe of all missing crops which may show up in dataset but were filtered due to associated "altered" non_string, or CDL dual name like "Winter Wheat"
missing_crops_ext <- scan(text = missing_crops, what = " ") # lengthened string which separates double worded crops like "Winter Wheat" to look for "Wheat"
missing_df <- data.frame(crop = character(), kcal = numeric()) # empty dataframe to fill
still_missing_crops = character() # crops that are still missing
extra_string_non_crop = character() # crops that are not crops but the artifact of the "Winter" & "Wheat" split -- winter is not a desired crop
for(crop in missing_crops_ext){
  # Get all FoodData which contains missing crop
  temp_df <- save[grepl(crop, save$crop), ]
  temp_df$crop_name <- str_extract_all(temp_df$crop, string_ext, simplify = TRUE)
  # If crop absent from FoodData, stop, else continue
  if(nrow(temp_df)==0){
    still_missing_crops <- paste(still_missing_crops, ", ", crop, sep = "")
  } else {
    # Check length of resulting matrix, only want simple crops, not mixes
    if(dim(temp_df$crop_name)[2]>1){
      # Loop to get only data which missing crop is primary descriptor, and only descriptor
      temp_df$rm = 0
      for(i in c(1:nrow(temp_df))){
        temp_df[i, ]$rm <- ifelse(temp_df[i, ]$crop_name[2]=="" & temp_df[i, ]$crop_name[1]==crop, 0, 1) }
      temp_df <- temp_df[which(temp_df$rm==0), ]
    } else {
      # Loop to get only data which missing crop is primary descriptor, and only descriptor
      temp_df$rm = 0
      for(i in c(1:nrow(temp_df))){
        temp_df[i, ]$rm <- ifelse(temp_df[i, ]$crop_name[1]==crop, 0, 1) }
      temp_df <- temp_df[which(temp_df$rm==0), ]
    }
    # If crop absent from FoodData, stop, else continue
    if(nrow(temp_df)==0){
      extra_string_non_crop <- paste(extra_string_non_crop, ", ", crop, sep = "")
    } else {
      # Get primary crop name
      temp_df$crop_name_1 = NA
      for(i in c(1:nrow(temp_df))){
        temp_df[i, ]$crop_name_1 <- temp_df[i, ]$crop_name[1] }
      temp_df$crop_name = NULL
      # Filter crop alteration words, and group
      temp_df <- temp_df[!grepl(non_string, temp_df$crop), ]
      temp_df <- temp_df %>% group_by(crop_name_1) %>% summarise(kcal = mean(amount, na.rm=TRUE))
      temp_df$crop <- temp_df$crop_name_1
      temp_df$crop_name_1 = NULL
      # Add avearage kcal to missing_df and crop name
      missing_df <- rbind(missing_df, temp_df)
    }
  }
} # For some reason this always reports and error, but it seems to be irrelevant

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Manual additions & equivalents

# Manually remove irrelevant artifacts from extention of missing crops string, and resulting duplicates, and duplicates from already derived crops
irrel_list = c("or", "other", "&", "tree", "seed", "sweet", "spring", "dry", "rape", "trees", "small", "winter", "grains", "herbs") # Grains & herbs is included here as this is too general of a term for FoodData
missing_df <- missing_df[which(!missing_df$crop %in% irrel_list & !missing_df$crop %in% df$crop), ]
missing_df <- missing_df %>% group_by(crop) %>% summarise(kcal = mean(kcal))
df <- rbind(df, missing_df)

# Make adjustments manually for certain crop names that were split
df$crop <- ifelse(df$crop=="durum", "durum wheat", df$crop)
df$crop <- ifelse(df$crop=="wheat", "winter wheat", df$crop)
df$crop <- ifelse(df$crop=="melons", "watermelons", df$crop)
df$crop <- ifelse(df$crop=="beans", "dry beans", df$crop) 
df$crop <- ifelse(df$crop=="fruits", "misc vegs & fruits", df$crop)
df$kcal <- ifelse(df$crop=="misc vegs & fruits", mean(df$kcal[which(df$crop %in% c("misc vegs & fruits", "broccoli"))]), df$kcal)

# Manual assumed equivalent comparisons kcal/kg
equiv_crop_func <- function(new_crop, equiv_crop){
  # Add double crop row by sum of irrig for both crops for both years of the study
  df <- add_row(df, crop = new_crop, kcal = df$kcal[which(df$crop==equiv_crop)])
  return(df) }

# Add equivalent crops (derived via logic, or a source)
df <- equiv_crop_func("spring wheat", "winter wheat") # Logic
df <- equiv_crop_func("other small grains", "winter wheat") # Logic
df <- equiv_crop_func("speltz", "oats") # https://hort.purdue.edu/newcrop/afcm/spelt.html
df <- equiv_crop_func("cantaloupes", "watermelons") # Logic
df <- equiv_crop_func("honeydew melons", "watermelons") # Logic
df <- equiv_crop_func("other hay/non alfalfa", "alfalfa") # Logic
df <- equiv_crop_func("grassland/pasture", "alfalfa") # Logic
df <- equiv_crop_func("clover/wildflowers", "alfalfa") # Logic
df <- equiv_crop_func("vetch", "alfalfa") # Logic
df <- equiv_crop_func("switchgrass", "alfalfa") # Logic
df <- equiv_crop_func("camelina", "canola") # Logic & https://catalog.extension.oregonstate.edu/sites/catalog/files/project/pdf/em8953.pdf
df <- equiv_crop_func("sugarbeets", "potatoes") # Logic
df <- equiv_crop_func("sweet corn", "corn") # Logic
df <- equiv_crop_func("pop or orn corn", "corn") # Logic
df <- equiv_crop_func("rape seed", "canola") # Logic and https://www.agmrc.org/commodities-products/grains-oilseeds/rapeseed

# Manual additions of kcal/kg -- include sources
df <- add_row(df, crop = "tobacco", kcal = 0) # Tobacco not consumed for calories
df <- add_row(df, crop = "sod/grass seed", kcal = 0) # sod grass not consumed
df <- add_row(df, crop = "christmas trees", kcal = 0) # christams trees not consumed
df <- add_row(df, crop = "sweet potatoes", kcal = 109) # https://fdc.nal.usda.gov/fdc-app.html#/food-details/1103233/nutrients
df <- add_row(df, crop = "chick peas", kcal = 210) # https://fdc.nal.usda.gov/fdc-app.html#/food-details/1100429/nutrients
df <- add_row(df, crop = "eggplants", kcal = 25) # https://fdc.nal.usda.gov/fdc-app.html#/food-details/1103353/nutrients
df <- add_row(df, crop = "gourds", kcal = 14) # https://fdc.nal.usda.gov/fdc-app.html#/food-details/169232/nutrients
df <- add_row(df, crop = "herbs", kcal = 23) # basil -- https://fdc.nal.usda.gov/fdc-app.html#/food-details/172232/nutrients
df <- add_row(df, crop = "caneberries", kcal = 57) # Raspberries -- https://fdc.nal.usda.gov/fdc-app.html#/food-details/2263888/nutrients
df <- add_row(df, crop = "other tree crops", kcal = mean(df$kcal[which(df$crop %in% c("almonds","citrus","cherries","nectarines","peaches"))])) #80%fruits,20% nuts--"FSA to CDL Crosswalk" https://www.nass.usda.gov/Research_and_Science/Cropland/sarsfaqs2.php
df <- add_row(df, crop = "greens", kcal = mean(df$kcal[which(df$crop %in% c("broccoli","lettuce","peas","cabbage","celery","cucumbers"))])) # Logic
df <- add_row(df, crop = "other crops", kcal = mean(df$kcal[which(df$kcal<400)])) # Grabs average of all non oil crops
# Raw sugarcane is tough to find data on. One pack of 20 raw sugar cane switzzle sticks is 0.35 lbs or 0.35/20/2.205 kg/stick or 0.0079kg/stick or 8g/stick, 12.5 sticks in 100g, 40 calories per stick, so calories per 100g is 500kcal
df <- add_row(df, crop = "sugarcane", kcal = 500) # https://www.amazon.com/Raw-Sugar-Cane-Swizzle-Sticks/dp/B001CDMN0E, https://www.livestrong.com/article/319462-how-many-calories-are-in-sugar-cane/


# Get still missing crops, this vector should be empty if everything is accounted for
remaining_crops <- setdiff(crop_names, df$crop)
wrong_crops <- setdiff(df$crop, crop_names)

# Create 100g to kg multipler
hundredG_kg <- 10
# Final df contains kcal/100g, so change to kcal/kg
df$kcal_per_kg <- df$kcal*hundredG_kg
df$kcal = NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Get yield to calculate kcal/kg to kcal/m2
#~~~~~~~~~~~~~~~~~~~~~~~ Yield for state

# Read in California yield df from https://quickstats.nass.usda.gov/#0FCCF70B-EA5A-3EA9-AE27-6A0D82A153B4
yield_conus <- crop_yield
yield_state <- yield_conus[which(yield_conus$State=="CALIFORNIA"), ]

# Get yield for state and conus
getYield_func <- function(yield){
  # Change column classes
  yield$Value <- as.numeric(gsub(",","",yield$Value))
  yield <- yield[!is.na(yield$Value), ]
  colnames(yield)[colnames(yield)=="Commodity"] <- "Crop"
  # Get unit column
  yield$Units <- sub(".*IN", "", yield$Data.Item)
  yield$Units <- sub("/.*", "", yield$Units)
  # Group by commodity title
  yield <- yield %>%
    group_by(Crop, Year, Units) %>%
    mutate() %>%
    summarise(Value = mean(Value))
  yield$Crop <- tolower(yield$Crop)
  
  # Convert yields from units/acre to kg/m2
  yield$Value <- ifelse(grepl("BU", yield$Units), yield$Value * BU_kg / acres_m2, 
                 ifelse(grepl("LB", yield$Units), yield$Value * LB_kg / acres_m2, 
                 ifelse(grepl("TONS", yield$Units), yield$Value * TONS_kg / acres_m2, 
                 ifelse(grepl("CWT", yield$Units), yield$Value * CWT_kg / acres_m2, 
                 ifelse(grepl("BOXES", yield$Units), yield$Value * BOXES_kg / acres_m2, 
                 ifelse(grepl("BARRELS", yield$Units), yield$Value * BARRELS_kg / acres_m2, "Not_Converted"))))))
  # kg / m2
  yield$Units = NULL
  names(yield)[names(yield) == "Value"] <- "kg_per_m2"
  yield$kg_per_m2 <- as.numeric(yield$kg_per_m2)
  # Group and average 2018, 2019, 2020, and 2021 yield values 
  yield <- yield %>% group_by(Crop) %>% summarise(kg_per_m2 = mean(kg_per_m2, na.rm = TRUE))
  return(yield)
}

# ~~~~~~~~~~~~~~~~ Combine missing yields at state level from CONUS 

# Run function 
yield_conus <- getYield_func(yield_conus)
yield_state <- getYield_func(yield_state)

# Get yield conus which not in conus tate
yield_conus <- yield_conus[which(!yield_conus$Crop %in% yield_state$Crop), ]

# Merge
yield <- rbind(yield_state, yield_conus)
rm(yield_state, yield_conus)

# ~~~~~~~~~~~~~~~~ Attach to kcal df

# Match yield of crops to kcal df
df$kg_per_m2 = 0
for(i in c(1:nrow(yield))){
  # Get forgone crop type and year from survey data
  crop <- yield[i, ]$Crop
  kg_per_m2 <- yield[i, ]$kg_per_m2 %>% as.numeric()
  # Get matching crop type from dataset
  df[which(grepl(crop, df$crop)), ]$kg_per_m2 <- kg_per_m2
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Make Assumptsion and equivalencies, this may change depending on yield input and annual FoodData

# Make assumptions about USDA reported crop yield types vs USDA NASS CDL croptypes
add_yield <- function(crop_1, crop_2){
  df$kg_per_m2 <- ifelse(df$crop==crop_1, yield$kg_per_m2[which(yield$Crop==crop_2)], df$kg_per_m2)
  return(df) }
# add yields
df <- add_yield("caneberries", "raspberries")
df <- add_yield("herbs", "mint")

# Make assumptions for all crops without data in USDA reported from yield types
add_yield <- function(crop_1, crop_2){
  df$kg_per_m2 <- ifelse(df$crop==crop_1, df$kg_per_m2[which(df$crop==crop_2)], df$kg_per_m2)
  return(df) }
# add yields
df <- add_yield("alfalfa", "winter wheat")
df <- add_yield("other small grains", "winter wheat")
df <- add_yield("switchgrass", "alfalfa")
df <- add_yield("vetch", "alfalfa")
df <- add_yield("cantaloupes", "watermelons")
df <- add_yield("citrus", "oranges")
df <- add_yield("clover/wildflowers", "other hay/non alfalfa")
df <- add_yield("sorghum", "winter wheat")
df <- add_yield("grassland/pasture", "other hay/non alfalfa")
df <- add_yield("peas", "broccoli")
df <- add_yield("triticale", "winter wheat")
df <- add_yield("camelina", "canola")
df <- add_yield("rape seed", "canola")
df <- add_yield("speltz", "oats")
# Couple of custom "mixes" from FSA to CDL Crosswalk
df$kg_per_m2 <- ifelse(df$crop=="other crops", df$kg_per_m2[which(df$kcal_per_kg<4000)], df$kg_per_m2) # non-oil crop kcal propostions
df$kg_per_m2 <- ifelse(df$crop=="other tree crops", mean(df$kg_per_m2[which(df$crop %in% c("almonds","citrus","cherries","nectarines","peaches"))], na.rm=TRUE), df$kg_per_m2)
df$kg_per_m2 <- ifelse(df$crop=="greens", mean(df$kg_per_m2[which(df$crop %in% c("broccoli","lettuce","peas","cabbage","celery","cucumbers"))], na.rm=TRUE), df$kg_per_m2)
df$kg_per_m2 <- ifelse(df$crop=="misc vegs & fruits", mean(df$kg_per_m2[which(df$crop %in% c("broccoli","lettuce","peas","cabbage","celery","cucumbers","citrus","cherries","nectarines","peaches"))], na.rm=TRUE), df$kg_per_m2)
df$kg_per_m2 <- ifelse(df$crop=="pomegranates", (mean(75,150,225,30,400) * BOXES_kg / acres_m2), df$kg_per_m2)  # https://coststudyfiles.ucdavis.edu/uploads/cs_public/d5/bd/d5bdaad2-b874-4b99-a3c2-cc7a89cfc72d/pomegranatevs2010.pdf -=- page 4 table c
df$kg_per_m2 <- ifelse(df$crop=="turnips", (12000 * LB_kg / acres_m2), df$kg_per_m2)# https://www.nrcs.usda.gov/wps/cmis_proxy/https/ecm.nrcs.usda.gov%3A443/fncmis/resources/WEBP/ContentStream/idd_10C6556A-0000-C213-A4DF-1845178DD030/0/ILGM-SPECIES_TurnipsForForage.pdf
df$kg_per_m2 <- ifelse(df$crop=="radishes", (12000 * LB_kg / acres_m2), df$kg_per_m2)# https://www.nrcs.usda.gov/wps/cmis_proxy/https/ecm.nrcs.usda.gov%3A443/fncmis/resources/WEBP/ContentStream/idd_10C6556A-0000-C213-A4DF-1845178DD030/0/ILGM-SPECIES_TurnipsForForage.pdf
df$kg_per_m2 <- ifelse(df$crop=="eggplants", (600 * BU_kg / acres_m2), df$kg_per_m2) # https://extension.missouri.edu/publications/g6369#:~:text=Eggplants%20can%20yield%20500%20to,pounds%20per%20bushel)%20per%20acre.
df$kg_per_m2 <- ifelse(df$crop=="gourds", (5468.75 * LB_kg / acres_m2), df$kg_per_m2) # https://www.uky.edu/ccd/sites/www.uky.edu.ccd/files/gourds.pdf -- page 2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Calcualte kcal/m2 and incorporate Calorie and mass efficiencies, then calculate double crops

# Account for dairy beef Calorie convertion efficiencies and regional dairy to beef ratios
feedsilage_cdl_crops_temp <- feedsilage_cdl_crops %>% tolower()
kcal_converstion_norm = (1/(regional_dairy_beef_ratio+1))*calorie_conversion_efficiency_beef + # Beef correction (gets total % beef and multiplies by converstion efficiency)
                        (regional_dairy_beef_ratio/(regional_dairy_beef_ratio+1))*calorie_conversion_efficiency_dairy # Dairy correction (gets total % dairy and multiplies by converstion efficiency)
df$kcal_per_kg <- ifelse(df$crop %in% feedsilage_cdl_crops_temp, df$kcal_per_kg*kcal_converstion_norm, df$kcal_per_kg)

# Account for mass proportions of seed oil crops (the oil being the primary method of consumption)
seedoil_cdl_crops_temp <- seedoil_cdl_crops %>% tolower()
df$kcal_per_kg <- ifelse(df$crop %in% seedoil_cdl_crops_temp, (df$kcal_per_kg*seedoil_mass_prop*seedoil_feed_prop*kcal_converstion_norm) + # Here we account for cotton seed used in dairy/beef
                        (df$kcal_per_kg*seedoil_mass_prop*(1-seedoil_feed_prop)), df$kcal_per_kg) # and gere we account for cotton seed used in human consumption

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Double crops

# Add double crop sums
dbl_crop_func <- function(crop_1, crop_2){
  # Add double crop row by sum of irrig for both crops for both years of the study
  df <- add_row(df, 
                crop = paste("dbl crop ", crop_1, "/", crop_2, sep = ""), 
                kcal_per_kg = as.numeric(df[which(df$crop==crop_1), ]['kcal_per_kg']) + as.numeric(df[which(df$crop==crop_2), ]['kcal_per_kg']), 
                kg_per_m2 = as.numeric(df[which(df$crop==crop_1), ]['kg_per_m2']) + as.numeric(df[which(df$crop==crop_2), ]['kg_per_m2']))
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
df$crop <- ifelse(df$crop=="dbl crop durum wheat/sorghum", "dbl crop durum wht/sorghum", df$crop)
df$crop <- ifelse(df$crop=="dbl crop lettuce/durum wheat", "dbl crop lettuce/durum wht", df$crop)
df$crop <- ifelse(df$crop=="dbl crop lettuce/cantaloupes", "dbl crop lettuce/cantaloupe", df$crop)
df$crop <- ifelse(df$crop=="dbl crop winter wheat/corn", "dbl crop winwht/corn", df$crop)
df$crop <- ifelse(df$crop=="dbl crop winter wheat/cotton", "dbl crop winwht/cotton", df$crop)
df$crop <- ifelse(df$crop=="dbl crop winter wheat/sorghum", "dbl crop winwht/sorghum", df$crop)
df$crop <- ifelse(df$crop=="dbl crop winter wheat/soybeans", "dbl crop winwht/soybeans", df$crop)

# Calcualte kcal/m2
df$kcal_per_m2 <- df$kcal_per_kg * df$kg_per_m2

# Write final kcal/m2 df
write.csv(df, "Data/Derived/kcal_per_m2.csv", row.names=FALSE)