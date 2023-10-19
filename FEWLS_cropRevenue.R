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
of crop revenue (price recieved) per kg from NASS Agricultual Yield Surveys. Only needs to 
be run if cleanFEWLS() is run, or if surveys are updated. 
"
# --------------------------- #
# Get USDA Survey Crop Prices # 
# --------------------------- #
# Read in California crop cost from https://quickstats.nass.usda.gov
cropdollar <- food_dollarValue

# Change column classes
cropdollar$Value <- as.numeric(gsub(",","",cropdollar$Value))
cropdollar <- cropdollar[!is.na(cropdollar$Value), ]
colnames(cropdollar)[colnames(cropdollar)=="Commodity"] <- "Crop"

# Get unit column
cropdollar$Units <- sub(".*IN", "", cropdollar$Data.Item)
cropdollar$Units <- sub(".*/", "", cropdollar$Units)

# Get only price recieved, on tree for citrus, and non dry basis for beans
cropdollar <- cropdollar[grepl("PRICE RECEIVED", cropdollar$Data.Item), ]
cropdollar <- cropdollar[which(!cropdollar$Units %in% c(" BOX, PHD EQUIV", " BOX, FOB", " TON, DRY BASIS")), ]

# Group by commodity title, year, and units
cropdollar <- cropdollar %>%
  group_by(Crop, Year, Units, .drop=FALSE) %>%
  summarise(Value = mean(Value))
cropdollar$Crop <- tolower(cropdollar$Crop)

# Convert cropdollars from units/acre to kg/m2 -- NOTE the altered arithmatic compared to other instances of 
cropdollar$Value <- ifelse(grepl("BU", cropdollar$Units), cropdollar$Value / BU_kg, 
                    ifelse(grepl("LB", cropdollar$Units), cropdollar$Value / LB_kg, 
                    ifelse(grepl("TON", cropdollar$Units), cropdollar$Value / TONS_kg, 
                    ifelse(grepl("CWT", cropdollar$Units), cropdollar$Value / CWT_kg, 
                    ifelse(grepl("BOXES", cropdollar$Units), cropdollar$Value / BOXES_kg, 
                    ifelse(grepl("BOX", cropdollar$Units), cropdollar$Value / BOX_kg, 
                    ifelse(grepl("BARRELS", cropdollar$Units), cropdollar$Value / BARRELS_kg, "Not_Converted")))))))

# Group by commodity title, year
cropdollar$Value <- as.numeric(cropdollar$Value)
cropdollar <- cropdollar %>%
  group_by(Crop, Year) %>%
  mutate() %>%
  summarise(dollar_kg = mean(Value, na.rm=TRUE))
cropdollar <- cropdollar[which(!is.na(cropdollar$dollar_kg)), ]

# Adjust for inflation -- prices reported as nominal
for(Year in c(unique(cropdollar$Year))){
  infl_rate_mult <- inflation_ratePPI$infl_rate_mult[which(inflation_ratePPI$Year==Year)]
  cropdollar$dollar_kg <- ifelse(cropdollar$Year==Year, cropdollar$dollar_kg * infl_rate_mult, cropdollar$dollar_kg)
}

# Get continous vector of crop prices for missing year (interpolate or extrapolate)
cropdollar$Crop <- cropdollar$Crop %>% as.character()
df <- cropdollar[which(cropdollar$Crop=="Non Existent"),] %>% as.data.frame()
unique_crops = unique(cropdollar$Crop) 
for(crop in unique_crops){
  # Subset for every crop
  cropdollar_temp <- cropdollar[which(cropdollar$Crop==crop), ] %>% as.data.frame()
  cropdollar_temp$Crop = NULL
  # Get continous year vectors -- Interpolate and extrapolate for start year to and year for years with missing data
  cropdollar_temp <- getContinous_installPeriodVector(cropdollar_temp)
  # If there remains any missing years, calculate the mean inflation adjusted value across all years
  cropdollar_temp$dollar_kg <- ifelse(is.na(cropdollar_temp$dollar_kg), mean(cropdollar_temp$dollar_kg, na.rm=TRUE), cropdollar_temp$dollar_kg)
  # Return crop name and save to df
  cropdollar_temp$Crop = crop
  df <- rbind(df, cropdollar_temp)
}
# Get results
cropdollar <- df

# Export df
write.csv(cropdollar, "Data/Derived/crop_revenue_kg.csv", row.names=FALSE)