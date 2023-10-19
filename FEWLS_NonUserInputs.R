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
This is an internal file (hence, the non-user input variables name), that is somewhat 
redundent, but helped solve issues with R saving environment variables between separate
R scripts. Within each function, these variables are also called, but to ensure they are
present, this function is also called within *FEWLS_REsourceModelResults.R* and
*FEWLS_EconModelResults.R*.
"

#////////////////////////#
#________________________#
#    SETUP vARIABLES     #
#________________________#
#\\\\\\\\\\\\\\\\\\\\\\\\#

# Non user variables 
start_install_year = in_df$Year %>% range() %>% min() %>% as.numeric() # first year of installation in dataset 
end_install_year  =  in_df$Year %>% range() %>% max() %>% as.numeric() # last year of installation in dataset (same as start year if only one array)

# SETUP Generated Variables -- no input necessary
install_period = end_install_year - start_install_year
model_length = install_period + system_lifespan 

# Set Variables for yield through time
food_yield_base = c( (1-yield_deficit_base) + yield_time_base * c(1:(model_length)))
food_yield_best = c( (1-yield_deficit_best) + yield_time_best * c(1:(model_length)))
food_yield_worst= c( (1-yield_deficit_worst) + yield_time_worst * c(1:(model_length)))