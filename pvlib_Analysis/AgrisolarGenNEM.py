# -*- coding: utf-8 -*-
"""
######################################################################################
#  ____        _              _____ _______        ___     ____    _____           _ #
# / ___|  ___ | | __ _ _ __  |  ___| ____\ \      / / |   / ___|  |_   _|__   ___ | |#
# \___ \ / _ \| |/ _` | '__| | |_  |  _|  \ \ /\ / /| |   \___ \    | |/ _ \ / _ \| |#
#  ___) | (_) | | (_| | |    |  _| | |___  \ V  V / | |___ ___) |   | | (_) | (_) | |#
# |____/ \___/|_|\__,_|_|    |_|   |_____|  \_/\_/  |_____|____/    |_|\___/ \___/|_|#
#                                                                                    #
######################################################################################

-- De"script"tion -- 
This is the execute file for the FEWLS tool sub-process of modeling electricity generation 
and net energy metering (NEM) returns. This code and data are self contained, and while 
not directly required by FEWLS to run, annual generation and returns ($) are in the 
example dataset as "Gen_20##" and "Econ_20##" respectively, and are required for 
FEWLS_ResourceModelResults.R and FEWLS_EconModelResults.R. For the default data, this 
file can be run as is. Make sure to specify the working directory (in_dir) that contains
the FEWLS_tool folder (should end with ...FEWLS_tool\pvlib_Analysis).

Created on Sun Feb 28 14:21:20 2021
Last edited on Sunday Oct 1 16:52:10 2023

@author: Siddharth Shukla (shuklas7@msu.edu) 
Contact: Siddharth Shukla (shuklas7@msu.edu) or Jacob Stid (stidjaco@msu.edu)
"""

import os
import glob
import pvlib
import pandas as pd
import numpy as np
import time 
#Remove all warnings from the kernel
import warnings
warnings.filterwarnings("ignore")

#start counting the time now 
t0=time.time()
from scipy.spatial import distance

# Designate folder paths 
#in_dir = r'S:\Users\stidjaco\R_files\FEWLS_tool\pvlib_Analysis'
out_dir = '\FEWLS_pvlibResults'

# Set year (for testing and debugging only)
# yr_inst = 2018

"################################################################################################################################################### -- Automate across install years"

# Funciton to run model for every installation year -- runs whole document
def getAgrisolarGenNEM(yr_inst):
    
    # Run for arrays with installation year of interest
    beginning_year = yr_inst # 2018 # installation year
    
    ""##############################################""
    
    #User defined variables related to PV
    gamma_s=180  #Azimuth for PV (usually 180 , please see pvlib documentation to change the value)
    albdo=0.2   #Albedo i.e. percentage reflection from the ground
    max_power_temp_coeff=0.5    #This is the value of maximum temperature power coefficient   
    last_year=2019 # Does not change, this is equivalent to Jan 1st, 2019, so model runs through Dec 31st, 2018
    capacity_threshold = 1 # MW, according to NEM guidelines
    #New degradation rates based on median rates Jordan et al., 2022, although, not statistical significant difference, so all degradation averages -0.75$/yr from slightly more conservative and updated (https://doi.org/10.1002/pip.3566)
    #used to be 0.6%/year for pre-2010 installations and 0.3%/year for post-2010 installations (Jordan et al. 2016). Jordan et al., 2022 also reports delineation by mount tech from -0.68%/yr for fixed and -0.76%/yr for single, withough sigdiff
    annual_deg = 0.0075 # -0.75%/yr: Efficiencies must be in decimal float, not %
    
    
    #FUNCTION TO GET NBC: Non-Bypassable charges are $0.02 to $0.03/kWh, we use $0.03/kWh to provide a conservative estimate
    #Important note: Keep the nbc as zero when the doing the analysis for modules from 2008-2016
    def getNBC(Year):
        if Year > 2016:
            nbc = 0.03  # $0.03/kWh -- This the non-bypassable charge the net metering availing cusotomers under NEM 2.0 - Jan 1, 2017
        else: 
            nbc = 0.00  # $0.00/kWh -- All installations prior to 2017 grandfathered into NEM 1.0 
        return nbc
    
    
    #FUNCTION TO GET EFFICIENCY: LBNL Tracking the Sun Report -- Really only reported efficiency is need, but others are included for future work
    eff_df_loc = '\efficiency_LBNL_08to18_monopoly.csv'
    eff_df = pd.read_csv(in_dir + eff_df_loc).to_numpy() #This factor accounts for all losses
    
    # Reported Efficiency
    def getReportedEfficiency(Year):
        median_pv_efficiency = eff_df[np.where(eff_df[:,0] == Year), 3].item(0) # Note, first column (index 0) must be Year, and column index 3 must be reported eff: Efficiencies must be in decimal float, not %
        return median_pv_efficiency
    
    
    # Monocrystalline efficiency
    def getMonoEfficiency(Year):
        mono_pv_efficiency = eff_df[np.where(eff_df[:,0] == Year), 1].item(0) # Note, first column (index 0) must be Year, and column index 1 must be monocrystalline eff: Efficiencies must be in decimal float, not %
        return mono_pv_efficiency
    
    
    # Polycrystalline efficiency
    def getMultiEfficiency(Year):
        multi_pv_efficiency = eff_df[np.where(eff_df[:,0] == Year), 2].item(0) # Note, first column (index 0) must be Year, and column index 2 must be multicrystalline eff: Efficiencies must be in decimal float, not %
        return multi_pv_efficiency
    
    
    #FUNCTION TO GET LOSS: From NREL Cost Benchmark, acquire temporally relevant pre-inverter, inverter efficiency, and soiling loss sum (%)
    derate_file='\Derating_Loss_Total_08to18.csv'
    derating_factorDF= pd.read_csv(in_dir + derate_file).to_numpy() #This factor accounts for all losses
    #Define function to get loss based on installation year
    def getDerate_factorLoss(Year):
        derating_factor=derating_factorDF[np.where(derating_factorDF[:,0] == Year), 1].item(0)  # Note, first column (index 0) must be Year, column index 1 must be derate: Losses in dataframe are in %, so are divided by 100 in next line
        loss=(100-derating_factor)/100 # End result in decimal, and is the proportion of generation NOT loss (ie. sum of loss of 15%, this variable will be 0.85)
        return loss
    
    
    # Depending on install year of interest, run internal variable functions
    loss = getDerate_factorLoss(beginning_year)
    nbc = getNBC(beginning_year)
    median_pv_efficiency = getReportedEfficiency(beginning_year)
    mono_pv_efficiency = getMonoEfficiency(beginning_year)
    multi_pv_efficiency = getMultiEfficiency(beginning_year)
    
    
    #create a numpy array which contains the yearly production and the type of tracker panel for each of the dataset
    matrix_size = 1000 # creates matrix of this many rows. Must be larger than size of dataset (ie number of arrays)   
    median_final_matrix=np.zeros((matrix_size,(last_year-beginning_year)+2))
    multi_final_matrix=np.zeros((matrix_size,(last_year-beginning_year)+2))
    cost_final_matrix=np.zeros((matrix_size,(last_year-beginning_year)+2))
    econ_final_matrix=np.zeros((matrix_size,(last_year-beginning_year)+2)) #July 19
    
    
    #time zone of the central valley (Adjusted according to the command (sign reversed)
    t_z='Etc/GMT+8'
    
    #Read the PV id datset csv file
    #pv_id_file_location='G:\Folder_for_PV_id_algo'
    
    #pv_id_file_location='D:\Dropbox\My Research stuff\INFEWS project related stuff\Data related to Jake project\Folder_for_PV_id_algo'
    pv_id_file='\SolarFEWE_generation_input.csv'
    pv_id_dataset_temp=pd.read_csv(in_dir + pv_id_file, index_col=False)
    
    # TEMP: Select array of an index for special case debugging
    #indx = 2
    #pv_id_dataset_temp=pv_id_dataset_temp.loc[pv_id_dataset_temp['Index'] == indx, :].reset_index(drop=True)
    
    # Subset by capacity threshold
    pv_id_dataset=pv_id_dataset_temp.loc[pv_id_dataset_temp['Capacity'] < capacity_threshold, :].reset_index(drop=True) # Returns only commercial-scale arrays (receiving NEM returns) or less than 1 MW installations, reset index solves subset issue return out of bounds index
    pv_id_workable=pv_id_dataset.to_numpy()
    
    #Defining a function to calculate miniumum euclidean distance between two set of points
    def closest_node(node, nodes):
        closest_index = distance.cdist([node], nodes).argmin()
        return nodes[closest_index]
    
    # Make a list of row indices 
    #pv_instYear_list = list(np.asarray(pv_id_workable[:,4] == beginning_year).nonzero()) # Creates numpy list of indices from workable 
    pv_instYear_list = pv_id_dataset.index[pv_id_dataset['Year'] == beginning_year].tolist() # Creates generic list of indices from pandas df
    
    # Get variables for progress reporting (number of arrays within installation year)
    numArrays = len(pv_instYear_list)
    
    #Running this loop the same number of times same as the number of points to be analyzed
    for k in pv_instYear_list: #range(43,45):
        
        # Get variables for progress reporting (array in subset)
        k_indx = pv_instYear_list.index(k)
        
        #Read the installation year , latitude , longitude, area and tilt of the PV ; Also read from the file weather you are dealing with fixed axis or tracker PV
        pv_inst_year=pv_id_workable[k,5]
        pv_lat=pv_id_workable[k,18]
        pv_long=pv_id_workable[k,19]
        pv_area=pv_id_workable[k,3]
        tilt=pv_id_workable[k,25]
        pv_type=pv_id_workable[k,1]
        utility_name=pv_id_workable[k,17]
        pv_index=pv_id_workable[k,0]


        #Tying the latitude and longitude in form of a dictionary so that closest node fucntion can work properly
        pv_coordinates=(pv_lat,pv_long)
        
        
        current_year=beginning_year #Intialize this variable with the year when we first started the analysis
        #Initialize these variables with efficiencies of the current year ; this number would be degraded every year 
        median_eff_current_year=median_pv_efficiency 
        multi_eff_current_year=multi_pv_efficiency
        
        #The below for loop does yearwise calculations i.e. 2017, 2018etc.
        for l in range(last_year-beginning_year):
            
            
            #Read the names of all the csv weather files available for that year
            weather_file_path=in_dir + '\Weather_files'
            weather_file_folder="\\"+str(current_year)
            weather_file_location=weather_file_path+weather_file_folder
            extension='csv'
            os.chdir(weather_file_location)
            result=glob.glob('*.{}'.format(extension))
            
            
            #Reading the latitude value from the csv file names and storing it in a list called lat_list (This would need to changed acccording to nomenclature of the new weather files)
            lat_list=[string[7:12] if string[0]=='1' else  string[6:11]  for string in result]
            
            #Reading the longitude value from the csv file names and storing it in a list called long_list
            long_list=[string[14:20] if string[0]=='1' else  string[13:19]  for string in result]
            
            #Converting the lat_list and long list values to floating point numbers for calculation
            lat_list=[float(i) for i in lat_list]
            long_list=[float(j) for j in long_list]
            
            #Multiplying longitude values by -1
            long_list=[neg*(-1) for neg in long_list]
            
            
            #Combining the two lists together to form a calculable 2-D array
            coordinate_list=np.zeros((len(lat_list),2))
            
            #Filling coordinae list with all the latitude and longitudes for which the weather files are available 
            for cod in range(len(lat_list)):
                coordinate_list[cod,0]=lat_list[cod]
                coordinate_list[cod,1]=long_list[cod]
            
            
            #Call the closest node function and find the coordinates available for the closest weather file
            chosen_weather_coord=closest_node(pv_coordinates,coordinate_list)
            
            #Convert the caclculated latitude and longitude values to string data type
            chosen_weather_coord=chosen_weather_coord.astype(str)
            
            #Find the weather file with the corresponding latitude and longitude values 
            
            for file in os.listdir(weather_file_location):
                if (chosen_weather_coord[0] in file and chosen_weather_coord[1] in file):
                    weather_file='\\'+file
           
            
            """========Reading the  DNI, GHI, DHI , temperature and elevation values from the weather files=========="""
            #Generating the time stamp excel to be used in thepvlib module
            tim_stmp=pd.date_range('2018-01-01',periods=8760,freq='H')
            tim_stmp=tim_stmp.tz_localize(t_z)
            
            #Caling in the values of direct normal irradiance from the selected weather file
            dni=pd.read_csv(weather_file_location+weather_file,usecols=[6],skiprows=2)
            #dni_vector=dni.as_matrix()
            dni_vector=dni.values
            
            #Calling in the values of direct horizontal irradiance
            dhi=pd.read_csv(weather_file_location + weather_file,usecols=[5],skiprows=2,)
            #dhi_vector=dhi.as_matrix()
            dhi_vector=dhi.values
            
            #Calling in the values of GHI
            ghi=pd.read_csv(weather_file_location + weather_file,usecols=[7],skiprows=2,)
            #ghi_vector=ghi.as_matrix()
            ghi_vector=ghi.values
            
            #Defining an extra DNI vector for the sake of model to run
            extra_dni_vector=np.zeros((8760,1))
            
            
            
            #Finding out the elevation value from the weather file 
            ele=pd.read_csv(weather_file_location + weather_file,usecols=[8],error_bad_lines=False,warn_bad_lines=False)
            elev=float(ele.values[0])
            
            #Finding out the temperature value from the weather file
            hourly_temperature=pd.read_csv(weather_file_location + weather_file,usecols=[9],skiprows=2,)
            hourly_temperature_vector=hourly_temperature.values
            """======================================================================================================================""" 
            
           
            """===========Calling pvlib modules to calaculate total irradaince on each type of PV modules================================"""
           
            #Empty list to store total irradiance from the pvlib modules
            tot_irrad=np.zeros((8760,1))
            
            #Create aother code to call PV generation for either fixed or single axis tracker systems
            #Finding out the hourly irradinace for a year with a Fixed PV system
            if(pv_type=='fixed_axis'):
                #This code calculates the hourly irradiance with a fixed PV system  
                #Finding out the hourly irradinace for a year with a Fixed PV system
                for i in range(8760):
                
                
                
                    my_azm=pvlib.solarposition.pyephem(tim_stmp[[i]],latitude=pv_lat,longitude=pv_long,altitude=elev,pressure=101325,temperature=hourly_temperature_vector[i],horizon='+0:00')
                
                
                    irrad_comp=pvlib.irradiance.get_total_irradiance(surface_tilt=tilt,surface_azimuth=gamma_s,solar_zenith=my_azm.zenith,solar_azimuth=my_azm.azimuth,dni=dni_vector[i],ghi=ghi_vector[i],dhi=dhi_vector[i],dni_extra=extra_dni_vector[i],airmass=None,albedo=albdo,surface_type=None,model='isotropic')                         
                    tot_irrad[i]=(loss*(irrad_comp.poa_global))/1000
                
                    pv_code=1  #Flag it as pv code 1 if it is a fixed PV system 
                
            else:
                
                #This code calculates the hourly irradiance with a single axis tracking system
                
                for i in range(8760):
                    
                    my_azm=pvlib.solarposition.pyephem(tim_stmp[[i]],latitude=pv_lat,longitude=pv_long,altitude=elev,pressure=101325,temperature=hourly_temperature_vector[i],horizon='+0:00')
                    
                    
                    tracker_func=pvlib.tracking.singleaxis(apparent_zenith=my_azm.apparent_zenith,apparent_azimuth=my_azm.apparent_azimuth,axis_tilt=tilt,axis_azimuth=180,max_angle=45,backtrack=True,)
                    
                    
                    #You will need to find a new method to find the total irradaince as the single axis tracker method has been deprectaed in v.10.1
                    #Calling instance of Single Axis tracking class
                    
                    irrad_comp=pvlib.irradiance.get_total_irradiance(surface_tilt=tracker_func.surface_tilt, surface_azimuth=tracker_func.surface_azimuth, solar_zenith=my_azm.zenith, solar_azimuth=my_azm.azimuth, dni=dni_vector[i], ghi=ghi_vector[i], dhi=dhi_vector[i])
                    
    # =============================================================================
    #                 SAT=pvlib.tracking.SingleAxisTracker()
    #                 
    #                 
    # 
    #                 #Calculating irradaince at each hour from a single axis tracking system
    #                 irrad_comp=SAT.get_irradiance(surface_tilt=tracker_func.surface_tilt,surface_azimuth=tracker_func.surface_azimuth,solar_zenith=my_azm.zenith,solar_azimuth=my_azm.azimuth,dni=dni_vector[i],ghi=ghi_vector[i],dhi=dhi_vector[i],dni_extra=extra_dni_vector[i],airmass=None,model='isotropic')
    # =============================================================================
                     
                    
                    
                    
                    tot_irrad[i]=(loss*(irrad_comp.poa_global))/1000
                    
                    
                
                    pv_code=0  #Flag it as PV code 0 if it is a single axis tracker 
            
            """========================================================================================================================================================================================="""
            
            
            #Reading windspeed for dataframe 
            wind_speed_dataframe=pd.read_csv(weather_file_location + weather_file,usecols=[8],skiprows=2,)
            #This vector contains hourly wind speed in the 
            wind_speed=wind_speed_dataframe.values
            
            #Rearranging all the vector in 365 by  24 format for easire calculation
            
            tot_irrad=np.reshape(tot_irrad,(24,365),order='F')
            tot_irrad=np.nan_to_num(tot_irrad)  #Convert nan value to zero if any 
            
            wind_speed=np.reshape(wind_speed,(24,365),order='F')
            temperature=np.reshape(hourly_temperature_vector,(24,365),order='F')
            
            term1=(np.exp(-0.075*wind_speed-3.56))*tot_irrad
            term2=temperature
            
            module_back_temp= term1+term2
            
            cell_temp=module_back_temp + (tot_irrad)*3
            
            
            #Check if the temperature is greater than 25 degrees or not and apply the condition for PV generation equation
            
            #this matrix wiil have  value of 1 if the temperature is less than 25 degree celsius and less than 1 whenever temperature goes more than 25
            
            matrix_a=np.where(cell_temp<25,1,cell_temp)
            
            temp_coeff_matrix=np.where(matrix_a>1,(1-((matrix_a-25)*(1-max_power_temp_coeff)/100)),matrix_a)
                
            #Calculating PV generation for this year for both median and multi si solar cells         
            median_pv_power_matrix=median_eff_current_year*pv_area*np.array(tot_irrad)*temp_coeff_matrix
            multi_pv_power_matrix=multi_eff_current_year*pv_area*np.array(tot_irrad)*temp_coeff_matrix
            
            
            #Code to find out the economic value of the electricity 
            if(utility_name=='Power and Water Resource Pooling Authority' or utility_name=='Pacific Gas & Electric Company' or utility_name=='Eastside Power Authority' or utility_name=='Lodi Electric Utility'):
                utility_key='pge'
            elif(utility_name=='Merced Irrigation District'):
                utility_key='mer'
            elif(utility_name=='Sacramento Municipal Utility District'):
                utility_key='smu'
            elif(utility_name=='Southern California Edison'):
                utility_key='sce'
            elif(utility_name=='Turlock Irrigation District'):
                utility_key='tur'
            elif(utility_name=='Modesto Irrigation District'):
                utility_key='mod'
            
            
            #I can add the folder name to the path rather than the csv file name itself
            price_file_path= in_dir + '\\utility_rates' 
            price_file_folder="\\"+str(current_year)
            price_file_utility_folder="\\"+utility_key
            
            
            
            # write a code here to choose the file name based on the 
            #energy_price_file_name=
            
    # =============================================================================
    #         price_file_path='D:\Dropbox\My Research stuff\INFEWS project related stuff\Data related to Jake project\Folder_for_PV_id_algo\\utility_rates - older rates'
    #         price_file_folder="\\"+str(current_year)
    #         price_file_name="\\"+'9999_'+utility_key+'_'+str(current_year)+'.xlsx'
    # =============================================================================
                    
            #price_file_use=pd.read_excel(price_file_path+price_file_folder+price_file_name, sheet_name=0,header=None) #reading the excel file for the concerned year and the utility
            
            #Rounding the area of panel upto three places of decimal for exact match in the econ matrix.
            pv_area=round(float(pv_area),3)
            
            #Read the econ dataframe with monthly energy demand for each panel and find the entry with the above panel area
            monthly_energy_file=in_dir + "\SolarFEWE_econ_inputs.csv"
            
            monthly_energy_dataframe=pd.read_csv(monthly_energy_file,index_col=False)
            
            #Rounding up the panel area upto three decimal palces for exact match with the earlier listed pv_area
            monthly_energy_dataframe['pnl_a']=monthly_energy_dataframe['pnl_a'].round(3)
            
            #Finding the row number in the econ dataframe with the selected panel area 
            row_number_load=monthly_energy_dataframe[monthly_energy_dataframe['pnl_a']==pv_area].index[0]
            
            
            #Pick out the maximum load for that farm (occurs for all in July) with the units of MWh/month
            max_monthly_load=monthly_energy_dataframe.iloc[row_number_load]['Jul']   
            
            
            """+++++++++++++++++++++Creating an hourly load dataframe based on the monthly data given by FEWLS Methods +++++++++++++++++++++++++++++++++++"""
            
            #Slicing off the monthly load dataframe and converting to a 1-D numpy array 
            select_monthly_load=monthly_energy_dataframe.iloc[row_number_load, 2:13] 
                  
            monthly_load=select_monthly_load.to_list()
            
            
            # Define the number of days in each month (assuming non-leap year)
            days_in_months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            
            
            # Create an empty NumPy array with the desired shape
            hourly_load = np.zeros((24, 365))
                   
            # Populate the array with the hourly load values
            
            day_index=0
                    
            #running a for loop that converts the yearly into hourly values        
            #The for loop below makes a tuple of the equally sized lists where "month_load" variable  goes through the "monthly_load" list and "days_in_month" variable goes through the list "days_in_months"
            for month_load, days_in_month in zip(monthly_load, days_in_months):
                
                # Calculate the daily hourly load in kWh
                per_hour_load=(month_load*1000) /(24* days_in_month)
                
                for day in range(days_in_month):
                    for hour in range(24):
                        
                        hourly_load[hour,day_index]=per_hour_load
                        
                    day_index +=1
                    
            """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++This part if for PGE utility++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
            
            #Here we would need to do multiple checks and selection based on the utility name, small vs. large field , utility rate structures etc.           
            if (utility_key=='pge'):
                
                #Converting the monthly MWh/month demand in kWh/h demand for easy check against utility large vs. small land classification
                hourly_demand=(max_monthly_load*1000)/(31*24)
                
                #Assign a field size tag based on the fact that whether this field had demand of more or less than 35kW
                if hourly_demand<35:
                    
                    field_size_tag='SmAg'
                    
                else:
                    
                    field_size_tag='LgAg'
                    
                    
                #Reading and lsiting all the rate files in the utility specfic folder for that year
                rate_files=os.listdir(price_file_path+price_file_folder+price_file_utility_folder)
                
                #Parse out the energy files based on the field size tag
                energy_file=[file for file in rate_files if field_size_tag in file and "energy" in file]
                
                #Parse out the demand files based on the field size tag
                demand_file=[file for file in rate_files if field_size_tag in file and "demand" in file]
                
                #Parse out the customer charge files based on the field size tag
                customer_charge_file=[file for file in rate_files if field_size_tag in file and "customer" in file]
                
                #Calling in the energy charge values
                energy_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"\\"+energy_file[0],header=None)            
                energy_charge_matrix=energy_charge_matrix.values
                
    # =============================================================================
    #             #Calling in the demand charge values 
    #             demand_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"\\"+demand_file[0],header=None)           
    #             demand_charge_matrix=demand_charge_matrix.values
    #             
    #             #Calling in the customer charge values            
    #             customer_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"\\"+customer_charge_file[0], header=None)
    #             customer_charge_matrix=customer_charge_matrix.values
    # =============================================================================
                           
    # =============================================================================
    #             #Converting the demand charges and customer charges to a per hour format for ease of combining the energy charges           
    #             counter=0            
    #             for day in days_in_months:                
    #                 for i in range(day):                    
    #                     demand_charge_matrix[:,counter]=demand_charge_matrix[:,counter]/(day*24)                    
    #                     customer_charge_matrix[:,counter]=customer_charge_matrix[:,counter]/(day*24)                    
    #                     counter+=1
    #             
    #             #Calculating the demand charges and customer charge  in the per hour format (assuming that load demand is uniform throughout the month)
    #             demand_charges=np.multiply(demand_charge_matrix,hourly_load)            
    #             fixed_monthly_charges=customer_charge_matrix
    # =============================================================================
            """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""" 
            
            """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++This part if for Merced irrigation district utility++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
            
            if (utility_key=='mer'):
                
                #No large or small utility demarcation required here 
                
                rate_files=os.listdir(price_file_path+price_file_folder+price_file_utility_folder)
                
                #Parse out the energy files
                energy_file=[file for file in rate_files if "energy_tier" in file]
                
                #parse out energy adjustment files
                energy_adj_file=[file for file in rate_files if "energy_adj" in file]
                            
                #parse out the demand charge file 
                demand_file=[file for file in rate_files if "demand" in file]
                
                #parse out the fixed charge file
                fixed_charge_file=[file for file in rate_files if "fixed" in file ]
                
                #Calling in the energy charge values
                energy_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+energy_file[0],header=None)
                energy_charge_matrix=energy_charge_matrix.values
                
                #Calling in the energy charge adj values
                energy_adj_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+energy_adj_file[0],header=None)
                energy_adj_charge_matrix=energy_adj_charge_matrix.values
                
    # =============================================================================
    #             #Calling in the demand charge values
    #             demand_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+demand_file[0], header=None)
    #             demand_charge_matrix=demand_charge_matrix.values
    #             
    #             #Calling in the fixed charge values
    #             fixed_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+fixed_charge_file[0], header=None)
    #             fixed_charge_matrix=fixed_charge_matrix.values
    # =============================================================================
                
                #Combining the energy and energy adjustment charges
                energy_charge_matrix=energy_charge_matrix+energy_adj_charge_matrix
                
                
    # =============================================================================
    #             #Converting the demand charges and customer charges to a per hour format for ease of combining the energy charges
    #             counter=0
    #             for day in days_in_months:
    #                 for i in range(day):
    #                     demand_charge_matrix[:,counter]=demand_charge_matrix[:,counter]/(day*24)
    #                     fixed_charge_matrix[:,counter]=fixed_charge_matrix[:,counter]/(day*24)
    #                     counter+=1
    #                     
    #             #Calculating the demand charges and customer charge  in the per hour format (assuming that load demand is uniform throughout the month)
    #             demand_charges=np.multiply(demand_charge_matrix,hourly_load)
    #             fixed_monthly_charges=fixed_charge_matrix
    #                     
    # =============================================================================
    # =============================================================================
    #             
    #         """+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++Calculate the money given to the utility+++++++++++++++++++++++++++++++++++++++"""
    #         
    #         #Find the charges to the customer based on hourly load and generation
    #         load_surplus_matrix=hourly_load-mono_pv_power_matrix
    #         
    #         #Put all the negative numbers in the above matrix as zero
    #         load_surplus_nonzero=np.where(load_surplus_matrix<0,0,load_surplus_matrix)
    #         
    #         #Calculate the energy charges that customer had to pay to the utility 
    #         energy_charges_lost=np.multiply(load_surplus_nonzero,energy_charge_matrix)
    #                     
    #         #This is the amount paid to utility
    #         money_to_utility=energy_charges_lost+demand_charges+fixed_monthly_charges
    #         
    #         """"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    #         
    # =============================================================================
                
            if (utility_key=='smu'):
                
                #Converting the monthly MWh/month demand in kWh/h demand for easy check against utility large vs. small land classification
                hourly_demand=(max_monthly_load*1000)/(31*24)
                
                #Assign a field size tag based on the fact that whether this field had demand of more or less than 30kW
                if hourly_demand<30:               
                    field_size_tag='SmAg'               
                else:                
                    field_size_tag='LgAg'
                    
                    
                #Reading and lsiting all the rate files in the utility specfic folder for that year
                rate_files=os.listdir(price_file_path+price_file_folder+price_file_utility_folder)
                
                #Parse out the energy files based on the field size tag
                energy_file=[file for file in rate_files if field_size_tag in file and "energy_20" in file]
                
                #Parse out energy adjustment files
                energy_adj_file=[file for file in rate_files if field_size_tag in file and "energy_adj" in file]
                
                #Parse out the energy addition files
                energy_add_file=[file for file in rate_files if field_size_tag in file and "energy_add" in file]
                
                
                #Calling in the energy charge values
                energy_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+energy_file[0],header=None)
                energy_charge_matrix=energy_charge_matrix.values
                
                #Calling in energy charge adj values
                energy_adj_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+energy_adj_file[0],header=None)
                energy_adj_charge_matrix=energy_adj_charge_matrix.values
                
                #Calling in energy charge add values
                energy_add_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+energy_add_file[0], header=None)
                energy_add_charge_matrix=energy_add_charge_matrix.values
                
                #Combining energy, energy adjustment and energy addition charges
                energy_charge_matrix=energy_charge_matrix+energy_adj_charge_matrix+energy_add_charge_matrix
                
                
            if (utility_key=='sce'):
                
                
                #Converting the monthly MWh/month demand in kWh/h demand for easy check against utility large vs. small land classification
                hourly_demand=(max_monthly_load*1000)/(31*24)
                
                #Assign a field size tag based on the fact that whether this field has demand of more or less than 200 kW
                if hourly_demand<200:
                    field_size_tag='SmAg'
                else:
                    field_size_tag='LgAg'
                    
                #Reading and lsiting all the rate files in the utility specfic folder for that year
                rate_files=os.listdir(price_file_path+price_file_folder+price_file_utility_folder)
                
                #Parse out the first part of energy charges files based on the field size tag
                energy_file_1=[file for file in rate_files if field_size_tag in file and "genUG" in file]
                
                #Parse out the second part of the energy charges files based on the field size tag
                energy_file_2=[file for file in rate_files if field_size_tag in file and "deliv" in file]
                
                #Calling in energy charges 1 values
                energy_charge_matrix_1=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"\\"+energy_file_1[0], header=None)
                energy_charge_matrix_1=energy_charge_matrix_1.values
                
                #Calling in energy charges 2 values 
                energy_charge_matrix_2=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"\\"+ energy_file_2[0], header=None)
                energy_charge_matrix_2=energy_charge_matrix_2.values
                
                
                energy_charge_matrix=energy_charge_matrix_1+energy_charge_matrix_2
                
                
            if (utility_key=='tur'):
                          
                #No large or small utility demarcation required here 
                
                rate_files=os.listdir(price_file_path+price_file_folder+price_file_utility_folder)
                
                #Parse out the energy files 
                energy_file=[file for file in rate_files if "energy_tier" in file ]
                
                #Parse out the energy adjustment files 
                energy_adj_file=[file for file in rate_files if "energy_adj" in file]
                
                #Calling in the energy charge values
                energy_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+ energy_file[0], header=None)
                energy_charge_matrix=energy_charge_matrix.values
                
                #Calling in the energy charge adj values 
                energy_adj_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+energy_adj_file[0], header=None)
                energy_adj_charge_matrix=energy_adj_charge_matrix.values
                
                #Calling the energy and energy adjustment charges
                energy_charge_matrix=energy_charge_matrix+energy_charge_matrix
                
                           
            if (utility_key=='mod'):
                
                #print("The Modesto part of the code is working perfectlys")
                
                #Converting the monthly MWh/month demand in kWh/h demand for easy check against utility large vs. small land classification for the demand
                hourly_demand=(max_monthly_load*1000)/(31*24)
                
                
                #Check if the maximum monthly load crosses 5000 kWh (or 5MWh) and the demand crosses 10hp (7.46 kW); assign the field size tag based on that
                if (max_monthly_load<5 and hourly_demand<7.46):
                    field_size_tag='tier1'
                else:
                    field_size_tag='tier2'
                    
                #Reading and lsiting all the rate files in the utility specfic folder for that year
                rate_files=os.listdir(price_file_path+price_file_folder+price_file_utility_folder)
                
                #Parse out the energy files based on the field size tag
                energy_file=[file for file in rate_files if field_size_tag in file and "energy_tier" in file]
                
                #Parse out the energy adjustment files based on the field size tag
                energy_adj_file=[file for file in rate_files if field_size_tag in file and "energy_adj" in file]
                
                #Calling in the energy charge values
                energy_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+energy_file[0], header=None)
                energy_charge_matrix=energy_charge_matrix.values
                
                
                #Calling in the energy charge adj values
                energy_adj_charge_matrix=pd.read_csv(price_file_path+price_file_folder+price_file_utility_folder+"//"+energy_adj_file[0], header=None)
                energy_adj_charge_matrix=energy_adj_charge_matrix.values
                
                
                energy_charge_matrix=energy_charge_matrix+energy_adj_charge_matrix
                
                
            """+++++++++++++++++++++++++++++++++++++++++++++++Calculate the money received from the utility ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
            
            #Filter out the locations where there was no PV generation whatsoever
            zero_pv_locations=median_pv_power_matrix==0
            
            #Making copy of the hourly load matrix so I maintain the original record 
            hourly_load_copy=np.copy(hourly_load)
            
            #Choose the load from those hours only where there was a non-zero PV generation
            hourly_load_copy[zero_pv_locations]=0
            
            #Find out the lower of the load and pv generation value as this value would be credited with energy, demand and customer charges 
            load_satisified_by_pv=np.minimum(median_pv_power_matrix,hourly_load_copy)
            
            #Putting customer charges for the times when the PV generation was present as zero
            
            #Find the money saved by utilizing the pv generation for satisfying the load; 
            money_saved=np.multiply(load_satisified_by_pv,energy_charge_matrix)
            
            
            #Find the hours where the PV generation was more than the load 
            pv_surplus_interim=median_pv_power_matrix-hourly_load
            
            
            pv_surplus_matrix=np.where(pv_surplus_interim<0,0,pv_surplus_interim)
                        
            
            energy_charges_minus_nbc=energy_charge_matrix-nbc
            
            #Find the economic gains due to the electricity sold back to the grid
            money_earned=np.multiply(pv_surplus_matrix,energy_charges_minus_nbc)
            
            
            money_from_utility= (money_saved + money_earned)
            
            """+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
            
            
            #elec_price=energy_charge_matrix  #Convert the called electricity price matrix in the readable form
                   
            #econ_value_year=np.multiply(mono_pv_power_matrix,elec_price) 
            
            econ_value_year=money_from_utility
            
            econ_final_matrix[k,0]=pv_code
            econ_final_matrix[k,1]=pv_index # pv_area
            
            econ_final_matrix[k,l+2]=np.sum(econ_value_year)
            
            
            #convert all the NAN values to zero so they can be summed up together .
            
            #For the next loop reduce the efficiency of PV by degradation rate
            #median_eff_current_year=(median_eff_current_year-(annual_deg*(l+1)*median_eff_current_year)/100) # Vestigial, from when annual_deg was in % instead of decimal in beginning variables
            #multi_eff_current_year=(multi_eff_current_year-(annual_deg*(l+1)*multi_eff_current_year)/100)# Vestigial, from when annual_deg was in % instead of decimal in beginning variables
            median_eff_current_year=(median_eff_current_year-(annual_deg*(l+1)*median_eff_current_year))
            multi_eff_current_year=(multi_eff_current_year-(annual_deg*(l+1)*multi_eff_current_year))
            
            #pv_code value of 1 means it is a fixed PV otherwise it is a single axis tracker PV
            median_final_matrix[k,0]=pv_code       
            multi_final_matrix[k,0]=pv_code
            
            #Store the area of the PV in the second column
            median_final_matrix[k,1]=pv_area
            multi_final_matrix[k,1]=pv_area
            
            #Final pv generation values are in kWh and stored in the final median and multi PV power  matrices 
            median_final_matrix[k,l+2]=np.sum(median_pv_power_matrix)
            multi_final_matrix[k,l+2]=np.sum(multi_pv_power_matrix)
            
            
            # Print progress 
            progress = (k_indx/numArrays)*100
            print(beginning_year, 'array model progress in', current_year, ':', round(progress, 2), '%')
            
            #After the loop has gone through this year , increase the current year's value by 1
            current_year+=1
            
            
            #Testing code
            #print("The panel area for this point is:", pv_area)
            #print("The corresponding utility code was", utility_key)
            #print("The panel was located in the following utility,", utility_name)
            #print("The type of panel was", pv_type)
            #print("The optimal tilt angle was", tilt)
            
                
    #Dump the final numpy matrices into a csv file to store the results
    #pd.DataFrame(median_final_matrix).to_csv("D:\Dropbox\My Research stuff\INFEWS project related stuff\Data related to Jake project\Code run results\median_PV_matrix.csv",header=None, index=None)
    #pd.DataFrame(multi_final_matrix).to_csv("D:\Dropbox\My Research stuff\INFEWS project related stuff\Data related to Jake project\Code run results\multi_PV_matrix.csv", header=None, index=None)
    #pd.DataFrame(econ_final_matrix).to_csv(in_dir + out_dir + "\economic_value_" + str(beginning_year) + ".csv", header=None, index=None)
    return pd.DataFrame(econ_final_matrix).to_csv(in_dir + out_dir + "\economic_value_" + str(beginning_year) + ".csv", header=None, index=None)

"################################################################################################################################################### -- Automate across install years"

# Run the function for each year
getAgrisolarGenNEM(2018)
getAgrisolarGenNEM(2017)
getAgrisolarGenNEM(2016)
getAgrisolarGenNEM(2015)
getAgrisolarGenNEM(2014)
getAgrisolarGenNEM(2013)
getAgrisolarGenNEM(2012)
getAgrisolarGenNEM(2011)
getAgrisolarGenNEM(2010)
getAgrisolarGenNEM(2009)
getAgrisolarGenNEM(2008)

#Save final time
elapsed_time = (time.time() - t0) / 3600
print('Execution time:', round(elapsed_time, 2), 'hours')
