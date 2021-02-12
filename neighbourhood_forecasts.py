'Code written by Carlo Cafaro in October 2020'

'for any questions e-mail me: c.cafaro@reading.ac.uk'

import nbood_code as nb
import plot_maps as plm
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np





def calculate_fractions(leadtime,varname_fcst,varname_obs,threshold,n,percentyes,percent):
    
    
    #get the data first
    
    rain_fcst,lat_f,lon_f=nb.get_rainfall_forecast_data(leadtime,varname_fcst)
    
    rain_obs,lat_obs,lon_obs=nb.get_observations(leadtime,varname_obs)
    
    #then convert it to binary fields and generating fractions
    
    
    if percentyes==1:
       thresh_fcst=np.percentile(rain_fcst,percent)
       thresh_obs=np.percentile(rain_obs,percent)
       
       print (thresh_fcst,'forecast threshold')
       print (thresh_obs,'obs threshold')
        
       obs_fractions,fcst_fractions=nb.generating_fractions(rain_fcst,rain_obs,thresh_fcst,thresh_obs,n)
       
    else:
        threshold1=threshold
        threshold2=threshold    
        obs_fractions,fcst_fractions=nb.generating_fractions(rain_fcst,rain_obs,threshold1,threshold2,n) 
        
    
    #save to netcdf files
    
    nb.create_netcdf_files_observations(obs_fractions,leadtime,n)
    
    nb.create_netcdf_files_forecasts(fcst_fractions,leadtime,n)
    
    return lat_f,lon_f,lat_obs,lon_obs,fcst_fractions,obs_fractions



def calculate_fss(fcst_fractions,obs_fractions,threshold_fcst,threshold_obs,n,lata,latb,lona,lonb):
    
    
    fss=nb.neighbourhood_verification(fcst_fractions,obs_fractions,threshold_fcst,threshold_obs,n,lata,latb,lona,lonb)
    
    
    return fss


def plot_fractions(lat_f,lon_f,lat_obs,lon_obs,lat1,lat2,lon1,lon2,n):
    
    
    
    obs=nc.Dataset('/home/users/ccafaro/fractions_observations'+'nbood_size'+str(n)+'.nc')
    fcst=nc.Dataset('/home/users/ccafaro/fractions_forecast'+'nbood_size'+str(n)+'.nc')
    
    obs_fractions=obs.variables['fcn']
    
    fcst_fractions=fcst.variables['fcn']
    
    
    plm.plot_fractional_maps(lat_f,lon_f,lat_obs,lon_obs,fcst_fractions,obs_fractions,lat1,lat2,lon1,lon2)
    
    return
        





def main():
    
    #these are the three parameter to choose:
    
    #leadtime of the forecast
    #rainfall threshold
    #neighbourhood size
    
    leadtime=6  #the time index of your file
    
    n=20        #this is the size of the neighbourhood you want to build the forecast on
    
    #specify the extent of the region you want to calculate the fractions skill score on
    
    lat1=-3
    lat2=20
    lon1=-20
    lon2=20
    
    
    
    
    percentyes=0  #if 1 means you use percentiles (defined below), if 0 you use physical threshold
    percent=95
    threshold=20 #this is in mm, accumulation over the time period you have the forecast on
    
    varname_obs=#specify the variable name of the observed field
    varname_fcst=#specify the variable name of the forecast field
    
    #this will generate the forecast fractions and the observed fractions
    
    lat_f,lon_f,lat_obs,lon_obs,fcst_fractions,obs_fractions=calculate_fractions(leadtime,varname_obs,varname_fcst,threshold,n,percentyes,percent)
    
    
    fss=calculate_fss(fcst_fractions,obs_fractions,threshold,threshold,n,lat1,lat2,lon1,lon2)
    
    
    
    #plot_fractions(lat_f,lon_f,lat_obs,lon_obs,lat1,lat2,lon1,lon2,n)
    
    print ('FSS=',fss)
    
if __name__ == "__main__":
    
    main()
