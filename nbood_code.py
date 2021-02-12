
'Code written by Carlo Cafaro in October 2020'

'for any questions e-mail me: c.cafaro@reading.ac.uk'

import numpy as np
from scipy import signal
import netCDF4 as nc
from numpy import dtype





#THIS IS FOR STEPS 1-2



#First of all we get the observations

def create_netcdf_files_observations(fractions,leadtime,varname,n):
    
    
    precip,latitude,longitude=get_observations(leadtime,varname)
    
    
    fn = '/home/users/ccafaro/fractions_observations'+'nbood_size'+str(n)+'.nc'
    ds = nc.Dataset(fn, 'w', format='NETCDF4')

    time = ds.createDimension('time', None)
    lats = ds.createDimension('lat', len(fractions[:,0]))
    lons = ds.createDimension('lon', len(fractions[0,:]))

     # create latitude axis
    lats = ds.createVariable('lat', dtype('double').char, ('lat'))
    lats.standard_name = 'latitude'
    lats.long_name = 'latitude'
    lats.units = 'degrees_north'
    lats.axis = 'Y'

    # create longitude axis
    lons = ds.createVariable('lon', dtype('double').char, ('lon'))
    lons.standard_name = 'longitude'
    lons.long_name = 'longitude'
    lons.units = 'degrees_east'
    lons.axis = 'X'

    value = ds.createVariable('fcn', dtype('double').char, ('time', 'lat', 'lon'))
    value.long_name = 'fractions'
    value.units = 'unknown'



    
    lats[:]= latitude[:]
    lons[:]= longitude[:]


    value[0, :, :] = fractions

    

    ds.close()
    
    return
    
def sub_domain(lata, latb, lona, lonb,varname):
    'here we specifiy the domain of our study, retrieving latitude and longitude indices from convperm, global models and gpm observations'

    
    precip,latobs,lonobs=get_observations(0,varname)
    

    indexlatobs = [i for i, e in enumerate(latobs) if lata <= e <= latb]
    indexlonobs = [i for i, e in enumerate(lonobs) if lona <= e <= lonb]

    

    return indexlatobs, indexlonobs





def create_netcdf_files_forecasts(fractions,leadtime,n):
    
    
    precip,latitude,longitude=get_rainfall_forecast_data(leadtime)
    
    
    fn = '/home/users/ccafaro/fractions_forecast'+'nbood_size'+str(n)+'.nc'
    ds = nc.Dataset(fn, 'w', format='NETCDF4')

    time = ds.createDimension('time', None)
    lats = ds.createDimension('lat', len(fractions[:,0]))
    lons = ds.createDimension('lon', len(fractions[0,:]))
    
    
   
    

     # create latitude axis
    lats = ds.createVariable('lat', dtype('double').char, ('lat'))
    lats.standard_name = 'latitude'
    lats.long_name = 'latitude'
    lats.units = 'degrees_north'
    lats.axis = 'Y'

    # create longitude axis
    lons = ds.createVariable('lon', dtype('double').char, ('lon'))
    lons.standard_name = 'longitude'
    lons.long_name = 'longitude'
    lons.units = 'degrees_east'
    lons.axis = 'X'
    
    
    value = ds.createVariable('fcn', dtype('double').char, ('time', 'lat', 'lon'))
    value.long_name = 'fractions'
    value.units = 'unknown'
    

    lats[:]= latitude[:]
    lons[:]= longitude[:]


    value[0, :, :] = fractions

    

    ds.close()
    
    return



def get_rainfall_forecast_data(leadtime,varname):
    
    dc=nc.Dataset('/your_path/wrfout_d03_2020-09-25_18_00_00.nc') #this need to be changed according to your file
    
    precip=dc.variables[varname][leadtime,:,:]
    
    lat=dc.variables['lat']
    
    lon=dc.variables['lon']
    
    return precip,lat,lon
    
def get_observations(leadtime,varname):
    
    
    dc=nc.Dataset('your_path/') #this could be GPM data
    
    precip=dc.variables[varname][leadtime,:,:]
    
    lat=dc.variables['lat']

    lon=dc.variables['lon']
    
    
    return precip,lat,lon
    



#THIS IS FOR STEPS 1-2
def produce_binary_fields(rainfall_field,threshold): #the threshold could be either physical or a percentile
    
    binary_rainfall=np.where(rainfall_field>=threshold,1,0) #where the thresholds is exceed you get 1, otherwise 0
    
    
    
    
    
    return binary_rainfall
    
    

    

#STEPS 3-5
def generating_fractions(rain_for,rain_obs,threshold_fcst,threshold_obs,n): #generating fractions applied to both forecast and observations on a given neighbourhood of size n and a given threshold
    
    
    #n must be an odd number (it is the number of grid points along the neighbourhood side
    
    binary_obs=produce_binary_fields(rain_obs,threshold_obs)
    binary_for=produce_binary_fields(rain_for,threshold_fcst)

    

    
        

    kernel = np.ones((n, n)) #this defines the neighbourhood size (a matrix of ones of size n x n
    
    mode='same'
    obs_fractions = fourier_filter(binary_obs, n, mode) / kernel.sum() #fourier filter function defined below is to compute the spatial mean within that neighbourhood

                
    prob_fractions = fourier_filter(binary_for, n, mode) / kernel.sum()
    
    
    
    return obs_fractions,prob_fractions
    
   
def fourier_filter(field, n, mode='same'): #this is the function to calculate fractions within a neighbourhood of size n applied to the binary field

    return signal.fftconvolve(field, np.ones((n, n)), mode)   
    
    
    

#THE FOLLOWING FUNCTIONS ARE TO PERFORM THE NEIGHBOURHOOD VERIFICATION
    
    
def fractions_skill_score(fcst_fractions,obs_fractions,threshold_fcst,threshold_obs,n,lata,latb,lona,lonb):  #this is to compute the metrics applied to the fractions
    
    
    indexlatobs, indexlonobs=sub_domain(lata, latb, lona, lonb)
    
    latmin=indexlatobs[0]
    latmax=indexlatobs[len(indexlatobs)-1]+1
    
    lonmin=indexlonobs[0]
    lonmax=indexlonobs[len(indexlonobs)-1]+1
    
    
    fcst_fractions=fcst_fractions[latmin:latmax,lonmin:lonmax]
    
    obs_fractions=obs_fractions[latmin:latmax,lonmin:lonmax]
    
    
    num = np.nanmean(np.power(fcst_fractions - obs_fractions, 2))

    denom = np.nanmean(np.power(fcst_fractions,2) + np.power(obs_fractions,2))

    
    
    
    return num,denom
    

def neighbourhood_verification(fcst_fractions,obs_fractions,threshold_fcst,threshold_obs,n,lata,latb,lona,lonb):
    
     #this defines the number of grid-points for the neighbourhood side, it's an odd number (1,3,5,7,...)
    
    
    #IN THIS CASE rain_fcst should be regridded to rain_obs first.
    
    
        
    num,denom=fractions_skill_score(fcst_fractions,obs_fractions,threshold_fcst,threshold_obs,n,lata,latb,lona,lonb)
        
    fss=1-num/denom
        

    return fss        
    


