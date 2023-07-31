#%%
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 15:43:49 2023

ASCII Reader for project 1
"""
__author__ = 'Wyatt Spies'
__email__ = 'wyatt.spies@colorado.edu'

import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from datetime import datetime

"""
Create a function to open ASCII files

Find a storm so we can compare JB2008 and TIEGCM densities

"""

#Create function to identify storm periods (DST < -100)
def storm_finder(filename, ind, starttime, endtime, storm_threshold):
    #Read date info and data
    with open(filename) as dat:
    
        #initialize data
        year = []
        day = []
        hour = []
        minute = []
        symh = []
        times = []
        
        down_time = []
        up_time = []
        
        storm_count = 0

        for line in dat:
            
            #skip line
            #print(dat.readline())
            
            tmp = line.split()
        
            #store variables to lists
            year.append(int(tmp[0]))
            day.append(int(tmp[1]))
            hour.append(int(tmp[2]))
            minute.append(int(tmp[3]))
            symh.append(int(tmp[ind]))
        
            if  len(year) >= 2: #search for 12 hour difference
                if symh[-1] <= storm_threshold and symh[-2] > storm_threshold:
                    down_time.append(int(tmp[1])*24*60 + int(tmp[1])*60 + int(tmp[3]))
                    
                    if len(up_time) == 0: 
                        storm_count = storm_count + 1
                        #print("start",storm_count, day[-1], hour[-1], minute[-1])
                        
                    elif (down_time[-1] - up_time[-1]) >= 720:
                        storm_count = storm_count + 1
                        #print("start",storm_count, day[-1], hour[-1], minute[-1])
                    
                if symh[-1] > storm_threshold and symh[-2] <= storm_threshold:
                    up_time.append(int(tmp[1])*24*60 + int(tmp[1])*60 + int(tmp[3]))
                    #print("stop",storm_count, day[-1], hour[-1], minute[-1])
                    
                    
            #convert from days to months
            datetime1 = dt.datetime(int(tmp[0]), 1, 1, int(tmp[2]), int(tmp[3])) + dt.timedelta(days=int(tmp[1])-1)
            #print(datetime1)
            times.append(datetime1)
            
    
        #Time to slice the data resolution to every hour
        times = times[0::60]
        symh = symh[0::60]
        #print(times)


    
            
            
        #Put all data and times into library
        all_data = {'Times':times,
                'SYMH':symh}
        
        #select important data range
        sel_data = {'Select Times':times,
                'Select SYMH':symh}
        
        #select stormy data range
        storm_data = {'Storm Times':times,
                'Storm SYMH':symh}

        #select specific time range
        time = np.array(all_data['Times'])
        print(np.shape(time))
        lp = (time>starttime)&(time<endtime)  #lp is a bool array
        time_sel = time[lp]
        
        symh = np.array(all_data['SYMH'])
        symh_sel = symh[lp]
        
        #include if start/stop time
        sel_data['Select Times'] = time_sel
        sel_data['Select SYMH'] = symh_sel

        #Identify storm periods
        #select specific time range
        symh = np.array(sel_data['Select SYMH'])
        lpS = (symh<storm_threshold)  #lp is a bool array
        symh_storm = symh[lpS]
        
        time = np.array(sel_data['Select Times'])
        time_storm = time[lpS]
        
        storm_data['Storm Times'] = time_storm
        storm_data['Storm SYMH'] = symh_storm
        
        

        #print(storm_count)
        return all_data, storm_data, storm_count, sel_data


    
 
    
 ##download 2002 and make and algorithm to find how many storms (check below 100 dst)
    
 
#Lets ID the files
directory = '/Users/wyatt/Documents/swsss2023/day_05/'
#Read data
dat_file = 'omniweb_2002.lst.txt'
#dat_file = 'omni_test.lst'
dat_location = directory + dat_file     # 

#data index  
index = -1

#choose start and stop times 
time1 = dt.datetime(2002, 4, 19, 12)
time2 = dt.datetime(2002, 4, 20, 23)


#call read function
data, stormy_data, counter, choose_data  = storm_finder(dat_location, index, time1, time2, -100)   #

counter = str(counter)

#plot data    
#plt.plot(data['Times'], data['SYMH'], 'blue', linewidth=.75)   #Plot
plt.plot(choose_data['Select Times'], choose_data['Select SYMH'], 'blue', linewidth=.75)
plt.plot(stormy_data['Storm Times'], stormy_data['Storm SYMH'], 'red', linewidth=.75)   #Plot
 
plt.grid(which='major', axis='both', color = 'coral')   #grid
plt.xlabel('Date')      #xlabel
plt.xticks(rotation=45)    #Rotate x axis labels to fit
plt.ylabel('SYMH $[nT]$')   #ylabel
plt.title('Geomagnetic Storms 2002: '+ counter)    #title
plt.show


#%%
"""
Plot and compare densities: Create a function to plot for every timestep
"""
# Import required packages
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat
import h5py
from datetime import datetime



######### JB2008 Data #########
dir_density_Jb2008 = 'C:/Users/wyatt/Documents/swsss2023/day_03/Data/JB2008/2002_JB2008_density.mat'

# Uses key to extract our data of interest
loaded_data = loadmat(dir_density_Jb2008)
JB2008_dens = loaded_data['densityData']
# Before we can visualize our density data, we first need to generate the 
# discretization grid of the density data in 3D space. We will be using 
# np.linspace to create evenly sapce data between the limits.
localSolarTimes_JB2008 = np.linspace(0,24,24)
latitudes_JB2008 = np.linspace(-87.5,87.5,20)
altitudes_JB2008 = np.linspace(100,800,36)
nofAlt_JB2008 = altitudes_JB2008.shape[0]
nofLst_JB2008 = localSolarTimes_JB2008.shape[0]
nofLat_JB2008 = latitudes_JB2008.shape[0]

# We can also impose additional constratints such as forcing the values to be integers.
time_array_JB2008 = np.linspace(0,8759,5, dtype = int)

# For the dataset that we will be working with today, you will need to reshape 
# them to be lst x lat x altitude
JB2008_dens_reshaped = np.reshape(JB2008_dens,(nofLst_JB2008,nofLat_JB2008,
                                               nofAlt_JB2008,8760), order='F')


######### TIEGCM Data #########
loaded_data = h5py.File('C:/Users/wyatt/Documents/swsss2023/day_03/Data_TIEGCM/TIEGCM/2002_TIEGCM_density.mat')

#discretization
tiegcm_dens = (10**np.array(loaded_data['density'])*1000).T # convert from g/cm3 to kg/m3
altitudes_tiegcm = np.array(loaded_data['altitudes']).flatten()
latitudes_tiegcm = np.array(loaded_data['latitudes']).flatten()
localSolarTimes_tiegcm = np.array(loaded_data['localSolarTimes']).flatten()
nofAlt_tiegcm = altitudes_tiegcm.shape[0]
nofLst_tiegcm = localSolarTimes_tiegcm.shape[0]
nofLat_tiegcm = latitudes_tiegcm.shape[0]

# We will be using the same time index as before.
time_array_tiegcm = time_array_JB2008

#resize to lst x lat x alt
tiegcm_dens_reshaped = np.reshape(tiegcm_dens,(nofLst_tiegcm,nofLat_tiegcm,nofAlt_tiegcm,8760), order='F')


######### Plot Data #########
#print('TIE-GCM altitude:', altitudes_tiegcm)
#print('JB2008 altitude:', altitudes_JB2008)


starttime_index = 109*24+12
endtime_index = starttime_index + 35
time_range = np.arange(starttime_index, endtime_index, 1)


#Plot for every time step
for time_index in time_range:
    
    day_index = time_index // 24
    hour_index = time_index % 24
    year = "2002"
    
    date_format = datetime.strptime(year + '-' + str(day_index) + " " + str(hour_index), "%Y-%j %H")
    #print(datetime.datetime.strptime(str(day_index) + " " + str(hour_index), "%d %H"))


    #generate 3d interpolate for JB2008
    JB2008_interpolator = RegularGridInterpolator((localSolarTimes_JB2008, latitudes_JB2008, altitudes_JB2008),
                                                  JB2008_dens_reshaped[:,:,:,time_index], bounds_error=False, fill_value=None)
    
    #idk what the hell this is? Make a grid in opposite space and evaluate JB2008 interp in this space at 400km
    JB_TIE_grid = np.zeros([len(localSolarTimes_tiegcm),len(latitudes_tiegcm)])
    for lst_i in range(len(localSolarTimes_tiegcm)):
        for lat_i in range(len(latitudes_tiegcm)):
            JB_TIE_grid[lst_i,lat_i] = JB2008_interpolator((localSolarTimes_tiegcm[lst_i], latitudes_tiegcm[lat_i],450))
    
    alt = 450
    hi2 = np.where(altitudes_tiegcm ==alt)
    hi1 = np.where(altitudes_JB2008 ==alt)
    
    #plot
    fig, axs = plt.subplots((4), figsize=(20, 20), sharex=True)
    
    #plot TIEGCM
    cs1 = axs[0].contourf(localSolarTimes_tiegcm, latitudes_tiegcm, JB_TIE_grid.T)
    axs[0].set_title('JB2008 Density at 450 km: ' + str(date_format), fontsize=18)
    axs[0].set_ylabel("Latitudes", fontsize=18)
    axs[0].tick_params(axis = 'both', which = 'major', labelsize = 16)
    
    #plot JB2008
    cs2 = axs[1].contourf(localSolarTimes_tiegcm, latitudes_tiegcm, 
                          tiegcm_dens_reshaped[:,:,hi2,time_index].squeeze().T)
    axs[1].set_title('TIEGCM Density at 450 km', fontsize=18)
    axs[1].set_ylabel("Latitudes", fontsize=18)
    axs[1].tick_params(axis = 'both', which = 'major', labelsize = 16)
    
    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig.colorbar(cs1,ax = axs[0])
    cbar.ax.set_ylabel('Density')
    
    cbar = fig.colorbar(cs2,ax = axs[1])
    cbar.ax.set_ylabel('Density')
    
    
    """
    Now, let's find the difference between both density models and plot out this 
    difference in a contour plot.
    """
    cs3 = axs[2].contourf(localSolarTimes_tiegcm, latitudes_tiegcm,
                          (JB_TIE_grid.T - tiegcm_dens_reshaped[:,:,hi2,time_index].squeeze().T))
    axs[2].set_title('Density Difference at 450 km', fontsize=18)
    axs[2].set_ylabel("Latitudes", fontsize=18)
    axs[2].tick_params(axis = 'both', which = 'major', labelsize = 16)
    
    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig.colorbar(cs3,ax = axs[2])
    cbar.ax.set_ylabel('Density')
    
    
    """
    APE = abs(tiegcm_dens-JB2008_dens)/tiegcm_dens
    """
    cs4 = axs[3].contourf(localSolarTimes_tiegcm, latitudes_tiegcm,
                          abs(JB_TIE_grid.T - tiegcm_dens_reshaped[:,:,hi2,time_index].squeeze().T)/JB_TIE_grid.T)
    axs[3].set_title('Density APE at 450 km', fontsize=18)
    axs[3].set_ylabel("Latitudes", fontsize=18)
    axs[3].tick_params(axis = 'both', which = 'major', labelsize = 16)
    
    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig.colorbar(cs4,ax = axs[3])
    cbar.ax.set_ylabel('Density')
    
    axs[3].set_xlabel("Local Solar Time", fontsize=18)
    
    
    #Save plots
    outfilename = 'Diff_plot_' + str(time_index) + '.png'
    plt.savefig(outfilename, dpi=150)

        
        
      
#%%
"""
calculate the hourly mean absolute percentage difference
for any arbitrary altitude over the same time period

arg(alt ('float'))
"""

def mean_calc(alt):
    
    mean_dens_JB2008 = []
    mean_dens_TIEGCM = []
    
    for time_index in time_range:
        
        #generate 3d interpolate for JB2008
        JB2008_interpolator = RegularGridInterpolator((localSolarTimes_JB2008, latitudes_JB2008, altitudes_JB2008),
                                                      JB2008_dens_reshaped[:,:,:,time_index], bounds_error=False, fill_value=None)
        
        JB_TIE_grid = np.zeros([len(localSolarTimes_tiegcm),len(latitudes_tiegcm)])
        for lst_i in range(len(localSolarTimes_tiegcm)):
            for lat_i in range(len(latitudes_tiegcm)):
                JB_TIE_grid[lst_i,lat_i] = JB2008_interpolator((localSolarTimes_tiegcm[lst_i], latitudes_tiegcm[lat_i],alt))
        
        
        #generate 3d interpolate for JB2008
        tiegcm_interpolator = RegularGridInterpolator((localSolarTimes_tiegcm, latitudes_tiegcm, altitudes_tiegcm),
                                                  tiegcm_dens_reshaped[:,:,:,time_index], bounds_error=False, fill_value=None)
        
        tie_JB_grid = np.zeros([len(localSolarTimes_JB2008),len(latitudes_JB2008)])
        for lst_i in range(len(localSolarTimes_JB2008)):
            for lat_i in range(len(latitudes_JB2008)):
                tie_JB_grid[lst_i,lat_i] = tiegcm_interpolator((localSolarTimes_JB2008[lst_i],latitudes_JB2008[lat_i],alt))
        
        #                                       lat/lon/alt/time
        mean_dens_JB2008.append(np.mean(JB_TIE_grid.T))
        mean_dens_TIEGCM.append(np.mean(tie_JB_grid.T))

    mean_APE = (abs(np.array(mean_dens_JB2008) - np.array(mean_dens_TIEGCM))/np.array(mean_dens_TIEGCM))   
    x = time_range
    y1 = mean_APE
    
    plt.plot(x, y1, 'blue')
    plt.xlabel('Time Steps')
    plt.ylabel('Mean APE')
    plt.title('Mean Density comparison by Altitude: Alt = ' + str(alt) + 'km')
    plt.show
    
mean_calc(250)
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    