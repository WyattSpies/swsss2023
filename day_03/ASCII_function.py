# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:43:49 2023

ASCII Reader
"""
__author__ = 'Wyatt Spies'
__email__ = 'wyatt.spies@colorado.edu'

import matplotlib.pyplot as plt
import numpy as np
import datetime as dt

"""
Create a function to open ASCII files

"""

#Create function to identify storm periods (DST < -100)
def storm_finder(filename, ind, storm_threshold):
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
                        print("start",storm_count, day[-1], hour[-1], minute[-1])
                        
                    elif (down_time[-1] - up_time[-1]) >= 720:
                        storm_count = storm_count + 1
                        print("start",storm_count, day[-1], hour[-1], minute[-1])
                    
                if symh[-1] > storm_threshold and symh[-2] <= storm_threshold:
                    up_time.append(int(tmp[1])*24*60 + int(tmp[1])*60 + int(tmp[3]))
                    print("stop",storm_count, day[-1], hour[-1], minute[-1])
                    
                    
            #convert from days to months
            datetime1 = dt.datetime(int(tmp[0]), 1, 1, int(tmp[2]), int(tmp[3])) + dt.timedelta(days=int(tmp[1])-1)
            #print(datetime1)
            times.append(datetime1)
            
            
            
        #Put all data and times into library
        all_data = {'Times':times,
                'SYMH':symh}
        #select important data range
        storm_data = {'Storm Times':times,
                'Storm SYMH':symh}

        #Identify storm periods
        #select specific time range
        symh = np.array(all_data['SYMH'])
        lp = (symh<storm_threshold)  #lp is a bool array
        symh_storm = symh[lp]
        
        time = np.array(all_data['Times'])
        time_storm = time[lp]
        
        storm_data['Storm Times'] = time_storm
        storm_data['Storm SYMH'] = symh_storm
        
        print(storm_count)
        return all_data, storm_data, storm_count


    
 
    
 ##download 2003 and make and algorithm to find how many storms (check below 100 dst)
    
 
#Lets ID the files
directory = '/Users/wyatt/Documents/swsss2023/day_03/'
#Read data
dat_file = 'omniweb_2003.lst'
#dat_file = 'omni_test.lst'
dat_location = directory + dat_file  

#data index  
index = -1

#choose start and stop times 
#time1 = dt.datetime(2003, 3, 17)
#time2 = dt.datetime(2003, 3, 18)


#call read function
data, choose_data, counter = storm_finder(dat_location, index, -100)   #, time1, time2

counter = str(counter)

#plot data    
#plt.plot(data['Times'], data['SYMH'], 'blue', linewidth=.75)   #Plot
plt.plot(choose_data['Storm Times'], choose_data['Storm SYMH'], 'red', linewidth=.75)   #Plot

plt.grid(which='major', axis='both', color = 'coral')   #grid
plt.xlabel('Date')      #xlabel
plt.xticks(rotation=45)    #Rotate x axis labels to fit
plt.ylabel('SYMH $[nT]$')   #ylabel
plt.title('Geomagnetic Storms 2003: '+ counter)    #title
plt.show


        
        
        

        
        
        
        
        
        
        