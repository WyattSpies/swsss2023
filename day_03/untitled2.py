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



#this is a function to open ASCII with inputs ("filename", index)
def read_ascii_file(filename, ind, starttime, endtime):         # include if wanted start/stop time    
    #Read date info and data
    with open(filename) as dat:
    
        #initialize data
        year = []
        day = []
        hour = []
        minute = []
        symh = []
        times = []
        
        

        for line in dat:
            
            tmp = line.split()
        
            #store variables to lists
            year.append(int(tmp[0]))
            day.append(int(tmp[1]))
            hour.append(int(tmp[2]))
            minute.append(int(tmp[3]))
            symh.append(int(tmp[ind]))
            
            #convert from days to months
            datetime1 = dt.datetime(int(tmp[0]), 1, 1, int(tmp[2]), int(tmp[3])) + dt.timedelta(days=int(tmp[1])-1)
            print(datetime1)
            times.append(datetime1)
            
            
        #Put all data and times into library
        all_data = {'Times':times,
                'SYMH':symh}
        #select important data range
        sel_data = {'Select Times':times,
                'Select SYMH':symh}
        
        #select specific time range
        time = np.array(all_data['Times'])
        lp = (time>starttime)&(time<endtime)  #lp is a bool array
        time_sel = time[lp]
        
        symh = np.array(all_data['SYMH'])
        symh_sel = symh[lp]
        
        #include if start/stop time
        sel_data['Select Times'] = time_sel
        sel_data['Select SYMH'] = symh_sel
        
        
        return all_data, sel_data       #include if start/stop time  , sel_data
    
 
#--- This is the main code ------------
# if __name__ == "__main__":
#     dat_location = '/Users/wyatt/Documents/swsss2023/day_02/omniweb_data.lst'
#     index = -1
#     time1 = (2013, 3, 16)
#     time2 = (2013, 3, 21)
#     data, choose_data = read_ascii_file(dat_location, index, time1, time2)
    

    
 
#Lets ID the files
directory = '/Users/wyatt/Documents/swsss2023/day_02/'
#Read data
dat_file = 'omniweb_data.lst'
#dat_file = 'omni_test.lst'
dat_location = directory + dat_file  

#data index  
index = -1

#choose start and stop times 
time1 = dt.datetime(2013, 3, 17)
time2 = dt.datetime(2013, 3, 19)


#call read function
data, choose_data = read_ascii_file(dat_location, index, time1, time2)   

#plot data    
plt.plot(data['Times'], data['SYMH'], 'blue', linewidth=.75)   #Plot
plt.plot(choose_data['Select Times'], choose_data['Select SYMH'], 'red', linewidth=.75)   #Plot
plt.grid(which='major', axis='both', color = 'coral')   #grid
plt.xlabel('Date')      #xlabel
plt.xticks(rotation=45)    #Rotate x axis labels to fit
plt.ylabel('SYMH $[nT]$')   #ylabel
plt.title('Geomagnetic Storms 03/17/2013: ')    #title
plt.show