# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:43:49 2023

ASCII Reader
"""
__author__ = 'Wyatt Spies'
__email__ = 'wyatt.spies@colorado.edu'

import matplotlib.pyplot as plt
import datetime as dt

"""
open ASCII files

Format file for variable descriptions
     omniweb_format.fmt   

    ITEMS                      FORMAT   
     
 1 Year                          I4        
 2 Day                           I4        
 3 Hour                          I3        
 4 Minute                        I3        
 5 SYM/H, nT                     I6   


Data file
    omniweb_data.lst
"""

#Lets open the files
directory = '/Users/wyatt/Documents/swsss2023/day_02/'

#Read data
dat_file = 'omniweb_data.lst '
#dat_file = 'omni_test.lst'
dat_location = directory + dat_file

#Read date info and data
with open(dat_location) as dat:
    
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
        symh.append(int(tmp[4]))
        
        #convert from days to months
        datetime1 = dt.datetime(int(tmp[0]), 1, 1, int(tmp[2]), int(tmp[3])) + dt.timedelta(days=int(tmp[1])-1)
        print(datetime1)
        times.append(datetime1)
        
#plot data        
plt.plot(times, symh, 'teal')   #Plot
plt.grid(which='major', axis='both', color = 'coral')   #grid
plt.xlabel('Date', color = 'red')      #xlabel
plt.xticks(rotation=45)    #Rotate x axis labels to fit
plt.ylabel('SYMH $[nT]$', color = 'yellow')   #ylabel
plt.title('Geomagnetic Storm on 03/17/2013', color = 'green')    #title
plt.show
        
        
        

        
        
        
        
        
        
        