# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 16:02:15 2023

Auroral Birthday Plot
"""
__author__ = 'wyatt spies'
__email__ = 'wyatt.spies@colorado.edu'


import matplotlib.pyplot as plt
from datetime import datetime
from swmfpy.web import get_omni_data

start_time = datetime(1999, 10, 26)
end_time = datetime(1999, 10, 27)
data = get_omni_data(start_time, end_time)
data.keys()

x = data['times']       #Time steps during my birthday
#print(x)
y = data['al']          #Auroral data during my birthday
#print(y)

plt.plot(x, y, 'coral')   #Plot y vs x with color "coral"
plt.xlabel('Time $[UTC]$')
plt.xticks(rotation=45)    #Rotate x axis labels to fit
plt.ylabel('AL $[nT]$')
plt.title('Auroral Data for 10/26/1999')
plt.show