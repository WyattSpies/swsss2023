# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 2023

netCDF4 Reader
"""

__author__ = 'Wyatt Spies'
__email__ = 'wyatt.spies@colorado.edu'

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np


#function to plot TEC on lat/lon grid
def plot_tec(filename, figsize=(12,6)):
    
    #read in data     //Users/wyatt/Documents/swsss2023/day_04/
    # dir = 'cdfData/'
    data = nc.Dataset(filename)
    
    #X = latitude, Y = longitude, Z = total electron counts (TEC)
    X, Y = np.meshgrid(data['lon'],data['lat'])
    Z = data['tec']
        
    #Plot TEC
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=figsize)
      
    plt.pcolormesh(X, Y, Z)
    plt.xlabel("Longitude: $($" + str(data['lon'].units) + '$)$')
    plt.ylabel("Latitude: $($" + str(data['lat'].units) + '$)$')
    plt.title("Total Electron Count", fontsize=15)
    plt.colorbar(label = str(data['tec'].units))
    
    return fig, ax

#function to save fig
def save_fig(infilename):
    
    #call plot function
    plot_tec(infilename)
    
    #Name png file
    outfilename = infilename + '.png'
    plt.savefig(outfilename, dpi=150)
    
    return














