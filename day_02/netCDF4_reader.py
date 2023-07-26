# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 11:39:23 2023

netCDF4 Reader
"""

__author__ = 'Wyatt Spies'
__email__ = 'wyatt.spies@colorado.edu'

import netCDF4 as nc

dataset = nc.Dataset('\\Users\wyatt\Documents\swsss2023\day_02\wfs.t06z.ipe05.20230725_112000.nc')

dataset['tec'][:]      #how to get numpy array of the data
dataset['tec'].units   #how you get units of data

print(dataset)
print(dataset['tec'][:])
print(dataset['tec'].units)