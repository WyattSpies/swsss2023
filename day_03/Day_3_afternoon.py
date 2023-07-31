#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Welcome to Space Weather Simulation Summer School Day 3

Today, we will be working with various file types, doing some simple data 
manipulation and data visualization

We will be using a lot of things that have been covered over the last two days 
with minor vairation.

Goal: Getting comfortable with reading and writing data, doing simple data 
manipulation, and visualizing data.

Task: Fill in the cells with the correct codes

@author: Peng Mun Siew
"""

#%% 
"""
This is a code cell that we can use to partition our data. (Similar to Matlab cell)
We hace the options to run the codes cell by cell by using the "Run current cell" button on the top.
"""
print ("Hello World")

#%%
"""
Creating a random numpy array
"""
# Importing the required packages
import numpy as np

# Generate a random array of dimension 10 by 5
data_arr = np.random.randn(10,5)
print(data_arr)

#%%
"""
TODO: Writing and reading numpy file
"""
# Save the data_arr variable into a .npy file
np.save('test_np_save.npy', data_arr)

# Load data from a .npy file
data_arr_loaded = np.load('test_np_save.npy')

# Verify that the loaded data matches the initial data
print(np.equal(data_arr, data_arr_loaded))
print(data_arr==data_arr_loaded)

#%%
"""
TODO: Writing and reading numpy zip archive/file
"""
# Generate a second random array of dimension 8 by 1
data_arr2 = np.random.randn(8,1)
print(data_arr2)

# Save the data_arr and data_arr2 variables into a .npz file
np.savez('test_savez.npz', data_arr, data_arr2)

# Load the numpy zip file
npzfile = np.load('test_savez.npz')
print('Variable names within this file: ', sorted(npzfile.files))
print(npzfile['arr_0'])

# Verify that the loaded data matches the initial data
print((data_arr==npzfile['arr_0']).all())
print((data_arr2==npzfile['arr_1']).all())

#%%
"""
Error and exception
"""
# Exception handling, can be use with assertion as well
try:
    # Python will try to execute any code here, and if there is an exception 
    # skip to below 
    print(np.equal(data_arr,npzfile).all())
except:
    # Execute this code when there is an exception (unable to run code in try)
    print("The codes in try returned an error.")
    print(np.equal(data_arr,npzfile['arr_0']).all())
    
#%%
"""
TODO: Error solving 1
"""
# What is wrong with the following line? 
np.equal(data_arr,data_arr2)
#incorrect array shapes
#%%
"""
TODO: Error solving 2
"""
# What is wrong with the following line? 
np.equal(data_arr2,npzfile['data_arr2'])
#wrong key reference

#%%
"""
TODO: Error solving 3
"""
# What is wrong with the following line? 
#numpy.equal(data_arr2,npzfile['arr_1'])
#numpy should be np
np.equal(data_arr2,npzfile['arr_1'])

#%%
"""
Loading data from Matlab
"""

# Import required packages
import numpy as np
from scipy.io import loadmat

dir_density_Jb2008 = 'Data/JB2008/2002_JB2008_density.mat'

# Load Density Data
#try:
loaded_data = loadmat(dir_density_Jb2008)
print (loaded_data)
#except:
 #   print("File not found. Please check your directory")

# Uses key to extract our data of interest
JB2008_dens = loaded_data['densityData']

# The shape command now works
print(JB2008_dens.shape)

#%%
"""
Data visualization I

Let's visualize the density field for 400 KM at different time.
"""
# Import required packages
import matplotlib.pyplot as plt

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
                                               nofAlt_JB2008,8760), order='F') # Fortran-like index order

#%%
"""
TODO: Plot the atmospheric density for 400 KM for the first time index in
      time_array_JB2008 (time_array_JB2008[0]).
"""

import matplotlib.pyplot as plt

# Look for data that correspond to an altitude of 400 KM
alt = 400
hi = np.where(altitudes_JB2008==alt)
ik = 4
#Time to plot
fig, axs = plt.subplots(1, figsize=(15, 4), sharex=True)

#       contour         LST                latitude         density at all lat/lon, 400km, time=0
cs = axs.contourf(localSolarTimes_JB2008,latitudes_JB2008,JB2008_dens_reshaped[:,:,hi,time_array_JB2008[0]].squeeze().T)
axs.set_title('JB2008 density at 400 km, t = {} hrs'.format(time_array_JB2008[ik]), fontsize=18)
axs.set_ylabel("Latitudes", fontsize=18)
axs.tick_params(axis = 'both', which = 'major', labelsize = 16)
    
# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig.colorbar(cs,ax=axs)
cbar.ax.set_ylabel('Density')

axs.set_xlabel("Local Solar Time", fontsize=18)  

#%%
"""
TODO: Plot the atmospheric density for 300 KM for all time indexes in
      time_array_JB2008
"""
import matplotlib.pyplot as plt

# Look for data that correspond to an altitude of 400 KM
alt = 300
hi = np.where(altitudes_JB2008==alt)

# Create a canvas to plot our data on. Here we are using a subplot with 5 spaces for the plots.
fig, axs = plt.subplots((5), figsize=(15, 10*2), sharex=True)

for ik in range (5):
    cs = axs[ik].contourf(localSolarTimes_JB2008, latitudes_JB2008, JB2008_dens_reshaped[:,:,hi,time_array_JB2008[ik]].squeeze().T)
    axs[ik].set_title('JB2008 density at 300 km, t = {} hrs'.format(time_array_JB2008[ik]), fontsize=18)
    axs[ik].set_ylabel("Latitudes", fontsize=18)
    axs[ik].tick_params(axis = 'both', which = 'major', labelsize = 16)
    
    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig.colorbar(cs,ax = axs[ik])
    cbar.ax.set_ylabel('Density')

axs[ik].set_xlabel("Local Solar Time", fontsize=18) 
#%%
"""
Assignment 1

Can you plot the mean density for each altitude at February 1st, 2002?
"""

# First identidy the time index that corresponds to  February 1st, 2002. 
# Note the data is generated at an hourly interval from 00:00 January 1st, 2002
time_index = 31*24
dens_data_feb1 = JB2008_dens_reshaped[:,:,:,time_index]
print('The dimension of the data are as followed (local solar time,latitude,altitude):', dens_data_feb1.shape)

mean_dens = [np.mean(dens_data_feb1[:,:,ik]) for ik in range(len(altitudes_JB2008))]

plt.semilogy(altitudes_JB2008, mean_dens)
plt.xlabel("Altitude $[km]$")
plt.ylabel("Density $[NA]$")
plt.title("Mean Density vs Altitude")
plt.show

#%%
"""
Data Visualization II

Now, let's us work with density data from TIE-GCM instead, and plot the density 
field at 310km

"""
# Import required packages
import h5py
loaded_data = h5py.File('Data_TIEGCM/TIEGCM/2002_TIEGCM_density.mat')

# This is a HDF5 dataset object, some similarity with a dictionary
print('Key within dataset:',list(loaded_data.keys()))

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

#%%
"""
TODO: Plot the atmospheric density for 310 KM for all time indexes in
      time_array_tiegcm
"""
import matplotlib.pyplot as plt

# Look for data that correspond to an altitude of 400 KM
alt = 310
hi = np.where(altitudes_tiegcm==alt)

# Create a canvas to plot our data on. Here we are using a subplot with 5 spaces for the plots.
fig, axs = plt.subplots((5), figsize=(15, 10*2), sharex=True)

for ik in range (5):
    cs = axs[ik].contourf(localSolarTimes_tiegcm, latitudes_tiegcm, tiegcm_dens_reshaped[:,:,hi,time_array_tiegcm[ik]].squeeze().T)
    axs[ik].set_title('JB2008 density at 310 km, t = {} hrs'.format(time_array_tiegcm[ik]), fontsize=18)
    axs[ik].set_ylabel("Latitudes", fontsize=18)
    axs[ik].tick_params(axis = 'both', which = 'major', labelsize = 16)
    
    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig.colorbar(cs,ax = axs[ik])
    cbar.ax.set_ylabel('Density')

axs[ik].set_xlabel("Local Solar Time", fontsize=18) 
#%%
print(JB2008_dens_reshaped.shape, tiegcm_dens_reshaped.shape)

#%%
"""
Assignment 1.5

Can you plot the mean density for each altitude at February 1st, 2002 for both 
models (JB2008 and TIE-GCM) on the same plot?
"""

# First identidy the time index that corresponds to  February 1st, 2002. 
# Note the data is generated at an hourly interval from 00:00 January 1st, 2002
time_index = 31*24


#%%
"""
Data Interpolation (1D)

Now, let's us look at how to do data interpolation with scipy
"""
# Import required packages
from scipy import interpolate
import matplotlib.pyplot as plt

# Let's first create some data for interpolation
x = np.arange(0, 10)
y = np.exp(-x/3.0)

plt.scatter(x,y, color = "blue")
# Generate 1D interpolant function
interp_func_1D = interpolate.interp1d(x, y)

# Let's select some new points
xnew = np.arange(0, 9, 0.1)
# use interpolation function returned by `interp1d`
ynew = interp_func_1D(xnew)
#plot linear interpolation
plt.plot(xnew,ynew, '--', color = "red")


interp_func_1D_cubic = interpolate.interp1d(x, y,kind='cubic')
# use interpolation function returned by `interp1d`
ycubic = interp_func_1D_cubic(xnew) 
#plot cubic interpolation
plt.plot(xnew,ycubic, '-.', color = "green")


interp_func_1D_quadratic = interpolate.interp1d(x, y,kind='quadratic')
# use interpolation function returned by `interp1d`
yquadratic = interp_func_1D_quadratic(xnew) 
#plot quadratic interpolation
plt.plot(xnew,ycubic, color = "black")


#%%
"""
Data Interpolation (3D)

Now, let's us look at how to do data interpolation with scipy
"""
# Import required packages
from scipy.interpolate import RegularGridInterpolator

# First create a set of sample data that we will be using 3D interpolation on
def function_1(x, y, z):
    return 2 * x**3 + 3 * y**2 - z

x = np.linspace(1, 4, 11)
y = np.linspace(4, 7, 22)
z = np.linspace(7, 9, 33)
xg, yg ,zg = np.meshgrid(x, y, z, indexing='ij', sparse=True)

sample_data = function_1(xg, yg, zg)

# Generate Interpolant (interpolating function)
interpolated_function_1 = RegularGridInterpolator((x, y, z), sample_data)

# Say we are interested in the points [[2.1, 6.2, 8.3], [3.3, 5.2, 7.1]]
pts = np.array([[2.1, 6.2, 8.3], [3.3, 5.2, 7.1]])
print('Using interpolation method:',interpolated_function_1(pts)) 
print('From true function:',function_1(pts[:,0],pts[:,1],pts[:,2]))

#%%
"""
Saving mat file

Now, let's us look at how to we can save our data into a mat file
"""
# Import required packages
from scipy.io import savemat

a = np.arange(20)
mdic = {"a": a, "label": "experiment"} # Using dictionary to store multiple variables
savemat("matlab_matrix.mat", mdic)

#%%
"""
Assignment 2 (a)

The two data that we have been working on today have different discretization 
grid.

Use 3D interpolation to evaluate the TIE-GCM density field at 400 KM on 
February 1st, 2002, with the discretized grid used for the JB2008 
((nofLst_JB2008,nofLat_JB2008,nofAlt_JB2008).
"""
print('TIE-GCM altitude:', altitudes_tiegcm)
#print('JB2008 altitude:', altitudes_JB2008)


time_index = 31*24
#generate 3d interpolate for tiegcm
tiegcm_interpolator = RegularGridInterpolator((localSolarTimes_tiegcm, latitudes_tiegcm, altitudes_tiegcm),
                                              tiegcm_dens_reshaped[:,:,:,time_index], bounds_error=False, fill_value=None)


#make a meshgrid
#xnew = 400km
#tie_interp = tiegcm_interpolator()

#idk what the hell this is? Make a grid in opposite space and evaluate tiegcm interp in this space at 400km
tie_JB_grid = np.zeros([len(localSolarTimes_JB2008),len(latitudes_JB2008)])
for lst_i in range(len(localSolarTimes_JB2008)):
    for lat_i in range(len(latitudes_JB2008)):
        tie_JB_grid[lst_i,lat_i] = tiegcm_interpolator((localSolarTimes_JB2008[lst_i],latitudes_JB2008[lat_i],400))

alt = 400
hi1 = np.where(altitudes_tiegcm ==alt)
hi2 = np.where(altitudes_JB2008 ==alt)

#plot
fig, axs = plt.subplots((4), figsize=(20, 20), sharex=True)

#plot TIEGCM
cs1 = axs[0].contourf(localSolarTimes_JB2008, latitudes_JB2008, tie_JB_grid.T)
axs[0].set_title('TIEGCM Density at 400 km', fontsize=18)
axs[0].set_ylabel("Latitudes", fontsize=18)
axs[0].tick_params(axis = 'both', which = 'major', labelsize = 16)

#plot JB2008
cs2 = axs[1].contourf(localSolarTimes_JB2008, latitudes_JB2008, JB2008_dens_reshaped[:,:,hi2,time_index].squeeze().T)
axs[1].set_title('JB2008 Density at 400 km', fontsize=18)
axs[1].set_ylabel("Latitudes", fontsize=18)
axs[1].tick_params(axis = 'both', which = 'major', labelsize = 16)

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig.colorbar(cs1,ax = axs[0])
cbar.ax.set_ylabel('Density')

cbar = fig.colorbar(cs2,ax = axs[1])
cbar.ax.set_ylabel('Density')

#axs[0].set_xlabel("Local Solar Time", fontsize=18) 
#axs[1].set_xlabel("Local Solar Time", fontsize=18)

"""
Assignment 2 (b)

Now, let's find the difference between both density models and plot out this 
difference in a contour plot.
"""
cs3 = axs[2].contourf(localSolarTimes_JB2008, latitudes_JB2008, (tie_JB_grid.T - JB2008_dens_reshaped[:,:,hi2,time_index].squeeze().T))
axs[2].set_title('Density Difference at 400 km', fontsize=18)
axs[2].set_ylabel("Latitudes", fontsize=18)
axs[2].tick_params(axis = 'both', which = 'major', labelsize = 16)

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig.colorbar(cs3,ax = axs[2])
cbar.ax.set_ylabel('Density')

#axs[2].set_xlabel("Local Solar Time", fontsize=18)
"""
Assignment 2 (c)

In the scientific field, it is sometime more useful to plot the differences in 
terms of absolute percentage difference/error (APE). Let's plot the APE 
for this scenario.

APE = abs(tiegcm_dens-JB2008_dens)/tiegcm_dens
"""
cs4 = axs[3].contourf(localSolarTimes_JB2008, latitudes_JB2008, abs(tie_JB_grid.T - JB2008_dens_reshaped[:,:,hi2,time_index].squeeze().T)/tie_JB_grid.T)
axs[3].set_title('Density APE at 400 km', fontsize=18)
axs[3].set_ylabel("Latitudes", fontsize=18)
axs[3].tick_params(axis = 'both', which = 'major', labelsize = 16)

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig.colorbar(cs4,ax = axs[3])
cbar.ax.set_ylabel('Density')

axs[3].set_xlabel("Local Solar Time", fontsize=18)

#%%






