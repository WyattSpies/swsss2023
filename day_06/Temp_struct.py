# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 15:17:13 2023

@author: wyatt
"""
__author__ = 'wyatt spies'
__email__ = 'wyatt.spies@colorado.edu'


import numpy as np
import matplotlib.pyplot as plt


"""
Conditions:
n_0 = 1E19
alt_0 = 100 km
alt_n = 500 km
nPts = 100
r = 6370 km
g = 3.99E14/((r+alt)*1000)**2
T_0 = 200
T_n = 1000  linear
m = 28*1.67E-27 kg
k = 1.38E-23 units


Steps:
1. Alt array
2. radius array
3. g array
4. T array
5. H array
6. n array
"""



def get_H(element, alt_0, alt_n, T_0, T_n, nPts):
    
    m = amu[element]*1.67e-27
    
    #1. Alt array (linear)
    delt_z = (alt_n - alt_0)/nPts
    #print(delt_z)
    alt = np.arange(alt_0, alt_n, delt_z)
    print('Alt_0 = ' + str(alt[0]))
    
    #2. Radius array
    rad = (6370 + alt)*1000
    print('Radius_0 = ' + str(rad[0]))
    
    #3. g array
    g = (3.99e14)/(rad)**2
    print('g_0 = ' + str(g[0]))
    
    #4. T array (linear)
    delt_T = (T_n-T_0)/nPts
    T = np.arange(T_0, T_n, delt_T)
    print('T_0 = ' + str(T[0]))
    
    #5. H array
    H = (k/m)*T/g
    print('H_0 = ' + str(H[0]), ', H_100 = ' + str(H[-1]))
    
    return alt, T, H, delt_z
    
def get_n(T, H, delt_z, n_0):
    #6. n array
    n = np.zeros(len(T))
    n[0] = n_0
    #print(n[0])
    
    #append nth terms
    for i in range(0, nPts-1):
        n[i+1] = (T[i]/T[i+1])*n[i]*np.exp(-delt_z*1000/H[i])
        
    print('n0 = ' + str(n[0]), 'n1 = ' + str(n[1]), 'n_100 = ' + str(n[-1]))
        
    return n

#Constants
k = 1.38e-23        #check units
#m = 28*1.67e-27
n_0_N2 = 1e19
n_0_O2 = 0.3e19
n_0_O = 1e18

#def parameters
alt_0 = 100
alt_n = 500
T_0 = 200
T_n = 1000
nPts = 100

amu = {'N2': 28,
       'O2': 32,
       'O': 16,
     }

#Choose Element
#element = '02'

#call get_H to retrieve H array
alt_N2, T_N2, H_N2, delt_z = get_H('N2', alt_0, alt_n, T_0, T_n, nPts)
alt_O2, T_O2, H_O2, delt_z = get_H('O2', alt_0, alt_n, T_0, T_n, nPts)
alt_O, T_O, H_O, delt_z = get_H('O', alt_0, alt_n, T_0, T_n, nPts)

#call get_n to retrieve n array
n_N2 = get_n(T_N2, H_N2, delt_z, n_0_N2)
n_O2 = get_n(T_O2, H_O2, delt_z, n_0_O2)
n_O = get_n(T_O, H_O, delt_z, n_0_O)

#plot n vs alt
plt.plot(n_N2, alt_N2, 'blue')
plt.plot(n_O2, alt_O2, 'red')
plt.plot(n_O, alt_O, 'green')
plt.xlabel('Density $(m^{-3})$')
plt.xscale('log')
plt.ylabel('Altitude $(km)$')
plt.title('1D Atmosphere Density')
plt.show()












