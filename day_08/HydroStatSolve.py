# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 15:17:13 2023

@author: wyatt
"""
__author__ = 'wyatt spies'
__email__ = 'wyatt.spies@colorado.edu'


import numpy as np


def get_H(element, alt, T, nPts):
    #Boltzmann constant
    k = 1.38e-23
    #Elements Dictionary for mass
    amu = {'N2': 28,
       'O2': 32,
       'O': 16,
       }
    m = amu[element]*1.67e-27
    
    #Radius array
    rad = (6370 + alt)*1000
    
    #g array
    g = (3.99e14)/(rad)**2
    
    #H array
    H = (k/m)*T/g
    
    return H
    
def get_n(T, H, delt_z, n_0, nPts):
    #n array
    n = np.zeros(len(T))
    n[0] = n_0
    
    #append nth terms
    for i in range(0, nPts-1):
        n[i+1] = (T[i]/T[i+1])*n[i]*np.exp(-delt_z*1000/H[i])
        
    return n












