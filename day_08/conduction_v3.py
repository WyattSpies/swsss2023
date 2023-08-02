#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from tridiagonal import solve_tridiagonal

from HydroStatSolve import get_H, get_n

# ----------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------

if __name__ == "__main__":
   
    #variables
    dx = 4.0

    # set x with 1 ghost cell on both sides:
    alt = 100 + np.arange(-dx, 400 + 2*dx, dx)
    Temp = []
    
    n_N2 = []
    n_O2 = []
    n_O = []
    
    nPts = len(alt)

    # set default coefficients for the solver:
    a = np.zeros(nPts) + 1
    b = np.zeros(nPts) - 2
    c = np.zeros(nPts) + 1
    d = np.zeros(nPts)
    
    #make Q
    Q_back = np.zeros(nPts)
    Q_EUV = np.zeros(nPts)
    
    AmpDi = 10
    AmpSd = 5
    PhaseDi = np.pi/2
    PhaseSd = 3*np.pi/2
    
    #make time
    ndays = 3
    dt = 1 #hours
    times = np.arange(0, ndays*24, dt)
    lon = -105.2705
    
    #sun intensity is lin rel to f10.7 
    f107 = 100 + 50/(24*365) + 25*np.sin(times*2*np.pi/(27*24))
    
    for hour in times:
        ut = hour % 24
        
        #Make solar dependent term
        loc_time = lon/15 + ut     
        fac = -np.cos(loc_time/24 * 2*np.pi)
        
        #temp
        T_lower = 200.0 + AmpDi*np.sin(loc_time*2*np.pi/24 + PhaseDi) + AmpSd*np.sin(loc_time*2*2*np.pi/24 + PhaseSd)
        #T_upper = 400.0
        
        if fac < 0:
            fac = 0
        
        sun_heat = f107[hour]* 0.4/100
        lam = 80.0
        dz = alt[1]-alt[0]
        dz2 = dz*dz
        
        # Add a source term:
        Q_back[(alt>200)&(alt<400)] = 0.4
        Q_EUV[(alt>200)&(alt<400)] = sun_heat * fac
        d = -(Q_back+Q_EUV)*dz2/lam 
        
        #Matrix boundary conditions (bottom - fixed):
        a[0] = 0
        b[0] = 1
        c[0] = 0
        d[0] = T_lower
    
        # top - fixed:
        a[-1] = 1
        b[-1] = -1
        c[-1] = 0
        d[-1] = 0    #T_upper
    
        # solve for Temperature:
        T = solve_tridiagonal(a, b, c, d)
        Temp.append(T)
        
        #Append densities
        n_0_N2 = 1e19
        n_0_O2 = 0.3e19
        n_0_O = 1e18
        
        H_N2 = get_H('N2', alt, T, nPts)
        H_O2 = get_H('O2', alt, T, nPts)
        H_O = get_H('O', alt, T, nPts)
        
        #Get_n    
        n_N2.append(get_n(T, H_N2, dx, n_0_N2, nPts))
        n_O2.append(get_n(T, H_O2, dx, n_0_O2, nPts))
        n_O.append(get_n(T, H_O, dx, n_0_O, nPts))
        
        
        #print(np.shape(Temp))
    Temp = np.array(Temp).T
    
    n_N2 = np.array(n_N2).T
    n_O2 = np.array(n_O2).T
    n_O = np.array(n_O).T
    
    print('Temp' + str(np.shape(Temp)))
    print(np.shape(n_N2))
    print(np.shape(n_O2))
    print(np.shape(n_O))
    
    
    # Make fig:
    fig, ax = plt.subplots(2, 2 , figsize = (15,10))
    
    #Plots
    #Plot Temperature
    imT = ax[0,0].contourf(times/24, alt, Temp)   #, cmap='spring'
    
    ax[0,0].set_xlabel('Time (Days)')
    ax[0,0].set_ylabel('Alt. (km)')
    ax[0,0].set_title('Daily Temperatures: Lon = ' + str(lon))
    
    fig.colorbar(imT, label='Temperatures (K)')
    
    #Plot Densities
    #N2
    imN2 = ax[0,1].contourf(times/24, alt, np.log(n_N2))
    
    ax[0,1].set_xlabel('Time (Days)')
    ax[0,1].set_ylabel('Alt. (km)')
    ax[0,1].set_title('Daily N2 Densities: Lon = ' + str(lon))
    
    fig.colorbar(imN2, label='Densities $m^{-3}$')

    #O2
    imO2 = ax[1,0].contourf(times/24, alt, np.log(n_O2))
    ax[1,0].set_xlabel('Time (Days)')
    ax[1,0].set_ylabel('Alt. (km)')
    ax[1,0].set_title('Daily O2 Densities: Lon = ' + str(lon))
    
    fig.colorbar(imO2, label='Densities $m^{-3}$')

    #O
    imO = ax[1,1].contourf(times/24, alt, np.log(n_O))
    ax[1,1].set_xlabel('Time (Days)')
    ax[1,1].set_ylabel('Alt. (km)')
    ax[1,1].set_title('Daily O Densities: Lon = ' + str(lon))
    
    fig.colorbar(imO, label='Densities $m^{-3}$')
    
    
    plotfile = 'conduction_v2.png'
    print('writing : ',plotfile) 
    plt.show(plotfile)
    fig.savefig(plotfile)
    plt.close()
    
    
    
