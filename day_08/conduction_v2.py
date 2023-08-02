#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from tridiagonal import solve_tridiagonal

# ----------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------

if __name__ == "__main__":
   
    
    #____
    dx = 0.25

    # set x with 1 ghost cell on both sides:
    x = np.arange(-dx, 10 + 2 * dx, dx)

    T_lower = 200.0
    #T_upper = 400.0
    Temp = []

    nPts = len(x)

    # set default coefficients for the solver:
    a = np.zeros(nPts) + 1
    b = np.zeros(nPts) - 2
    c = np.zeros(nPts) + 1
    d = np.zeros(nPts)
    
    #make Q
    Q_back = np.zeros(nPts)
    Q_EUV = np.zeros(nPts)
    
    #make time
    ndays = 3
    dt = 1 #hours
    times = np.arange(0, ndays*24, dt)
    lon = -105.2705
    for hour in times:
        ut = hour % 24
        
        #Make solar dependent term
        loc_time = lon/15 + ut     
        fac = -np.cos(loc_time/24 * 2*np.pi)
        
        if fac < 0:
            fac = 0
        
        sun_heat = 100
        lam = 10.0
        dz = x[1]-x[0]
        dz2 = dz*dz
        
        # Add a source term:
        Q_back[(x>3)&(x<7)] = 100
        Q_EUV[(x>3)&(x<7)] = sun_heat * fac
        d = -(Q_back+Q_EUV)*dz2/lam 
        
        # boundary conditions (bottom - fixed):
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
        #print(np.shape(Temp))
    
    Temp = np.array(Temp).T
    print('shape = ' + str(np.shape(Temp)))
    alt = 100 + 40*x
    
    # Make fig:
    fig, ax = plt.subplots()
    im = ax.contourf(times/24, alt, Temp, cmap='spring')
    
    ax.set_xlabel('Time (Days)')
    ax.set_ylabel('Alt. (km)')
    ax.set_title('Daily Temperatures: Lon = ' + str(lon))
    fig.colorbar(im, label='Temperature (K)')
    #ax.legend()
    
    plotfile = 'conduction_v2.png'
    print('writing : ',plotfile) 
    plt.show(plotfile)
    fig.savefig(plotfile)
    plt.close()
    
    
    
