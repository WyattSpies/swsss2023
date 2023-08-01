#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

__author__ = 'wyatt spies'
__email__ = 'wyatt.spies@colorado.edu'

# ----------------------------------------------------------------------
# Take first derivative of a function
# ----------------------------------------------------------------------

def first_derivative(f, x, f_L, f_R):

    """ Function that takes the first derivative

    Parameters
    ----------
    f - values of a function that is dependent on x -> f(x)
    x - the location of the point at which f(x) is evaluated
    
    f(x) = e^-3x * sin(2/x^2)

    Returns
    -------
    dfdx - the first derivative of f(x)
    
    Analytical : f'(x) = (-3e^-3x)*sin(2/x^2) - (4x^-3)*(e^-3x)*cos(2/x^2)
    
    Numerical : df/dx = [e^-3(x+dx) * sin(2/(x+dx)^2) - e^-3(x-dx) * sin(2/(x-dx)^2)]/[2dx]

    Notes
    -----
    take the first derivative of f(x) here

    """

    nPts = len(f)
    
    dfdx = np.zeros(nPts)
    
    # dx = 0.1

    # do calculation here - need 3 statements:
    #  1. left boundary ( dfdx(0) = ...)
    dfdx[0] = (f[(1)] - f_L)/(2*dx)
    
    #  2. central region (using spans, like dfdx(1:nPts-2) = ...)
    dfdx[1:-1] = [(f[(a + 1)] - f[(a - 1)])/(2*dx) for a in range(1, nPts-1)]
    #print(dfdx)
    
    #  3. right boundary ( dfdx(nPts-1) = ... )
    dfdx[-1] = (f_R - f[-2])/(2*dx)
    
    return dfdx

# ----------------------------------------------------------------------
# Take second derivative of a function
# ----------------------------------------------------------------------

def second_derivative(f, x, f_L, f_R):

    """ Function that takes the second derivative

    Parameters
    ----------
    f - values of a function that is dependent on x -> f(x)
    x - the location of the point at which f(x) is evaluated

    Returns
    -------
    d2fdx2 - the second derivative of f(x)
    
    analytical : f"(x) = (9e^-3x)*sin(2/x^2) + (36x^-3)*(e^-3x)*cos(2/x^2) + (4x^-4)*(e^-3x)*cos(2/x^2)
    
    numerical : d2f/dx2 = [(e^-3(x+dx)*sin(2/(x+dx^2)) + (e^-3(x-df)*sin(2/(x-df)^2) - 2(e^-3x*sin(2/x^2))]/dx^2

    Notes
    -----
    take the second derivative of f(x) here

    """

    nPts = len(f)
    
    d2fdx2 = np.zeros(nPts)

    # do calculation here - need 3 statements:
    #  1. left boundary ( dfdx(0) = ...)
    d2fdx2[0] = (f[1] + f_L - 2*f[0])/(dx**2)
    
    #  2. central region (using spans, like dfdx(1:nPts-2) = ...)
    d2fdx2[1:-1] = [(f[(a + 1)] + f[(a - 1)] - 2*f[a])/(dx**2) for a in range(1, nPts-1)]
    print(d2fdx2)
    
    #  3. right boundary ( dfdx(nPts-1) = ... )
    d2fdx2[-1] = (f_R + f[-2] - 2*f[-1])/(dx**2)
    
    return d2fdx2






# ----------------------------------------------------------------------
# Get the analytic solution to f(x), dfdx(x) and d2fdx2(x)
# ----------------------------------------------------------------------

def analytic(x):

    """ Function that gets analytic solutions

    Parameters
    ----------
    x - the location of the point at which f(x) is evaluated

    Returns
    -------
    f - the function evaluated at x
    dfdx - the first derivative of f(x)
    d2fdx2 - the second derivative of f(x)

    Notes
    -----
    These are analytic solutions!

    """
    #x = np.arange(1,3,0.1)
    f = np.exp(-3*x)*np.sin(2/x**2)
    f_L = np.exp(-3*(x[0]-dx))*np.sin(2/(x[0]-dx)**2)
    f_R = np.exp(-3*(x[-1]+dx))*np.sin(2/(x[-1]+dx)**2)
    
    dfdx = (-3*np.exp(-3*x))*np.sin(2/x**2) - (4*x**(-3))*(np.exp(-3*x))*np.cos(2/x**2)
    d2fdx2 = 9*np.exp(-3*x)*np.sin(2/x**2) + 36*x**(-3)*np.exp(-3*x)*np.cos(2/x**2) + 4*x**(-4)*np.exp(-3*x)*np.cos(2/x**2)

    return f, dfdx, d2fdx2, f_L, f_R


# ----------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------

if __name__ == "__main__":

    # get figures:
    fig = plt.figure(figsize = (10,10))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    # define dx:
    dx = 0.1
    
    # arange doesn't include last point, so add explicitely:
    x = np.arange(1, 3 + dx, dx)

    # get analytic solutions:
    f, a_dfdx, a_d2fdx2, f_L, f_R = analytic(x)

    # get numeric first derivative:
    n_dfdx = first_derivative(f, x, f_L, f_R)

    # get numeric first derivative:
    n_d2fdx2 = second_derivative(f, x, f_L, f_R)

    # plot:
    ax1.plot(x, f)

    # plot first derivatives:
    error = np.sum(np.abs(n_dfdx - a_dfdx)) / len(n_dfdx)
    sError = ' (Err: %5.1f)' % error
    ax2.plot(x, a_dfdx, color = 'black', label = 'Analytic')
    ax2.plot(x, n_dfdx, color = 'red', label = 'Numeric'+ sError)
    ax2.scatter(x, n_dfdx, color = 'red')
    ax2.legend()

    # plot second derivatives:
    ax3.plot(x, a_d2fdx2, color = 'black', label = 'Analytic')
    ax3.plot(x, n_d2fdx2, color = 'red', label = 'Numeric'+ sError)
    ax3.scatter(x, n_d2fdx2, color = 'red')
    ax3.legend()
    
    plotfile = 'plot.png'
    print('writing : ',plotfile)    
    #fig.savefig(plotfile)
    #plt.close()
    plt.show()
