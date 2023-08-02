#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 16:41:10 2023

@author: holtorf
"""

import numpy as np
import matplotlib.pyplot as plt

# Use Euler's method with different stepsizes to solve the IVP:
# dx/dt = -2*x, with x(0) = 3 over the time-horizon [0,2]


def f(x):
    return -2*x

#time step
h = 0.2

t_0 = 0
t_f = 2

x_0 = 3

nPts = int(np.ceil(t_f/h))

t = np.zeros(nPts+1)
x = np.zeros(nPts+1)

x[0] = x_0
t[0] = t_0

for i in range(nPts):
    x[i+1] = x[i] + h*f(x[i])
    t[i+1] = t[i] + h

fig, ax = plt.subplots()
ax.scatter(t,x,linewidth=1, color='black', label='Explicit Euler')

t_range = np.linspace(0, t_f, 1000)

ax.plot(t_range, 3*np.exp(-2*t_range), linewidth=1, color='red', label='Analytical')
ax.legend()

# Compare the numerical approximation of the IVP solution to its analytical
# solution by plotting both solutions in the same figure. 

