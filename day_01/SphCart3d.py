# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 14:59:01 2023

3d plot script for sperical coords
"""

__author__ = 'Wyatt Spies'
__email__ = 'wyattspies@colorado.edu'


import numpy as np
import matplotlib.pyplot as plt

"""
theta = angle from x-axis  (0,2pi)   (radians)
phi = angle from z-axis   (0,pi)    (radians)
r = distance from origin      (units)
"""

def sph_to_cart (r, theta, phi):
    x = r*np.sin(phi)*np.cos(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(phi)
#    print(x, y, z)
    return x,y,z


#sph_to_cart(1, np.pi/2, np.pi/2)


# Make function test scenarios
test1 = sph_to_cart(1,0,0)
test2 = sph_to_cart(1,np.pi,np.pi)
test3 = sph_to_cart(1,2*np.pi,2*np.pi)
test4 = sph_to_cart(1,-np.pi/2,-2*np.pi)
test5 = sph_to_cart(1,-2*np.pi,np.pi/2)

#test function at different points
assert np.allclose(test1, (0,0,1), rtol = 1e-5), "uh oh, look who just blew in from stupid town"
assert np.allclose(test2, (0,0,-1), rtol = 1e-5), "uh oh, look who just blew in from stupid town"
assert np.allclose(test3, (0,0,1), rtol = 1e-5), "uh oh, look who just blew in from stupid town"
assert np.allclose(test4, (0,0,1), rtol = 1e-5), "uh oh, look who just blew in from stupid town"
assert np.allclose(test5, (1,0,0), rtol = 1e-5), "uh oh, look who just blew in from stupid town"



#plot a function from sph to cart
fig = plt.figure()  #better control
ax = plt.axes(projection='3d')   # make 3d axes

r = np.linspace(0,1)
theta = np.linspace(0, 2*np.pi)
phi = np.linspace(0, 2*np.pi)

x, y, z = sph_to_cart(r,theta,phi)
ax.plot(x,y,z)



