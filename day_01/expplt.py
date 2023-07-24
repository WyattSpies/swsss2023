# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:50:58 2023

@author: wyatt
"""
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,1)
plt.plot(x, np.exp(x))
plt.plot(x,x)
plt.xlabel('$0 \leq x < 1$')       #r stands for rich text
plt.ylabel('$e^x$')
plt.title('Exponential Function')
plt.show