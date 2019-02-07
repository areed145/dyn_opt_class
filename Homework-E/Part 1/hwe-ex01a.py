# -*- coding: utf-8 -*-
"""
Dynamic Optimization - Benchmarks
Created on Thu Feb  7 13:57:21 2019
see - http://apmonitor.com/do/index.php/Main/DynamicOptimizationBenchmarks

    minimize: x2[tf]
    
    subject to:
        dx1/dt = u
        dx2/dt = x1^2 + u^2
        x(0) = [1,0]T
        tf = 1

"""

#%%Import packages
import numpy as np
from random import random
from gekko import GEKKO
import matplotlib.pyplot as plt

#%% Non-model parameters
rmt = False  # Solve local or remote

tf = 1.0     # Final time
npts = 100   # time steps

#%% Specify model
m = GEKKO(remote=rmt)

m.time = np.linspace(0, tf, npts+1)  # Model Timeline

# Define Variables
x1 = m.Var(value=1.0)
x2 = m.Var(value=0.0)
u = m.Var(value = 0.0)

# Set up objective function to use only last value of x2
lastval = m.Param(value = np.array([0.0]*npts + [1.0]))
m.Obj(lastval*x2)

# Model equations
m.Equations([x1.dt() == u, \
             x2.dt() == x1**2 + u**2 ])

# Solver options
m.options.IMODE = 6  # Dynamic Optimization

# Solve the model
m.solve()

#%% Display the results
plt.figure('Example 1a')
plt.plot(m.time,x1.value,'b-', label=r'$x_1$')
plt.plot(m.time,x2.value,'b--', label=r'$x_2$')
plt.plot(m.time,u.value,'g-.', label=r'$u$')
plt.xlabel('time')
plt.ylabel('values')
plt.legend(loc='best')
plt.show()