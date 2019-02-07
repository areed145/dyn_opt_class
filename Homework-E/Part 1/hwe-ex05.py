# -*- coding: utf-8 -*-
"""
Dynamic Optimization - Benchmarks
Created on Thu Feb  7 13:57:21 2019
see - http://apmonitor.com/do/index.php/Main/DynamicOptimizationBenchmarks

Example 5

    maximize: 1 - x1(tf) - x2(tf)
    
    subject to:
        dx1/dt = u*(10*x2-x1)
        dx2/dt = -u*(10*x2-x1) - (1-u)*x2
        x(0) = [1,0]T
        0 <= u <= 1
        tf = 12

"""

#%%Import packages
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

#%% Non-model parameters
rmt = False  # Solve local or remote

tf = 12.0     # Final time
npts = 1000   # time steps

#%% Specify model
m = GEKKO(remote=rmt)

m.time = np.linspace(0, tf, npts+1)  # Model Timeline

# Define Variables
x1 = m.Var(value=1.0)
x2 = m.Var(value=0.0)

u = m.MV(value = 0.0, lb=0, ub=1)
u.STATUS = 1
u.DCOST = 0

# Set up objective function to use only last value of x2
lastval = m.Param(value = np.array([0.0]*npts + [1.0]))
m.Obj(-lastval*(1.0-x1-x2))

# Model equations
m.Equations([x1.dt() == u * (10*x2-x1), \
             x2.dt() == -u*(10*x2-x1) - (1-u)*x2 ])

# Solver options
m.options.IMODE = 6  # Dynamic Optimization

# Solve the model
m.solve()

#%% Display the results
plt.figure('Example 5')
plt.subplot(2,1,1)
plt.plot(m.time,x1.value,'b-', label=r'$x_1$')
plt.plot(m.time,x2.value,'b--', label=r'$x_2$')
plt.xlabel('time')
plt.ylabel('States')
plt.legend(loc='best')

plt.subplot(2,1,2)
plt.plot(m.time,u.value,'g-.', label=r'$u$')
plt.xlabel('time')
plt.ylabel('MV')
plt.legend(loc='best')
plt.show()