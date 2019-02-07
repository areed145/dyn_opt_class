# -*- coding: utf-8 -*-
"""
Dynamic Optimization - Benchmarks
Created on Thu Feb  7 13:57:21 2019
see - http://apmonitor.com/do/index.php/Main/DynamicOptimizationBenchmarks

Example 2

    minimize: x4[tf]
    
    subject to:
        dx1/dt = x2
        dx2/dt = -x3*u + 16*t - 8
        dx3/dt = u
        dx4/dt = x1^2 + x2^2 + 0.0005*(x2+16t-8-0.1x3*u^2)^2
        x(0) = [0,-1,-sqrt(5),0]T
        -4 <= u <= 10
        tf = 1

"""

#%%Import packages
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

#%% Non-model parameters
rmt = False  # Solve local or remote

tf = 1.0     # Final time
npts = 1000   # time steps

#%% Specify model
m = GEKKO(remote=rmt)

m.time = np.linspace(0, tf, npts+1)  # Model Timeline

# Define Variables
x1 = m.Var(value=0.0)
x2 = m.Var(value=-1.0)
x3 = m.Var(value=-np.sqrt(5))
x4 = m.Var(value=0.0)
t = m.Var(value=0.0)

u = m.MV(value = 0.0, lb=-4, ub=10)
u.STATUS = 1
u.DCOST = 0

# Set up objective function to use only last value of x2
lastval = m.Param(value = np.array([0.0]*npts + [1.0]))
m.Obj(lastval*x4)

# Model equations
m.Equations([x1.dt() == x2, \
             x2.dt() == -x3*u + 16*t -8.0, \
             x3.dt() == u, \
             x4.dt() == x1**2 + x2**2 + 0.0005*(x2+16*t-8-0.1*x3*u**2)**2, \
             t.dt() == 1.0 ])

# Solver options
m.options.IMODE = 6  # Dynamic Optimization

# Solve the model
m.solve()

#%% Display the results
plt.figure('Example 2')
plt.subplot(2,1,1)
plt.plot(m.time,x1.value,'b-', label=r'$x_1$')
plt.plot(m.time,x2.value,'b--', label=r'$x_2$')
plt.plot(m.time,x3.value,'b-.', label=r'$x_3$')
plt.plot(m.time,x4.value,'b:', label=r'$x_4$')
plt.xlabel('time')
plt.ylabel('States')
plt.legend(loc='best')

plt.subplot(2,1,2)
plt.plot(m.time,u.value,'g-.', label=r'$u$')
plt.xlabel('time')
plt.ylabel('MV')
plt.legend(loc='best')
plt.show()