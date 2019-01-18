#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 16:52:12 2019
Orthogonal Collocation on Finite Elements
Homework Exercise #2

Objective: Compare orthogonal collocation on finite elements with 3 nodes 
with a numerical integrator (e.g. ODE15s in MATLAB or ODEINT in Python). 
Calculate the error at each of the solution points for the equation 
(same as for Exercise 1) 

@author: rob
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from gekko import GEKKO

# Initial conditions
x0 = 0.0
tf = 1.0
ans = np.zeros(5)
i = 0
t = np.linspace(0,tf,20)
x_a = np.zeros(20)

plt.figure(2, clear=True)

# Analytic solution
for i in range(20):
    x_a[i] = x0 + 2.0*np.tanh(t[i]*2.0/5.0)
plt.plot(t,x_a,'r--',label='Analytic')
    
# solve with ODEINT 
def model(x,t):
    return (-x**2 + 4.0)/5.0

x_ode = odeint(model,x0,t)
plt.plot(t,x_ode,'g.',label='ODEINT')

# Solve with 3-node collocation in GEKKO
# Create model
m = GEKKO(remote=False)
u = m.Param(value=4)
x1 = m.Var()
x2 = m.Var()
dx1 = m.Var()
dx2 = m.Var()
# Equations (use N array)
m.Equations([5 * dx1 == -x1**2 + u, \
             5 * dx2 == -x2**2 + u, \
             0.75*dx1-0.25*dx2==x1-x0, \
             1.00*dx1+0.0*dx2==x2-x0])
m.options.IMODE = 1       # Non-Dynamic Simulation
m.solve(disp=False)
t = np.array([0.0,0.5,1.0])
x_g =  np.array([x0,x1[-1],x2[-1]])
plt.plot(t,x_g,'ro',label='3-node Collocation')

plt.legend(loc='best')
plt.ylabel('x(t)')
plt.xlabel('time')

plt.show()