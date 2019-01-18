#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Orthogonal Collocation on Finite Elements
Homework Exercise #1

Objective: Solve a differential equation with orthogonal collocation on 
finite elements. Create a MATLAB or Python script to simulate and display
the results.

Solve the following differential equation from time 0 to 1 with orthogonal 
collocation on finite elements with 4 nodes for discretization in time.

 5 dx/dt = -x2 + u

Specify the initial condition for x as 0 and the value of the input, u, as 4. 
Compare the solution result with 2-6 time points (nodes). Report the solution 
at the final time for each and comment on how the solution changes with an 
increase in the number of nodes. 


Created on Thu Jan 17 20:14:29 2019
@author: rob
"""

import numpy as np
import matplotlib.pyplot as plt
from gekko import GEKKO

# Initial conditions
x0 = 0.0
tf = 1.0
nodes = [2,3,4,5,6]
ans = np.zeros(5)
i = 0

# Analytic solution
x_a = x0 + 2.0*np.tanh(tf*2.0/5.0)
print('Analytic: x =', x_a)
print()

for node in nodes:

    # Create model
    m = GEKKO(remote=False)
    
    # Timeline
    m.time = [0, tf]
    
    # Parameters
    u = m.Param(value=4)
    
    # Intermediates
    
    # Variables
    x = m.Var(value=x0)
    
    # Equations
    m.Equation(5 * x.dt() == -x**2 + u)
    
    # Global Solution Options
    m.options.IMODE = 4       # Dynamic Simultaneous ODE Solver
    m.options.TIME_SHIFT = 0  # Not really doing dynamics, no previous sol'n
    m.options.NODES = node    # Set number of collocation nodes
    
    # Solve model
    m.solve(disp=False)
    ans[i] = np.abs(x[-1] - x_a)
    i += 1
    
    print('Nodes =', node, ', x =', x[-1], ',error =', x[-1]-x_a, \
          ', percent =', 100*(x[-1]-x_a)/x_a)

plt.plot(nodes,ans,'ro')
plt.ylabel('|Error|')
plt.xlabel('Nodes')
plt.show()
