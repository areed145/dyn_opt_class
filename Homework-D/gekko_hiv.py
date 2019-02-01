# -*- coding: utf-8 -*-
"""
Model Estimation Exercise: HIV
http://apmonitor.com/do/index.php/Main/EstimatorObjective
Created on Wed Jan 30 14:11:57 2019
@author: Rob Hawkins

Problem Description:
    
Initial Conditions
 H = healthy cells = 1,000,000
 I = infected cells = 0
 V = virus = 100
 LV = log virus = 2

Equations
 dH/dt = kr1 - kr2 H - kr3 H V
 dI/dt = kr3 H V - kr4 I
 dV/dt = -kr3 H V - kr5 V + kr6 I
 LV = log10(V)
 
Parameters
 kr1 = new healthy cells
 kr2 = death rate of healthy cells
 kr3 = healthy cells converting to infected cells
 kr4 = death rate of infected cells
 kr5 = death rate of virus
 kr6 = production of virus by infected cells
 
Initial Parameter Values: (from hiv.py model provided)
 kr[0] = 1e5
 kr[1] = 0.1
 kr[2] = 2e-7
 kr[3] = 0.5
 kr[4] = 5
 kr[5] = 100
"""

# Import necessary modules
from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt

# Load data into Array
data = np.loadtxt('hiv_data.csv', delimiter=',', skiprows=1)
LV_meas = data[:,1]  # log10(V) is in second column

# Create a gekko model
m = GEKKO(remote=False)

# Set timeline to match data in CSV file
m.time = np.linspace(0,15,np.size(data,0))

# Define Parameters
kr_init = [1e4, 0.1, 2e-7, 0.5, 5.0, 100.0]
lkr = [m.FV(value=np.log10(kr_init[i])) for i in range(6)]

# Define Variables
kr =[m.Var(value=kr_init[i]) for i in range(6)]  # Raw Parameters
H = m.Var(value=1000000)                         # healthy cells
I = m.Var(value=0)                               # infected cells
V = m.Var(value=100)                              # virus

# Define measurements
LV = m.CV(value=np.log10(100))                    # log virus (raw data)

# Define Equations
# dH/dt = kr1 - kr2 H - kr3 H V
# dI/dt = kr3 H V - kr4 I
# dV/dt = -kr3 H V - kr5 V + kr6 I
# LV = log10(V)
m.Equations([10**lkr[i] == kr[i] for i in range(6)])
m.Equations([H.dt() == kr[0] - kr[1]*H - kr[2]*H*V,\
            I.dt() == kr[2]*H*V - kr[3]*I,\
            V.dt() == -kr[2]*H*V -kr[4]*V + kr[5]*I,\
            LV == m.log10(V)])

# Solver Options
m.options.IMODE = 5         # dynamic estimation
m.options.TIME_SHIFT = 0    # don't timeshift on new solve
m.options.EV_TYPE = 2       # l2 norm
m.options.COLDSTART = 2
m.options.SOLVER = 1        # APOPT
m.options.MAX_ITER = 2000

# Solve the initial problem
m.solve()

# Set parameters to be optimized
for i in range(5):
    lkr[i].STATUS = 1 # Allow optimizer to fit
    lkr[i].DMAX = 2
    lkr[i].LOWER = -20
    lkr[i].UPPER = 20
    
# Connect data to LV variable
LV.FSTATUS = 1  # Connect to measurement
LV.STATUS = 1   # Use pred-meas in ObjFcn
LV.value = LV_meas

# Solve the fitting problem
m.solve()

# Display the results
plt.figure(1)
plt.semilogy(m.time,H,'b-')
plt.semilogy(m.time,I,'g:')
plt.semilogy(m.time,V,'r--')
plt.semilogy(data[:,][:,0],np.power(10, LV_meas),'r.')
plt.xlabel('Time (yr)')
plt.ylabel('Cell & Virus Count (log scale)')
plt.legend(['H','I','V'])
plt.show()

for i in range(6):
    print ('kr[',i,'] =', kr[i].value[-1])