#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 17:13:36 2019

@author: areed145
"""

import numpy as np
import time
import matplotlib.pyplot as plt
import random
# get gekko package with:
#   pip install gekko
from gekko import GEKKO
import pandas as pd

# import data
try:
    # read data file if available
    data = pd.read_csv('tclab_dyn_data2.csv')
except:
    # retrieve data file from Internet source
    url = 'http://apmonitor.com/do/uploads/Main/tclab_dyn_data2.txt'
    data = pd.read_csv(url) 
    data.to_csv('tclab_dyn_data2.csv')

tm = data['Time'].values
Q1s = data['H1'].values  # heater 1
Q2s = data['H2'].values  # heater 2
T1s = data['T1'].values
T2s = data['T2'].values

#########################################################
# Initialize Model as Estimator
#########################################################
m = GEKKO(remote=True)

m.time = tm

# Parameters to Estimate
K1 = m.FV(value=0.5,lb=0.1,ub=1.0)
K2 = m.FV(value=0.3,lb=0.1,ub=1.0)
K3 = m.FV(value=0.1,lb=0.0001,ub=1.0)
tau12 = m.FV(value=150,lb=50,ub=250)
tau3 = m.FV(value=150,lb=50,ub=250)

K1.STATUS = 1
K2.STATUS = 1
K3.STATUS = 1
tau12.STATUS = 1
tau3.STATUS = 1

# Measured inputs
Q1 = m.MV(value=Q1s)
Q2 = m.MV(value=Q2s)

# use measurements
Q1.FSTATUS = 1 # measured
Q2.FSTATUS = 1 # measured

# Ambient temperature
Ta = m.Param(value=19.0) # degC

TH1 = m.SV(value=T1s[0])
TH2 = m.SV(value=T2s[0])

# Measurements for model alignment
TC1 = m.CV(value=T1s)
TC1.STATUS = 1     # minimize error between simulation and measurement
TC1.FSTATUS = 1    # receive measurement
TC1.MEAS_GAP = 0.1 # measurement deadband gap

TC2 = m.CV(value=T2s)
TC2.STATUS = 1     # minimize error between simulation and measurement
TC2.FSTATUS = 1    # receive measurement
TC2.MEAS_GAP = 0.1 # measurement deadband gap

# Heat transfer between two heaters
DT = m.Intermediate(TC2-TC1)

# Empirical correlations
m.Equation(tau12 * TH1.dt() + (TH1-Ta) == K1*Q1 + K3*DT)
m.Equation(tau12 * TH2.dt() + (TH2-Ta) == K2*Q2 - K3*DT)
m.Equation(tau3 * TC1.dt() + TC1 == TH1)
m.Equation(tau3 * TC2.dt() + TC2 == TH2)

# Global Options
m.options.IMODE   = 5 # MHE
m.options.EV_TYPE = 2 # Objective type
m.options.NODES   = 3 # Collocation nodes
m.options.SOLVER  = 3 # IPOPT

# Predict Parameters and Temperatures
m.solve() 

# Create plot
plt.figure(figsize=(10,7))

ax=plt.subplot(2,1,1)
ax.grid()
plt.plot(tm,T1s,'ro',label=r'$T_1$ measured')
plt.plot(tm,TC1.value,'k-',label=r'$T_1$ predicted')
plt.plot(tm,T2s,'bx',label=r'$T_2$ measured')
plt.plot(tm,TC2.value,'k--',label=r'$T_2$ predicted')
plt.ylabel('Temperature (degC)')
plt.legend(loc=2)
ax=plt.subplot(2,1,2)
ax.grid()
plt.plot(tm,Q1s,'r-',label=r'$Q_1$')
plt.plot(tm,Q2s,'b:',label=r'$Q_2$')
plt.ylabel('Heaters')
plt.xlabel('Time (sec)')
plt.legend(loc='best')

# Print optimal values
print('K1: ' + str(K1.newval))
print('K2: ' + str(K2.newval))
print('K3: ' + str(K3.newval))
print('tau12: ' + str(tau12.newval))

# Save and show figure
plt.savefig('tclab_2nd_order.png')
plt.show()