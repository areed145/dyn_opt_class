#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code from:
    http://apmonitor.com/pdc/index.php/Main/ArduinoEstimation2
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from gekko import GEKKO

# Import or generate data
filename = 'mimo_test_data.csv'
data = pd.read_csv(filename)

# Create GEKKO Model
m = GEKKO()
m.time = data['tm'].values

# Parameters to Estimate
U = m.FV(value=10,lb=1,ub=20)
alpha1 = m.FV(value=0.01,lb=0.003,ub=0.03)  # W / % heater
alpha2 = m.FV(value=0.005,lb=0.002,ub=0.02) # W / % heater

tau1 = m.FV(value=1,lb=0,ub=100)
tau2 = m.FV(value=1,lb=0,ub=200)
As = m.FV(value=2.0e-2,lb=0.0,ub=1.0)    # Area between heaters in m^2

# STATUS=1 allows solver to adjust parameter
U.STATUS = 1  
alpha1.STATUS = 1 
alpha2.STATUS = 1 

tau1.STATUS = 1  
tau2.STATUS = 1  
As.STATUS = 1  

# Measured inputs
Q1 = m.MV(value=data['Q1'].values)
Q2 = m.MV(value=data['Q2'].values)

# State variables
Tmea1 = m.CV(value=data['T1'].values)
Tmea1.FSTATUS = 1    # minimize fstatus * (meas-pred)^2
#TC1 = m.CV(value=data['T1'].values)
#TC1.FSTATUS = 1    # minimize fstatus * (meas-pred)^2
Tmea2 = m.CV(value=data['T2'].values)
Tmea2.FSTATUS = 1    # minimize fstatus * (meas-pred)^2
#TC2 = m.CV(value=data['T2'].values)
#TC2.FSTATUS = 1    # minimize fstatus * (meas-pred)^2

TC1 = m.SV(value=Tmea1[0])
TC2 = m.SV(value=Tmea2[0])

Ta = m.Param(value=22.5+273.15)     # K
mass = m.Param(value=4.0/1000.0)    # kg
Cp = m.Param(value=0.5*1000.0)      # J/kg-K    
A = m.Param(value=10.0/100.0**2)    # Area not between heaters in m^2
#As = m.Param(value=2.0/100.0**2)    # Area between heaters in m^2
eps = m.Param(value=0.9)            # Emissivity
sigma = m.Const(5.67e-8)            # Stefan-Boltzmann

# Heater temperatures in Kelvin
T1 = m.Intermediate(TC1+273.15)
T2 = m.Intermediate(TC2+273.15)

# Heat transfer between two heaters
Q_C12 = m.Intermediate(U*As*(T2-T1)) # Convective
Q_R12 = m.Intermediate(eps*sigma*As*(T2**4-T1**4)) # Radiative

# Semi-fundamental correlations (energy balances)
m.Equation(TC1.dt() == (1.0/(mass*Cp))*(U*A*(Ta-T1) \
                    + eps * sigma * A * (Ta**4 - T1**4) \
                    + Q_C12 + Q_R12 \
                    + alpha1*Q1))

m.Equation(TC2.dt() == (1.0/(mass*Cp))*(U*A*(Ta-T2) \
                    + eps * sigma * A * (Ta**4 - T2**4) \
                    - Q_C12 - Q_R12 \
                    + alpha2*Q2))

m.Equation(tau1*Tmea1.dt() == 1.0 * (TC1 - Tmea1))
m.Equation(tau2*Tmea2.dt() == 1.0 * (TC2 - Tmea2))

# Options
m.options.IMODE   = 5 # MHE
m.options.EV_TYPE = 2 # Objective type
m.options.NODES   = 2 # Collocation nodes
m.options.SOLVER  = 3 # IPOPT

# Solve
m.solve(disp=True)

# Parameter values
print('U     : ' + str(U.value[0]))
print('alpha1: ' + str(alpha1.value[0]))
print('alpha2: ' + str(alpha2.value[0]))
print('tau1  : ' + str(tau1.value[0]))
print('tau2  : ' + str(tau2.value[0]))
print('As    : ' + str(As.value[0]))

# Create plot
plt.figure(8)
ax=plt.subplot(2,1,1)
ax.grid()
plt.plot(data['tm'],data['T1'],'r-',label=r'$T_1$ measured')
plt.plot(m.time,Tmea1.value,color='red',linestyle='--',\
         linewidth=2,label=r'$T_1$ pred meas')
plt.plot(m.time,TC1.value,color='red',linestyle='--',\
         linewidth=1,label=r'$T_1$ pred element')
plt.plot(data['tm'],data['T2'],'b-',label=r'$T_2$ measured')
plt.plot(m.time,Tmea2.value,color='blue',linestyle='--',\
         linewidth=2,label=r'$T_2$ pred meas')
plt.plot(m.time,TC2.value,color='blue',linestyle='--',\
         linewidth=1,label=r'$T_2$ pred element')
plt.ylabel('Temperature (degC)')
plt.legend(loc=2)
ax=plt.subplot(2,1,2)
ax.grid()
plt.plot(data['tm'],data['Q1'],'r-',\
         linewidth=3,label=r'$Q_1$')
plt.plot(data['tm'],data['Q2'],'b:',\
         linewidth=3,label=r'$Q_2$')
plt.ylabel('Heaters')
plt.xlabel('Time (sec)')
plt.legend(loc='best')
plt.show()
