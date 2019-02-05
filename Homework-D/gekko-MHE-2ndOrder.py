# -*- coding: utf-8 -*-
"""
2nd Order System Identification with MHE

Based on code example from:
http://apmonitor.com/do/index.php/Main/TCLabD
"""
import numpy as np
import matplotlib.pyplot as plt
from gekko import GEKKO
import pandas as pd

rmt = False
fig = 'tclab_mhe_03.png'
skip = 3

# import data
data = pd.read_csv('mimo_test_data_6.csv')

raw_tm = data['tm'].values
raw_Q1m = data['H1'].values  # heater 1
raw_Q2m = data['H2'].values  # heater 2
raw_T1m = data['T1'].values
raw_T2m = data['T2'].values
raw_n = len(raw_tm)

tm = np.zeros(int(raw_n/skip)+1)
Q1m = np.zeros(int(raw_n/skip)+1)
Q2m = np.zeros(int(raw_n/skip)+1)
T1m = np.zeros(int(raw_n/skip)+1)
T2m = np.zeros(int(raw_n/skip)+1)
n = len(tm)

# down sample the data
for i in range(n):
    tm[i] = raw_tm[skip*i]
    Q1m[i] = raw_Q1m[skip*i]
    Q2m[i] = raw_Q2m[skip*i]
    T1m[i] = raw_T1m[skip*i]
    T2m[i] = raw_T2m[skip*i]

# Store MHE values for plots
Tmhe1 = T1m[0] * np.ones(n)
Tmhe2 = T2m[0] * np.ones(n)
K1s = 0.5 * np.ones(n)
K2s = 0.3 * np.ones(n)
K3s = 0.005 * np.ones(n)
tau12s = 150.0 * np.ones(n)
tau3s = 5.0 * np.ones(n)

#########################################################
# Initialize Model as Estimator
#########################################################
m = GEKKO(name='tclab-mhe',remote=rmt)

# 120 second time horizon, 40 steps
m.time = np.linspace(0,120,int(120/skip+1))

# Parameters to Estimate
K1 = m.FV(value=0.5)
K1.STATUS = 0
K1.FSTATUS = 0
K1.DMAX = 0.1
K1.LOWER = 0.1
K1.UPPER = 10.0

K2 = m.FV(value=0.3)
K2.STATUS = 0
K2.FSTATUS = 0
K2.DMAX = 0.1
K2.LOWER = 0.1
K2.UPPER = 10.0

K3 = m.FV(value=0.2)
K3.STATUS = 0
K3.FSTATUS = 0
K3.DMAX = 0.01
K3.LOWER = 0.1
K3.UPPER = 10.0

tau12 = m.FV(value=150)
tau12.STATUS = 0
tau12.FSTATUS = 0
tau12.DMAX = 5.0
tau12.LOWER = 50.0
tau12.UPPER = 1000

tau3 = m.FV(value=15)
tau3.STATUS = 0
tau3.FSTATUS = 0
tau3.DMAX = 1
tau3.LOWER = 10
tau3.UPPER = 200

# Measured inputs
Q1 = m.MV(value=0)
Q1.FSTATUS = 1 # measured

Q2 = m.MV(value=0)
Q2.FSTATUS = 1 # measured

# State variables
TH1 = m.SV(value=T1m[0])
TH2 = m.SV(value=T2m[0])

# Measurements for model alignment
TC1 = m.CV(value=T1m[0])
TC1.STATUS = 1     # minimize error
TC1.FSTATUS = 1    # receive measurement
TC1.MEAS_GAP = 0.1 # measurement deadband gap

TC2 = m.CV(value=T2m[0])
TC2.STATUS = 1     # minimize error
TC2.FSTATUS = 1    # receive measurement
TC2.MEAS_GAP = 0.1 # measurement deadband gap

Ta = m.Param(value=26.0) # degC

# Heat transfer between two heaters
DT = m.Intermediate(TH2-TH1)

# Empirical correlations
m.Equation(tau12 * TH1.dt() + (TH1-Ta) == K1*Q1 + K3*DT)
m.Equation(tau12 * TH2.dt() + (TH2-Ta) == K2*Q2 - K3*DT)
m.Equation(tau3 * TC1.dt()  + TC1 == TH1)
m.Equation(tau3 * TC2.dt()  + TC2 == TH2)

# Global Options
m.options.IMODE   = 5 # MHE
m.options.EV_TYPE = 1 # Objective type
m.options.NODES   = 3 # Collocation nodes
m.options.SOLVER  = 3 # IPOPT
m.options.COLDSTART = 1 # COLDSTART on first cycle
m.options.TIME_SHIFT = skip
##################################################################
# Create plot
plt.figure(figsize=(10,7))
plt.ion()
plt.show()

# Main Loop
start_time = tm[0]
prev_time = start_time

try:
    for i in range(1, n):
        # Record time and change in time
        #dt = tm[i] - tm[i-1]

        # Insert measurements
        TC1.MEAS = T1m[i]
        TC2.MEAS = T2m[i]
        Q1.MEAS = Q1m[i-1]
        Q2.MEAS = Q2m[i-1]

        # Start estimating U after ~30 sec
        if tm[i] > 30:
            K1.STATUS = 1
            K2.STATUS = 1
            K3.STATUS = 1
            tau12.STATUS = 1
            tau3.STATUS = 1

        # Predict Parameters and Temperatures with MHE
        m.solve(disp=False) 

        if m.options.APPSTATUS == 1:
            # Retrieve new values
            Tmhe1[i]  = TC1.MODEL
            Tmhe2[i]  = TC2.MODEL
            K1s[i]    = K1.NEWVAL
            K2s[i]    = K2.NEWVAL
            K3s[i]    = K3.NEWVAL
            tau12s[i] = tau12.NEWVAL
            tau3s[i]  = tau3.NEWVAL
        else:
            # Solution failed, copy prior solution
            Tmhe1[i]  = Tmhe1[i-1]
            Tmhe2[i]  = Tmhe1[i-1]
            K1s[i]    = K1s[i-1]   
            K2s[i]    = K2s[i-1]   
            K3s[i]    = K3s[i-1]   
            tau12s[i] = tau12s[i-1]
            tau3s[i]  = tau3s[i-1] 

        # Plot
        plt.clf()
        ax=plt.subplot(4,1,1)
        ax.grid()
        plt.plot(tm[0:i],T1m[0:i],'r.',label=r'$T_1$ measured')
        plt.plot(tm[0:i],Tmhe1[0:i],'k-',label=r'$T_1$ MHE')
        plt.plot(tm[0:i],T2m[0:i],'b.',label=r'$T_2$ measured')
        plt.plot(tm[0:i],Tmhe2[0:i],'k--',label=r'$T_2$ MHE')
        plt.ylabel('Temperature (degC)')
        plt.legend(loc=2)
        ax=plt.subplot(4,1,2)
        ax.grid()
        plt.plot(tm[0:i],K1s[0:i],'k-',label='K1')
        plt.plot(tm[0:i],K2s[0:i],'g:',label='K2')        
        plt.plot(tm[0:i],K3s[0:i],'r--',label='K3')
        plt.ylabel('Gains')
        plt.legend(loc='best')
        ax=plt.subplot(4,1,3)
        ax.grid()
        plt.plot(tm[0:i],tau12s[0:i],'b-',label=r'$\tau_{12}$')
        plt.plot(tm[0:i],tau3s[0:i]*10,'r--',label=r'$\tau_3$ x 10')
        plt.ylabel('Time constant')
        plt.legend(loc='best')
        ax=plt.subplot(4,1,4)
        ax.grid()
        plt.plot(tm[0:i],Q1m[0:i],'r-',label=r'$Q_1$')
        plt.plot(tm[0:i],Q2m[0:i],'b:',label=r'$Q_2$')
        plt.ylabel('Heaters')
        plt.xlabel('Time (sec)')
        plt.legend(loc='best')
        plt.draw()
        plt.pause(0.05)

    # Save figure
    plt.savefig(fig)

# Allow user to end loop with Ctrl-C           
except KeyboardInterrupt:
    print('Shutting down')
    plt.savefig(fig)

