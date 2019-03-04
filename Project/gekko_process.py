# -*- coding: utf-8 -*-
"""
Dynamic Optimization - Group 17 Project
Managed Pressure Drilling Model in GEKKO


"""

#%%Import packages
import numpy as np
import pandas as pd
from gekko import GEKKO
import matplotlib.pyplot as plt

#%% Non-model parameters
rmt = True  # Solve local or remote

#%% Process Model functions
#%% Reservoir Model
res_data = pd.read_csv('reservoir_v1.csv')
res_MD  = res_data['MD'].values
res_TVD = res_data['TVD'].values
res_ROP = res_data['ROP'].values
res_PF  = res_data['Pf'].values

def reservoir(depth):
    tvd = 0
    rop = 0
    pf = 0
    for i in range(len(res_MD)):
        if (depth >= res_MD[i]):
            if (i+1<len(res_MD)):
                tvd_slope = (res_TVD[i+1]-res_TVD[i]) / (res_MD[i+1]-res_MD[i])
            else:
                tvd_slope = (res_TVD[i]-res_TVD[i-1]) / (res_MD[i]-res_MD[i-1])
        tvd =  res_TVD[i] + (depth-res_MD[i]) * tvd_slope
        rop = res_ROP[i]
        pf = res_PF[i]
    return (tvd, rop, pf)

#%% Drillstring Model
def drillstring():
    # Specify Model
    d = GEKKO(remote=rmt)
    d.time = np.linspace(0, 20, 21)

    # Model Constants
    g = 9.81    # Gravitational constant (kg-m/s^2)
    h = 1951.0  # Height of wellbore (m)

    betaD = d.Const(value = 2.0e+9)   # bulk modulus of the drill string
    betaA = d.Const(value = 1.0e+9)   # bulk modulus of the annulus

    Ma    = d.FV(value = 0.0)      # mass of fluid in the annulus
    Md    = d.FV(value = 0.0)      # mass of fluid in the drill string
    M     = d.FV(value = 4.3e+8)   # Total mass of fluid (Ma + Md)

    rhoA  = d.FV(value = 1580.0)   # density of the fluid in the drill string
    rhoD  = d.FV(value = 1580.0)   # density of the fluid in the annulus
    rhoM  = d.FV(value = 1580.0)   # density of the mud used

    Fa    = d.Const(value = 2.0e+9)   # friction coefficient of the annulus
    Fd    = d.Const(value = 5.0e+9)   # friction coefficient of the drill string
    
    Qres  = d.FV(value = 0.0)      # volume flow rate from the reservoir
    
    # Define Variables
    Qpump = d.FV(value = 2000.0*(1e-3/60))   # volume flow rate through the pump
    Qback = d.FV(value = 800.0*(1e-3/60))    # volume flow rate through the back pressure pump
    Vd    = d.FV(value = 17.0)     # volume of the drill string
    Va    = d.FV(value = 48.0)     # volume of the annulus
    
    Pp    = d.Var(value = 40.0e+5)  # pump pressure inside the drill string
    Pc    = d.Var(value = 10.0e+5)  # choke pressure inside the annulus
    Qbit  = d.Var(value = 2000.0*(1e-3/60))   # volume flow rate through the drill bit
     # volume flow rate through the choke (Qpump+Qback)
    Qchoke = d.Intermediate( Qbit + Qback )
    #dPa = m.Var(value = 0.0)
      
    
    Pbit  = d.Var(value = 0.0)      # Pressure at the drill bit / BHP
    #dPa   = m.Var(value = 0.0)      # differential pressure across the annulus (?)
    
    # Model equations
    # Equation 5.1
    #m.Equation( Pp.dt() == (betaD/Vd) * (Qpump - Qbit - Vd.dt()) )
    d.Equation( Vd*Pp.dt() == (betaD) * (Qpump - Qbit) )
    
    # Equation 5.2
    #m.Equation( Pc.dt() == (betaA/Va) * (Qres + Qbit + Qback - Qchoke - Va.dt()) )
    d.Equation( Va*Pc.dt() == (betaA) * (Qres + Qbit + Qback - Qchoke) )
    
    # Equation 5.3
    #m.Equation( Qbit.dt() == (1/M) * (Pp - Pc - dPa + (rhoD-rhoA)*g*h) )
    #Using the equation by replacing dPa by the pressure balance in drillpipe
    d.Equation( Qbit.dt() == (1/M) * (Pp - Pbit - Fd*(Qbit**2) + rhoD*g*h ) )
    
    # Equation 5.6
    d.Equation( Pbit == Pc + rhoA*g*h + Fa*(Qbit**2) )
    
    # Other equations for ... BALANCE?
    #m.Equation( dPa == Pbit - Pc )
    
    # Solver options
    d.options.IMODE = 4  # Dynamic Simulation
    
    
    # Solve the model
    d.solve()

    return


#%% Specify controller model
tf = 100        # control horizon
npts = 10       # condensation points

c = GEKKO(remote=rmt)
c.time = np.linspace(0, tf, npts+1)  # Controller horizon




#%% Display results
#plt.figure(3)
#plt.subplot(3,1,1)
#plt.plot(m.time, Pbit.VALUE, 'r', label='BHP(Pa)')
#plt.plot(m.time, Pp.VALUE, 'b', label='Pump Pressure(Pa)')
#plt.legend(loc='best')
#plt.subplot(3,1,2)
#plt.plot(m.time, Pc.VALUE, 'g', label ='Choke Pressure(Pa)')
#plt.legend(loc='best')
#plt.subplot(3,1,3)
#plt.plot(m.time, np.array(Qbit.VALUE)*60000, 'k', label ='Qbit (l/min)')
#plt.plot(m.time, np.array(Qchoke.VALUE)*60000, 'r', label ='Qchoke (l/min)')
#plt.legend(loc='best')
#plt.show()

