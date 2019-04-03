# -*- coding: utf-8 -*-
"""
Dynamic Optimization - Group 17 Project
Managed Pressure Drilling Model in GEKKO


"""

#%%Import packages
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

#%% Non-model parameters
rmt = True  # Solve local or remote

tf = 20.0     # Final time
npts = 20   # time steps

#%% Specify model
m = GEKKO(remote=rmt)
m.time = np.linspace(0, tf, npts+1)  # Model Timeline

# Model Constants
g = 9.81    # Gravitational constant (kg-m/s^2)
h = 1951.0  # Height of wellbore (m)

betaD = m.Const(value = 2.0e+9)   # bulk modulus of the drill string
betaA = m.Const(value = 1.0e+9)   # bulk modulus of the annulus

Ma    = m.FV(value = 0.0)      # mass of fluid in the annulus
Md    = m.FV(value = 0.0)      # mass of fluid in the drill string
M     = m.FV(value = 4.3e+8)   # Total mass of fluid (Ma + Md)

rhoA  = m.FV(value = 1580.0)   # density of the fluid in the drill string
rhoD  = m.FV(value = 1580.0)   # density of the fluid in the annulus
rhoM  = m.FV(value = 1580.0)   # density of the mud used

Fa    = m.Const(value = 2.0e+9)   # friction coefficient of the annulus
Fd    = m.Const(value = 5.0e+9)   # friction coefficient of the drill string

Qres  = m.FV(value = 0.0)      # volume flow rate from the reservoir

# Define Variables
Qpump = m.FV(value = 2000.0*(1e-3/60))   # volume flow rate through the pump
Qback = m.FV(value = 800.0*(1e-3/60))    # volume flow rate through the back pressure pump
Vd    = m.FV(value = 17.0)     # volume of the drill string
Va    = m.FV(value = 48.0)     # volume of the annulus

Pp    = m.Var(value = 40.0e+5)  # pump pressure inside the drill string
Pc    = m.Var(value = 10.0e+5)  # choke pressure inside the annulus
Qbit  = m.Var(value = 2000.0*(1e-3/60))   # volume flow rate through the drill bit
 # volume flow rate through the choke (Qpump+Qback)
Qchoke = m.Intermediate( Qbit + Qback )
#dPa = m.Var(value = 0.0)
  

Pbit  = m.Var(value = 0.0)      # Pressure at the drill bit / BHP
#dPa   = m.Var(value = 0.0)      # differential pressure across the annulus (?)

# Model equations
# Equation 5.1
#m.Equation( Pp.dt() == (betaD/Vd) * (Qpump - Qbit - Vd.dt()) )
m.Equation( Vd*Pp.dt() == (betaD) * (Qpump - Qbit) )

# Equation 5.2
#m.Equation( Pc.dt() == (betaA/Va) * (Qres + Qbit + Qback - Qchoke - Va.dt()) )
m.Equation( Va*Pc.dt() == (betaA) * (Qres + Qbit + Qback - Qchoke) )

# Equation 5.3
#m.Equation( Qbit.dt() == (1/M) * (Pp - Pc - dPa + (rhoD-rhoA)*g*h) )
#Using the equation by replacing dPa by the pressure balance in drillpipe
m.Equation( Qbit.dt() == (1/M) * (Pp - Pbit - Fd*(Qbit**2) + rhoD*g*h ) )

# Equation 5.6
m.Equation( Pbit == Pc + rhoA*g*h + Fa*(Qbit**2) )

# Other equations for ... BALANCE?
#m.Equation( dPa == Pbit - Pc )

# Solver options
m.options.IMODE = 4  # Dynamic Simulation


# Solve the model
m.solve()

plt.figure(3)
plt.subplot(3,1,1)
plt.plot(m.time, Pbit.VALUE, 'r', label='BHP(Pa)')
plt.plot(m.time, Pp.VALUE, 'b', label='Pump Pressure(Pa)')
plt.legend(loc='best')
plt.subplot(3,1,2)
plt.plot(m.time, Pc.VALUE, 'g', label ='Choke Pressure(Pa)')
plt.legend(loc='best')
plt.subplot(3,1,3)
plt.plot(m.time, np.array(Qbit.VALUE)*60000, 'k', label ='Qbit (l/min)')
plt.plot(m.time, np.array(Qchoke.VALUE)*60000, 'r', label ='Qchoke (l/min)')
plt.legend(loc='best')
plt.show()
