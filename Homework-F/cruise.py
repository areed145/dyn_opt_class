"""
MPC - Cruise Control Example
http://apmonitor.com/do/index.php/Main/DynamicControl
Solve with GEKKO

Problem Definition:
 Constants
   m = 500 ! Mass (kg)
 Parameters
   b = 50  ! Resistive coefficient (N-s/m)  
   K = 0.8 ! Gain (m/s-%pedal)
   p = 0 >= 0 <= 100  ! Gas pedal position (%)
 Variables
   v = 0 ! initial condition
 Equations
   m * $v = -v * b + K * b * p

"""

#%%Import packages
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

#%% Non-model parameters
rmt = False  # Solve local or remote

tf = 25.0     # Final time
npts = 50     # time steps

#%% Specify model
m = GEKKO(remote=rmt)
m.time = np.linspace(0, tf, npts+1)  # Model Timeline

# Model Constants and Parameters
mass = 500      # vehicle mass

b = m.Param(value=50.0)     # resistive coefficient
k = m.Param(value=0.8)      # gain

# Define Variables
p = m.MV(value=0.0, lb=0.0, ub=100.0)   # gas pedal
v = m.CV(value=0.0)                     # vehicle velocity

# Model Equations
m.Equation( mass*v.dt() == -v*b + k*b*p )


# Tuning

# MV tuning parameters
p.STATUS = 1        # turn MV ON
p.DCOST  = 0.5      # move penalty
p.DMAX   = 100.0    # maximum move

# CV tuning parameters
v.STATUS = 1        # turn CV ON
v.SP   = 40.0       # setpoint for L2 norm
#v.SPLO = 40.0       # low setpoint for L1 norm
#v.SPHI = 45.0       # high setpoint for L1 norm
v.TR_INIT = 1        # initial equal to the current value on coldstart
v.TAU     = 4.0     # speed of SP response


# Solver options
m.options.IMODE = 6     # Dynamic Optimization (Control)
m.options.CV_TYPE = 2   # L1 or L2 Norm

# Solve the model
m.solve()

#%% Display the results
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(m.time,p.value,'k-',label=r'$p$')
plt.ylabel('gas')
plt.subplot(2,1,2)
plt.plot(m.time,v.value,'r--',label=r'$v$')
plt.ylabel('velocity')
plt.xlabel('time')
plt.show()