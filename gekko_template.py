"""
Describe this program...

Problem Definition...

"""

#%%Import packages
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

#%% Non-model parameters
rmt = False  # Solve local or remote
tf = 10.0     # Final time
npts = 100    # time steps
tmax = 8.0    # end point for "terminal" kinds of limits 

#%% Specify model
m = GEKKO(remote=rmt)
m.time = np.linspace(0, tf, npts+1)  # Model Timeline

# Model Constants and Parameters
c1 = 1.0       # constant #1

# Model Parameters
p1 = m.Param(value=0.0)     # parameter #1

# Model Variables
y = m.Var(value=-1.0)       # general variable

u = m.MV(value=1.0)         # MV

x = m.CV(value=1.0)         # CV

# Objective
term = m.Param(value=np.array([int(t>=tmax) for t in m.time]))
m.Obj(term*y*y)
m.Obj(term*x*x)
m.Obj(term*u*u)

# Model Equations
m.Equation( y.dt() == -y + u )
m.Equation( 5.0*x.dt() == -x + u )

# Tuning

# MV tuning parameters
u.STATUS = 1        # turn MV ON
u.DCOST  = 0.01     # move penalty
u.DMAX   = 100.0    # maximum move

# CV tuning parameters
x.STATUS = 1        # turn CV ON
x.SP   = 0.0        # setpoint for L2 norm
x.SPLO = -1.0       # low setpoint for L1 norm
x.SPHI = 1.0       # high setpoint for L1 norm
x.TR_INIT = 1       # initial equal to the current value on coldstart
x.TAU     = 2.0     # speed of SP response

# Solver options
m.options.IMODE = 6     # Dynamic Optimization (Control)
m.options.CV_TYPE = 2   # L1 or L2 Norm

# Solve the model
m.solve()

#%% Display the results
plt.figure()

plt.subplot(2,1,1)
plt.plot(m.time,u.value,'k-',label=r'$u$')
plt.legend(loc='best')
plt.ylabel('MV')
plt.subplot(2,1,2)
plt.plot(m.time,y.value,'r--',label=r'$y$')
plt.plot(m.time,x.value,'g--',label=r'$x$')
plt.legend(loc='best')
plt.ylabel('CV')
plt.xlabel('time')
plt.show()
