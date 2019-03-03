"""
MPC - Pendulum Example
http://apmonitor.com/do/index.php/Main/ControlTypes
Solve with GEKKO

Problem Definition:

A pendulum is described by the following dynamic equations:

   | distdot  |   | 0 1  0  0 | | dist  |   |  0 |
   | velDot   | = | 0 0 eps 0 |*| vel   | + |  1 | * U
   | thetaDot |   | 0 0  0  1 | | theta |   |  0 |
   | rocDot   |   | 0 0 -1  0 | | ROC   |   | -1 |

 Parameters
   m1 = 10
   m2 = 1
   eps = m2/(m1+m2)

 Variables
   dist  = cart position
   vel   = vart velocity
   theta = pendulum angle
   ROC   = rate of change of pendulum angle

"""

#%%Import packages
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

#%% Non-model parameters
rmt = False  # Solve local or remote

tf = 10.0     # Final time
npts = 100    # time steps

#%% Specify model
m = GEKKO(remote=rmt)
m.time = np.linspace(0, tf, npts+1)  # Model Timeline

# Model Constants and Parameters
m1 = 10.0       # cart mass
m2 = 5.0        # pendulum mass
tmax = 6.2      # terminal time

eps = m.Param(value=m2/(m1+m2))

# Define Variables
y = m.Var(value=-1.0)       # cart position
v = m.Var(value=0.0)        # cart velocity
theta = m.Var(value=0.0)    # pendulum angle
q = m.Var(value=0.0)        # pendulum angle rate of change

u = m.MV(value=0.0)         # force on cart

# Objective
term = m.Param(value=np.array([int(t>=tmax) for t in m.time]))
m.Obj(term*(y*y + v*v + theta*theta + q*q))

# Model Equations
m.Equation( y.dt() == v )
m.Equation( v.dt() == eps * theta + u )
m.Equation( theta.dt() == q )
m.Equation( q.dt() == -theta - u )

# Tuning

# MV tuning parameters
u.STATUS = 1        # turn MV ON
#u.DCOST  = 0.001    # move penalty

# CV tuning parameters

# Solver options
m.options.IMODE = 6     # Dynamic Optimization (Control)
m.options.CV_TYPE = 2   # L1 or L2 Norm

# Solve the model
m.solve()

#%% Display the results
plt.figure(2)
plt.subplot(4,1,1)
plt.plot(m.time,u.value,'r-',label=r'$u$')
plt.ylabel('Force')
plt.legend(loc='best')
plt.subplot(4,1,2)
plt.plot(m.time,v.value,'b--',label=r'$v$')
plt.legend(loc='best')
plt.ylabel('Velocity')
plt.subplot(4,1,3)
plt.plot(m.time,y.value,'g:',label=r'$y$')
plt.legend(loc='best')
plt.ylabel('Position')
plt.subplot(4,1,4)
plt.plot(m.time,theta.value,'m-',label=r'$\theta$')
plt.plot(m.time,q.value,'k.-',label=r'$q$')
plt.legend(loc='best')
plt.ylabel('Angle')
plt.xlabel('Time')
plt.show()