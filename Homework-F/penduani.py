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
import matplotlib.animation as anim

#%% Non-model parameters
rmt = False  # Solve local or remote

tf = 10.0     # Final time
npts = 100    # time steps

#%% Specify model
m = GEKKO(remote=rmt)
m.time = np.linspace(0, tf, npts+1)  # Model Timeline

# Model Constants and Parameters
m1 = 10.0       # cart mass
m2 = 50.0       # pendulum mass
tmax = 8.2      # terminal time

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
u.DCOST  = 0.01     # move penalty

# CV tuning parameters

# Solver options
m.options.IMODE = 6     # Dynamic Optimization (Control)
m.options.CV_TYPE = 2   # L1 or L2 Norm

# Solve the model
m.solve()

#%% Display the results
plt.figure()
#plt.figure(figsize=(12,10))

plt.subplot(2,2,1)
plt.plot(m.time,u.value,'m',lw=2)
plt.legend([r'$u$'],loc=1)
plt.ylabel('Force')
plt.xlabel('Time')
plt.xlim(m.time[0],m.time[-1])

plt.subplot(2,2,2)
plt.plot(m.time,v.value,'g',lw=2)
plt.legend([r'$v$'],loc=1)
plt.ylabel('Velocity')
plt.xlabel('Time')
plt.xlim(m.time[0],m.time[-1])

plt.subplot(2,2,3)
plt.plot(m.time,y.value,'r',lw=2)
plt.legend([r'$y$'],loc=1)
plt.ylabel('Position')
plt.xlabel('Time')
plt.xlim(m.time[0],m.time[-1])

plt.subplot(2,2,4)
plt.plot(m.time,theta.value,'y',lw=2)
plt.plot(m.time,q.value,'c',lw=2)
plt.legend([r'$\theta$',r'$q$'],loc=1)
plt.ylabel('Angle')
plt.xlabel('Time')
plt.xlim(m.time[0],m.time[-1])

plt.rcParams['animation.html'] = 'html5'

x1 = y.value
y1 = np.zeros(len(m.time))

#suppose that l = 1
x2 = 1*np.sin(theta.value)+x1
x2b = 1.05*np.sin(theta.value)+x1
y2 = -1*np.cos(theta.value)+y1
y2b = -1.05*np.cos(theta.value)+y1

fig = plt.figure(figsize=(8,6.4))
ax = fig.add_subplot(111,autoscale_on=False,\
                     xlim=(-1.5,0.5),ylim=(-1.2,0.4))
ax.set_xlabel('position')
ax.get_yaxis().set_visible(False)

crane_rail, = ax.plot([-1.5,0.5],[0.2,0.2],'k-',lw=4)
start, = ax.plot([-1,-1],[-1.5,1],'k:',lw=2)
objective, = ax.plot([0,0],[-1.5,1],'k:',lw=2)
mass1, = ax.plot([],[],linestyle='None',marker='s',\
                 markersize=40,markeredgecolor='k',\
                 color='orange',markeredgewidth=2)
mass2, = ax.plot([],[],linestyle='None',marker='o',\
                 markersize=20,markeredgecolor='k',\
                 color='orange',markeredgewidth=2)
line, = ax.plot([],[],'o-',color='orange',lw=4,\
                markersize=6,markeredgecolor='k',\
                markerfacecolor='k')
time_template = 'time = %.1fs'
time_text = ax.text(0.05,0.9,'',transform=ax.transAxes)
start_text = ax.text(-1.06,-1.1,'start',ha='right')
end_text = ax.text(0.06,-1.1,'objective',ha='left')

def init():
    mass1.set_data([],[])
    mass2.set_data([],[])
    line.set_data([],[])
    time_text.set_text('')
    return line, mass1, mass2, time_text

def animate(i):
    mass1.set_data([x1[i]],[y1[i]+0.1])
    mass2.set_data([x2b[i]],[y2b[i]])
    line.set_data([x1[i],x2[i]],[y1[i],y2[i]])
    time_text.set_text(time_template % m.time[i])
    return line, mass1, mass2, time_text

ani_a = anim.FuncAnimation(fig, animate, \
         np.arange(1,len(m.time)), \
         interval=40,blit=False,init_func=init)

#ani_a.save('Pendulum_Control.mp4',fps=30)
ani_a.save('Pendulum_Control.htm',fps=30)

plt.show()