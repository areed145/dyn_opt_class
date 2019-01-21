# -*- coding: utf-8 -*-
import tclab  # pip install tclab
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from gekko import GEKKO

#initialize GEKKO model
m = GEKKO(remote=False)

# Run time in minutes
run_time = 10.0

#model discretized time
t_max = 60*run_time+1  # Number of second time points (10min)
m.time = np.linspace(0,t_max-1,t_max) # Time vector

# Parameters
Tamb = 23               # K
U = 5.0                 # W/m^2-K
mass = 4.0/1000.0          # kg
Cp = 0.5 * 1000.0       # J/kg-K    
A = 10.0 / 100.0**2     # Area in m^2
As = 4.0 / 100.0**2     # Area in m^2
alpha1 = 0.0100         # W / % heater 1
alpha2 = 0.0075         # W / % heater 2
eps = 0.9               # Emissivity
sigma = 5.67e-8         # Stefan-Boltzman

# Percent Heater (0-100%)
Q1to = 10               # at 10 sec
Q2to = 300              # at 5.0 min (300 sec)

# Heaters as time-varying inputs
Q1v = np.zeros(int(t_max))
Q1v[Q1to:] = 100.0 
Q2v = np.zeros(int(t_max))
Q2v[Q2to:] = 100.0
Q1_ = m.Param(value=Q1v) # Percent Heater (0-100%)
Q2_ = m.Param(value=Q2v) # Percent Heater (0-100%)

Tamb_ = m.Param(value=Tamb)
U_ =  m.Param(value=U)
mass_ = m.Param(value=mass)
Cp_ = m.Param(value=Cp)
A_ = m.Param(value=A)
As_ = m.Param(value=As)
alpha1_ = m.Param(value=alpha1)
alpha2_ = m.Param(value=alpha2)
eps_ = m.Param(value=eps)
sigma_ = m.Const(sigma)

# Temperature states as GEKKO variables
T1_ = m.Var(value=Tamb_+273.15)
T2_ = m.Var(value=Tamb_+273.15)

# Between two heaters
Q_C12_ = m.Intermediate(U_*As_*(T2_-T1_)) # Convective
Q_R12_ = m.Intermediate(eps_*sigma_*As_*(T2_**4-T1_**4)) # Radiative

# Nonlinear Energy Balances
m.Equation(T1_.dt() == (1.0/(mass_*Cp_))*(U_*A_*(Tamb_-T1_) \
                 + eps_*sigma_*A_*(Tamb_**4-T1_**4) \
                 + Q_C12_+Q_R12_+alpha1_*Q1_))

m.Equation(T2_.dt() == (1.0/(mass_*Cp_))*(U_*A_*(Tamb_-T2_) \
                 + eps_*sigma_*A_*(Tamb_**4-T2_**4) \
                 - Q_C12_-Q_R12_+alpha2_*Q2_))

#simulation mode
m.options.IMODE = 4

#simulation model
m.solve()

# define energy balance model
def heat(x,t,Q1,Q2):

    # Temperature States 
    T1 = x[0]
    T2 = x[1]

    # Heat Transfer Exchange Between 1 and 2
    conv12 = U*As*(T2-T1)
    rad12  = eps*sigma*As*(T2**4-T1**4)

    # Nonlinear Energy Balances
    dT1dt = (1.0/(mass*Cp))*(U*A*(Tamb-T1) \
            + eps*sigma*A*(Tamb**4-T1**4) \
            + conv12+rad12+alpha1*Q1)
    dT2dt = (1.0/(mass*Cp))*(U*A*(Tamb-T2) \
            + eps*sigma*A*(Tamb**4-T2**4) \
            - conv12-rad12+alpha2*Q2)
    return [dT1dt,dT2dt]

# save txt file
def save_txt(t,u1,u2,y1,y2):
    data = np.vstack((t,u1,u2,y1,y2))  # vertical stack
    data = data.T                 # transpose data
    top = 'Time (sec), Heater 1 (%), Heater 2 (%), ' \
        + 'Temperature 1 (degC), Temperature 2 (degC)'
    np.savetxt('data.txt',data,delimiter=',',header=top,comments='')
    
# Connect to Arduino
a = tclab.TCLab()

# Turn LED on
print('LED On')
a.LED(100)

# initialize lists
ltm = [] # time
lQ1 = [] # set point 1
lQ2 = [] # set point 2
lT1p = [] # pred T1
lT2p = [] # pred T2
lT1m = [] # measured T1
lT2m = [] # measured T2
lerr = [] # error

print('Running Main Loop. Ctrl-C to end.')
print('  Time   Q1     Q2    T1     T2')

# Create plot
plt.figure(figsize=(10,7))
plt.ion()
plt.show()

# Main Loop
start_time = time.time()
prev_time = start_time

t_elapsed = 0
Q1 = 0
Q2 = 0
T1pp = a.T1
T2pp = a.T2
T1mp = a.T1
T2mp = a.T2
errp = 0

try:
    while t_elapsed <= t_max:
        # Sleep time
        sleep_max = 1.0
        sleep = sleep_max - (time.time() - prev_time)
        if sleep >= 0.01:
            time.sleep(sleep-0.01)
        else:
            time.sleep(0.01)
            
        # Record time and change in time
        t = time.time()
        dt = t - prev_time
        prev_time = t
        t_elapsed = t - start_time
        
        # Heater input
        if t_elapsed >= Q1to:
            Q1 = 100.0
        if t_elapsed >= Q2to:
            Q2 = 100.0        

        # Read temperatures in Kelvin 
        T1m = a.T1
        T2m = a.T2

        # Simulate one time step with Energy Balance
        Tnext = odeint(heat,[T1pp+273.15,T2pp+273.15],[0,dt],args=(Q1,Q2))
        T1p = Tnext[1,0]-273.15
        T2p = Tnext[1,1]-273.15
        err = errp+(abs(T1p-T1m)+abs(T2p-T2m))*dt
                         
        ltm.append(t_elapsed)
        lQ1.append(Q1)
        lQ2.append(Q2)
        lT1p.append(T1p)
        lT2p.append(T2p)
        lT1m.append(T1m)
        lT2m.append(T2m)
        lerr.append(err)
                         
        T1pp = T1p
        T2pp = T2p
        T1mp = T1m
        T2mp = T2m
        errp = err

        # Write output (0-100)
        a.Q1(Q1)
        a.Q2(Q2)

        # Print line of data
        print('{:6.1f} {:6.2f} {:6.2f} {:6.2f} {:6.2f}'.format(t_elapsed, Q1, Q2, T1m, T2m))
        
        # Plot
        plt.clf()
        ax=plt.subplot(3,1,1)
        ax.grid()
        plt.plot(ltm,lT1m,'rx',label=r'$T_1$ measured')
        plt.plot(ltm,lT1p,'r-',label=r'$T_1$ energy balance')
        plt.plot(m.time,np.array(T1_.value)-273.15,'r:',label=r'$T_1$ GEKKO')
        plt.plot(ltm,lT2m,'bx',label=r'$T_2$ measured')
        plt.plot(ltm,lT2p,'b-',label=r'$T_2$ energy balance')
        plt.plot(m.time,np.array(T2_.value)-273.15,'b:',label=r'$T_2$ GEKKO')
        plt.ylabel('Temperature (degC)')
        plt.legend(loc=2)
        ax=plt.subplot(3,1,2)
        ax.grid()
        plt.plot(ltm,lerr,'k-',label='Energy Balance Error')
        plt.ylabel('Cumulative Error')
        plt.legend(loc='best')
        ax=plt.subplot(3,1,3)
        ax.grid()
        plt.plot(ltm,lQ1,'r-',label=r'$Q_1$')
        plt.plot(ltm,lQ2,'b-',label=r'$Q_2$')
        plt.plot(m.time,np.array(Q1_.value),'r:',label=r'$Q_1$ GEKKO')
        plt.plot(m.time,np.array(Q2_.value),'b:',label=r'$Q_2$ GEKKO')
        plt.ylabel('Heaters')
        plt.xlabel('Time (sec)')
        plt.legend(loc='best')
        plt.draw()
        plt.pause(0.05)

    # Turn off heaters
    a.Q1(0)
    a.Q2(0)
    # Save text file and plot at end
    save_txt(t_elapsed, Q1, Q2, T1m, T2m)
    # Save figure
    plt.savefig('test_Models.png')

# Allow user to end loop with Ctrl-C           
except KeyboardInterrupt:
    # Disconnect from Arduino
    a.Q1(0)
    a.Q2(0)
    print('Shutting down')
    a.close()
    save_txt(t_elapsed, Q1, Q2, T1m, T2m)
    plt.savefig('test_Models.png')

# Make sure serial connection still closes when there's an error
except:           
    # Disconnect from Arduino
    a.Q1(0)
    a.Q2(0)
    print('Error: Shutting down')
    a.close()
    save_txt(t_elapsed, Q1, Q2, T1m, T2m)
    plt.savefig('test_Models.png')
    raise
