import numpy as np
import time
import matplotlib.pyplot as plt
import random
# get gekko package with:
#   pip install gekko
from gekko import GEKKO
# get tclab package with:
#   pip install tclab
from tclab import TCLab

# Connect to Arduino
a = TCLab()

# Turn LED on
print('LED On')
a.LED(100)

# Final time
tf = 10 # min

# number of data points (1 pt every 3 seconds)
n = tf * 20 + 1

# Configure heater levels
# Percent Heater (0-100%)
Q1s = np.zeros(n)
Q2s = np.zeros(n)

# Heater random steps every 30 steps
# Alternate steps by Q1 and Q2
for i in range(0,130):
    if i%30==0:
        Q1s[i:i+30] = random.random() * 100.0
    if (i+5)%30==0:
        Q2s[i:i+30] = random.random() * 100.0

# Heater random steps every 10 steps
for i in range(130,200):
    if i%10==0:
        Q1s[i:i+10] = random.random() * 100.0
    if (i+5)%10==0:
        Q2s[i:i+10] = random.random() * 100.0
        
# Record initial temperatures (degC)
T1m = a.T1 * np.ones(n)
T2m = a.T2 * np.ones(n)

# Store MHE values for plots
Tmhe1 = T1m[0] * np.ones(n)
Tmhe2 = T2m[0] * np.ones(n)
Umhe = 10.0 * np.ones(n)
amhe1 = 0.01 * np.ones(n)
amhe2 = 0.0075 * np.ones(n)
solve = 0 * np.ones(n)
err1 = 0 * np.ones(n)
err2 = 0 * np.ones(n)
e1 = 0
e2 = 0

#########################################################
# Initialize Model as Estimator
#########################################################

m = GEKKO(name='tclab-mhe',remote=True)

# 60 second time horizon, 20 steps
m.time = np.linspace(0,60,21)

# Parameters to Estimate
U = m.FV(value=10,name='u')
U.STATUS = 1  # don't estimate initially
U.LOWER = 5
U.UPPER = 15

alpha1 = m.FV(value=0.01,name='a1')   # W / % heater
alpha1.STATUS = 1  # don't estimate initially
alpha1.LOWER = 0.003
alpha1.UPPER = 0.03

alpha2 = m.FV(value=0.0075,name='a2') # W / % heater
alpha2.STATUS = 1  # don't estimate initially
alpha2.LOWER = 0.002
alpha2.UPPER = 0.02

# Measured inputs
Q1 = m.MV(value=0,name='q1')
Q1.STATUS = 0  # don't estimate
Q1.FSTATUS = 1 # receive measurement

Q2 = m.MV(value=0,name='q2')
Q2.STATUS = 0  # don't estimate
Q2.FSTATUS = 1 # receive measurement

# Measurements for model alignment
TC1 = m.CV(value=T1m[0],name='tc1')
TC1.STATUS = 1     # minimize error between simulation and measurement
TC1.FSTATUS = 1    # receive measurement
TC1.MEAS_GAP = 0.1 # measurement deadband gap
TC1.LOWER = 0
TC1.UPPER = 200

TC2 = m.CV(value=T2m[0],name='tc2')
TC2.STATUS = 1     # minimize error between simulation and measurement
TC2.FSTATUS = 1    # receive measurement
TC2.MEAS_GAP = 0.1 # measurement deadband gap
TC2.LOWER = 0
TC2.UPPER = 200

Ta = m.Param(value=23.0+273.15)     # K
mass = m.Param(value=4.0/1000.0)    # kg
Cp = m.Param(value=0.5*1000.0)      # J/kg-K    
A = m.Param(value=10.0/100.0**2)    # Area not between heaters in m^2
As = m.Param(value=2.0/100.0**2)    # Area between heaters in m^2
eps = m.Param(value=0.9)            # Emissivity
sigma = m.Const(5.67e-8)            # Stefan-Boltzmann

# Heater temperatures
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

# Global Options
m.options.IMODE   = 5 # MHE
m.options.EV_TYPE = 2 # Objective type
m.options.NODES   = 2 # Collocation nodes
m.options.SOLVER  = 3 # IPOPT

##################################################################

# Create plot
plt.figure(figsize=(10,7))
plt.ion()
plt.show()

# Main Loop
start_time = time.time()
prev_time = start_time
tm = np.zeros(n)

try:
    for i in range(1,n):
        # Sleep time
        sleep_max = 3.0
        sleep = sleep_max - (time.time() - prev_time)
        if sleep>=0.01:
            time.sleep(sleep-0.01)
        else:
            time.sleep(0.01)

        # Record time and change in time
        t = time.time()
        dt = t - prev_time
        prev_time = t
        tm[i] = t - start_time

        # Read temperatures in Celsius 
        T1m[i] = a.T1
        T2m[i] = a.T2

        # Insert measurements
        TC1.MEAS = T1m[i]
        TC2.MEAS = T2m[i]
        Q1.MEAS = Q1s[i-1]
        Q2.MEAS = Q2s[i-1]

        # Start estimating U after 10 cycles (20 sec)
        if i==10:
            U.STATUS = 1
            alpha1.STATUS = 1
            alpha2.STATUS = 1
        
        # Predict Parameters and Temperatures with MHE
        m.solve() 

        if m.options.APPSTATUS == 1:
            # Retrieve new values
            Tmhe1[i] = TC1.MODEL
            Tmhe2[i] = TC2.MODEL
            Umhe[i]  = U.NEWVAL
            amhe1[i] = alpha1.NEWVAL
            amhe2[i] = alpha2.NEWVAL
            solve[i] = 1
        else:
            # Solution failed, copy prior solution
            Tmhe1[i] = Tmhe1[i-1]
            Tmhe2[i] = Tmhe1[i-1]
            Umhe[i]  = Umhe[i-1]
            amhe1[i] = amhe1[i-1]
            amhe2[i] = amhe2[i-1]
            solve[i] = 0
        
        # Write new heater values (0-100)
        a.Q1(Q1s[i])
        a.Q2(Q2s[i])

        # determine cumulative error
        e1 += abs(T1m[i] - Tmhe1[i])
        e2 += abs(T2m[i] - Tmhe2[i])
        err1[i] = e1
        err2[i] = e2

        # Plot
        plt.clf()
        ax=plt.subplot(4,1,1)
        ax.grid()
        plt.plot(tm[0:i],T1m[0:i],'rx',label=r'$T_1$ measured')
        plt.plot(tm[0:i],Tmhe1[0:i],'r--',label=r'$T_1$ MHE')
        plt.plot(tm[0:i],T2m[0:i],'bx',label=r'$T_2$ measured')
        plt.plot(tm[0:i],Tmhe2[0:i],'b--',label=r'$T_2$ MHE')
        plt.ylabel('Temperature (degC)')
        plt.legend(loc=2)
        ax=plt.subplot(4,1,2)
        ax.grid()
        plt.plot(tm[0:i],Umhe[0:i],'k:',label='Heat Transfer Coeff')       
        plt.plot(tm[0:i],amhe1[0:i]*1000,'r--',label=r'$\alpha_1$x1000')
        plt.plot(tm[0:i],amhe2[0:i]*1000,'b--',label=r'$\alpha_2$x1000')
        plt.ylabel('Parameters')
        plt.legend(loc='best')
        ax=plt.subplot(4,1,3)
        ax.grid()
        plt.plot(tm[0:i],Q1s[0:i],'r-',label=r'$Q_1$')
        plt.plot(tm[0:i],Q2s[0:i],'b-',label=r'$Q_2$')
        plt.ylabel('Heaters')
        plt.xlabel('Time (sec)')
        plt.legend(loc='best')
        ax=plt.subplot(4,1,4)
        ax.grid()
        plt.plot(tm[0:i],err1[0:i],'r-',label='')
        plt.plot(tm[0:i],err2[0:i],'b-',label='')
        plt.plot(tm[0:i],solve[0:i],'k:',label='')
        plt.ylabel('Error / Status')
        plt.xlabel('Time (sec)')
        plt.legend(loc='best')
        plt.draw()
        plt.pause(0.05)

    # Turn off heaters
    a.Q1(0)
    a.Q2(0)
    # Save figure
    plt.savefig('tclab_mhe_C.png')

    # Turn LED off
    print('LED Off')
    a.LED(0)
    
# Allow user to end loop with Ctrl-C           
except KeyboardInterrupt:
    # Disconnect from Arduino
    a.Q1(0)
    a.Q2(0)
    print('Shutting down')
    a.close()
    plt.savefig('tclab_mhe_C.png')

    # Turn LED off
    print('LED Off')
    a.LED(0)
    
# Make sure serial connection still closes when there's an error
except:           
    # Disconnect from Arduino
    a.Q1(0)
    a.Q2(0)
    print('Error: Shutting down')
    a.close()
    plt.savefig('tclab_mhe_C.png')

    # Turn LED off
    print('LED Off')
    a.LED(0)

    raise