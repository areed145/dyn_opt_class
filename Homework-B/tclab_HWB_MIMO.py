import tclab
import numpy as np
import time
import matplotlib.pyplot as plt
from gekko import GEKKO

# Connect to Arduino
a = tclab.TCLab()

# Get Version
print(a.version)

# Turn LED on
print('LED On')
a.LED(100)

# Set up the heater model in GEKKO
m = GEKKO(remote=False)
print(m.path)

# Set up the model time array (10 minutes, 1 second intervals)
loops = 60*10+1
m.time = np.linspace(0, loops-1, loops)

# Constants
t1_0 = a.T1 + 273  # Initial T1 temperature
t2_0 = a.T2 + 273  # Initial T2 temperature
tAmb = 23 + 273  # Ambient temperature
area = 1.0e-3      # Area not between heaters
shArea = 2.0e-4    # Shared area between heaters
mass = 4.0e-3      # Mass
Cp = 500.0         # Heat capacity
sig = 5.67e-8    # Stefan-Boltzmann constant
eps = 0.9        # Emissivity 
alpha = 0.01     # Heater scaling factor / mCp

U = 8         # Initial estimate of HT coefficient

# Parameters:
q1_0 = np.zeros(loops)
q1_0[:30]  = 0.0
q1_0[30:]  = 100.0
Q1 = m.Param(value=q1_0)

q2_0 = np.zeros(loops)
q2_0[:300]  = 0.0
q2_0[300:]  = 100.0
Q2 = m.Param(value=q2_0)


# Variables:
T1 = m.Var(value=t1_0)
T2 = m.Var(value=t2_0)

# Intermediates
#Qc11 = m.Intermediate(U*area*(tAmb-T1))
#Qc12 = m.Intermediate(U*shArea*(T2-T1))
#Qc22 = m.Intermediate(U*area*(tAmb-T2))
#Qr11 = m.Intermediate(eps*sig*area*(tAmb**4-T1**4))
#Qr12 = m.Intermediate(eps*sig*shArea*(T2**4-T1**4))
#Qr22 = m.Intermediate(eps*sig*area*(tAmb**4-T2**4))

Qc11 = m.Var(value=U*area*(tAmb-T1))
Qc12 = m.Var(value=U*shArea*(T2-T1))
Qc22 = m.Var(value=U*area*(tAmb-T2))
Qr11 = m.Var(value=eps*sig*area*(tAmb**4-T1**4))
Qr12 = m.Var(value=eps*sig*shArea*(T2**4-T1**4))
Qr22 = m.Var(value=eps*sig*area*(tAmb**4-T2**4))

# Equations
m.Equations([mass*Cp*T1.dt()==Qc11 + Qr11 + Qc12 + Qr12 + alpha*Q1, \
             mass*Cp*T2.dt()==Qc22 + Qr22 - Qc12 - Qr12 + alpha*Q2,
             Qc11==U*area*(tAmb-T1), \
             Qc12==U*shArea*(T2-T1), \
             Qc22==U*area*(tAmb-T2), \
             Qr11==eps*sig*area*(tAmb**4-T1**4), \
             Qr12==eps*sig*shArea*(T2**4-T1**4), \
             Qr22==eps*sig*area*(tAmb**4-T2**4)     ])

# Simulation mode
m.options.IMODE = 4

# Solve simulation model
m.solve()

# Create plot
plt.figure(1)
plt.ion()

#r = input('press a key to continue...')

aT1 = np.zeros(loops)
aT1[0] = a.T1
aT2 = np.zeros(loops)
aT2[0] = a.T2
aQ1 = np.zeros(loops)
aQ1[0] = Q1[0]
aQ2 = np.zeros(loops)
aQ2[0] = Q2[0]

tm = np.zeros(loops)

# Main Loop
start_time = time.time()
prev_time = start_time
try:
    for i in range(1,loops):
        # Sleep time
        sleep_max = 1.0
        sleep = sleep_max - (time.time() - prev_time)
        if sleep>=0.01:
            time.sleep(sleep)
        else:
            time.sleep(0.01)

        # Record time and change in time
        t = time.time()
        dt = t - prev_time
        prev_time = t
        tm[i] = t - start_time
                    
        # Read temperatures
        aT1[i] = a.T1
        aT2[i] = a.T2

        # Write outputs (0-100)
        a.Q1(Q1[i])
        a.Q2(Q2[i])

        aQ1[i] = Q1[i]
        aQ2[i] = Q2[i]
		
        # Plot
        plt.clf()
        plt.subplot(2,1,1)
        plt.plot(tm[0:i]/60.0,aT1[0:i],'r--',label='T1 Act')
        plt.plot(tm[0:i]/60.0,aT2[0:i],'b--',label='T2 Act')

        plt.subplot(2,1,2)
        plt.plot(tm[0:i]/60.0,aQ1[0:i],'r--',label='Q1 Act')
        plt.plot(tm[0:i]/60.0,aQ2[0:i],'b--',label='Q2 Act')
		
        plt.draw()
        plt.pause(0.05)

    # Turn off heaters
    a.Q1(0)
    a.Q2(0)

    # Turn LED on
    print('LED Off')
    a.LED(0)

    plt.subplot(2,1,1)
    #ax.grid()
    plt.plot(m.time/60.0,np.array(T1.value)-273.15,'r-',label='T1 Sim')
    plt.plot(m.time/60.0,np.array(T2.value)-273.15,'b-',label='T2 Sim')
    plt.ylabel('Temperature (degC)')
    plt.legend(loc='best')

    plt.subplot(2,1,2)
    plt.plot(m.time/60.0,np.array(Q1.value),'r-',label='Q1 Sim')
    plt.plot(m.time/60.0,np.array(Q2.value),'b-',label='Q2 Sim')
    plt.ylabel('Heater (%)')
    plt.xlabel('Time (min)')
    plt.legend(loc='best')

    plt.show()

    # Save figure
    plt.savefig('TCLab_MIMO.png')
        
# Allow user to end loop with Ctrl-C           
except KeyboardInterrupt:
    # Disconnect from Arduino
    a.Q1(0)
    a.Q2(0)
    print('Shutting down')
    
    # Turn LED off
    print('LED Off')
    a.LED(0)

    a.close()
    plt.savefig('TCLab_MIMO.png')
    
# Make sure serial connection still closes when there's an error
except:           
    # Disconnect from Arduino
    a.Q1(0)
    a.Q2(0)
    print('Error: Shutting down')

    # Turn LED off
    print('LED Off')
    a.LED(0)

    a.close()
    plt.savefig('TCLab_MIMO.png')
    raise
