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

# Set up the model time array (10 minutes, 1 second intervals)
loops = 60*10+1
m.time = np.linspace(0, loops-1, loops)

# Constants
t0 = 23 + 273    # Initial temperature
tAmb = 23 + 273  # Ambient temperature
AmCp = 6.0e-4    # A/mCp from notes
sig = 5.67e-8    # Stefsn-Boltzmann constant
eps = 0.9        # Emissivity 
alpha = 0.005    # Heater scaling factor / mCp

U = 8         # Initial estimate of HT coefficient

# Parameters:
q0 = np.zeros(loops)
q0[:60]     = 0.0
q0[60:240]  = 100.0
q0[240:300] = 0.0
q0[300:360] = 100.0
q0[360:]    = 0.0
Q = m.Param(value=q0)

# Variables:
T = m.Var(value=t0)

# Equation
m.Equation(T.dt()==U*AmCp*(tAmb-T) \
	              + eps*sig*AmCp*(tAmb**4-T**4) \
	              + alpha*Q)

# Simulation mode
m.options.IMODE = 4

# Solve simulation model
m.solve()

# Create plot
plt.figure(1)
plt.ion()
plt.show()

T1 = np.zeros(loops)
T1[0] = a.T1

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
                    
        # Read temperature 
        T1[i] = a.T1

        # Write output (0-100)
        a.Q1(Q[i])

        # Plot
        plt.clf()
        ax=plt.subplot(2,1,1)
        ax.grid()

        plt.plot(m.time/60.0,np.array(T.value)-273.15,'r-')
        plt.plot(tm[0:i]/60.0,T1[0:i],'g--')
        
        plt.ylabel('Temperature (degC)')
        plt.legend(['Simulation', 'Actual'])

        plt.subplot(2,1,2)
        plt.plot(m.time/60.0,np.array(Q.value),'b-')
        plt.ylabel('Heater (%)')
        plt.xlabel('Time (min)')

        plt.draw()
        plt.pause(0.05)

    # Turn off heaters
    a.Q1(0)

    # Turn LED on
    print('LED Off')
    a.LED(0)

    # Save figure
    plt.savefig('TCLab_eval.png')
        
# Allow user to end loop with Ctrl-C           
except KeyboardInterrupt:
    # Disconnect from Arduino
    a.Q1(0)
    print('Shutting down')
    
    # Turn LED off
    print('LED Off')
    a.LED(0)

    a.close()
    plt.savefig('TCLab_eval.png')
    
# Make sure serial connection still closes when there's an error
except:           
    # Disconnect from Arduino
    a.Q1(0)
    print('Error: Shutting down')

    # Turn LED off
    print('LED Off')
    a.LED(0)

    a.close()
    plt.savefig('TCLab_eval.png')
    raise
