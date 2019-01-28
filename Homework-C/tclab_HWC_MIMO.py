import tclab
import numpy as np
import time
import matplotlib.pyplot as plt
#from gekko import GEKKO

# Constants
tAmb = 24.5 + 273.15  # Ambient temperature
area = 1.0e-3         # Area not between heaters
shArea = 4.0e-4       # Shared area between heaters
mass = 4.0e-3         # Mass
Cp = 500.0            # Heat capacity
sig = 5.67e-8         # Stefan-Boltzmann constant
eps = 0.9             # Emissivity 
alpha = 0.01          # Heater scaling factor / mCp
q2factor = 0.75       # Heater 2 factor

U = 12                # Initial estimate of HT coefficient
U12 = 30              # Initial estimate of T1-T2 HT coefficient

loops = 60*11+1

# Create plot
plt.figure(2)
plt.subplot(2,1,1)
plt.ylabel('Temperature (degC)')
plt.subplot(2,1,2)
plt.ylabel('Heater (%)')
plt.xlabel('Time (min)')
plt.ion()

# Connect to Arduino
a = tclab.TCLab()

# Turn LED on
print('LED On')
a.LED(100)

# Let the TIs cool off
start_time = time.time()
a.Q1(0.0)
a.Q2(0.0)
while a.T1>25 or a.T2>25:
    print('\r', int(time.time()-start_time), '-- T1 = ', a.T1, 'T2 =', a.T2, \
          end='', flush=True)
    time.sleep(5)
    a.LED(0)
    time.sleep(5)
    a.LED(100)

aT1 = np.zeros(loops)
aT1[0] = a.T1

aT2 = np.zeros(loops)
aT2[0] = a.T2

# test steps
aQ1 = np.zeros(loops)
aQ1[10:] = 100.0
aQ1[200:] = 5.0
aQ1[400:] = 70.0

aQ2 = np.zeros(loops)
aQ2[100:] = 50.0
aQ2[300:] = 100.0
aQ2[500:] = 10.0

tm = np.zeros(loops)

# Main Loop
t = time.time()
start_time = t
end_time = t + loops
prev_time = start_time
i = 1
try:
    #for i in range(1,loops):
    while t < end_time and i < loops:
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
        a.Q1(aQ1[i])
        a.Q2(aQ2[i])
        
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

        i += 1

    # Turn off heaters
    a.Q1(0)
    a.Q2(0)

    # Turn LED on
    print('LED Off')
    a.LED(0)

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
    plt.savefig('TCLab_MIMO_fail.png')
    
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
    plt.savefig('TCLab_MIMO_fail.png')
    raise

# save training data to file
data = np.vstack((tm,aT1,aT2,aQ1,aQ2)).T
np.savetxt('mimo_test_data.csv',data,header='tm,T1,T2,Q1,Q2', \
           comments='',delimiter=',')
    
input('Press any key to exit.')
