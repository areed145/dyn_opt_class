
import tclab
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time

# Connect to Arduino
a = tclab.TCLab()

zeroK = 273.1
ambTemp = 23.0

# Defining energy balance model
def heat(x,t,Q1,Q2):
    #Parameters
    Ta = ambTemp + zeroK    #Ambient Temperature, K
    #U = 10.0                #Overall heat transfer coefficient, W/m^2-K
    U = 10.0                 #Overall heat transfer coefficient, W/m^2-K
    m = 4.0/1000            #Mass, kg
    Cp = 500.0              #Heat Capacity, J/kg-K
    #Cp = 550.0              #Heat Capacity, J/kg-K
    A = 10.0/100.0**2            #Surface area not between heat sinks,m^2
    alpha = 0.01            #Heater output coefficient, W/% heater
    eps = 0.9               #Emissivity
    sigma = 5.67e-8         #Stefan Boltzmann constant,W/m^2-K^4
    As = 2.0/100.0**2         #Surface area between heat sinks, m^2
    
    #Temperature state
    T1 = x[0]
    T2 = x[1]
    
    #Non-linear Energy Balance
    dT1dt = (1.0/(m*Cp))*(U*A*(Ta-T1) \
            + eps*sigma*A*(Ta**4-T1**4) \
            + U*As*(T2-T1) \
            + eps*sigma*As*(T2**4-T1**4) \
            + alpha*Q1)
    dT2dt = (1.0/(m*Cp))*(U*A*(Ta-T2) \
            + eps*sigma*A*(Ta**4-T2**4) \
            + U*As*(T1-T2) \
            + eps*sigma*As*(T1**4-T2**4) \
            + alpha*Q2)
    
    return [dT1dt,dT2dt]

# Run time in minutes
run_time = 10.0

# Number of cycles
loops = int(60.0*run_time)
tm = np.zeros(loops)

# Temperature (K)
Tsp1 = np.ones(loops) * ambTemp     # set point (degC)
T1 = np.ones(loops) * a.T1          # measured T (degC)
Tsp2 = np.ones(loops) * ambTemp     # set point (degC)
T2 = np.ones(loops) * a.T2          # measured T (degC)

# Predictions
Tp1 = np.ones(loops) * a.T1
error_eb1 = np.zeros(loops)
Tp2 = np.ones(loops) * a.T2
error_eb2 = np.zeros(loops)

# impulse tests (0 - 100%)
Q1 = np.ones(loops) * 0.0
Q1[6:] = 100.0                      # step up Heater 1 @ .1 min
Q2 = np.ones(loops) * 0.0
Q2[300:] = 100.0                    # step up Heater 2 @ 5 min

# Create plot
plt.figure(figsize=(10,7))
plt.ion()
plt.show()

#Main Loop
startTime = time.time()
previousTime = startTime

#Turn LED ON
a.LED(50)

try:
    for i in range(1,loops):
        #Sleep time
        sleepMax = 1.0
        sleep = sleepMax - (time.time() - previousTime)
        if (sleep > 0.01):
            time.sleep(sleep - 0.01)
        else:
            time.sleep(0.01)
        
        # Record the time and the change in time, dt
        timeNow = time.time()
        dt = timeNow - previousTime
        previousTime = timeNow
        
        tm[i] = timeNow - startTime
        
        #Read the temperatures in Celcius
        T1[i] = a.T1
        T2[i] = a.T2

        #Simulate one time step with energy balance
        Tnext = odeint(heat, [Tp1[i-1]+zeroK,Tp2[i-1]+zeroK]\
                       , [0,dt], args=(Q1[i-1],Q2[i-1]))
        Tp1[i] = Tnext[1, 0] - zeroK
        Tp2[i] = Tnext[1, 1] - zeroK
        error_eb1[i] = error_eb1[i-1] + abs(Tp1[i] - T1[i])*dt
        error_eb2[i] = error_eb2[i-1] + abs(Tp2[i] - T2[i])*dt
        
        # Write heater output
        a.Q1(Q1[i])
        a.Q2(Q2[i])

        #Plot
        plt.clf()
        ax = plt.subplot(3,1,1)
        ax.grid()
        plt.plot(tm[0:i], Tp1[0:i], 'b-', label=r'$T_1$ Energy balance')
        plt.plot(tm[0:i], T1[0:i], 'r+', label=r'$T_1$ Measured')
        plt.plot(tm[0:i], Tp2[0:i], 'g-', label=r'$T_2$ Energy balance')
        plt.plot(tm[0:i], T2[0:i], 'p', label=r'$T_2$ Measured')
        plt.ylabel('Temperature (degC)')
        plt.legend(loc=2)
        
        ax = plt.subplot(3,1,2)
        ax.grid()
        plt.plot(tm[0:i],error_eb1[0:i],'k-',label='Energy Balance 1')
        plt.plot(tm[0:i],error_eb2[0:i],'m-',label='Energy Balance 2')
        plt.ylabel('Cumulative Error')
        plt.legend(loc='best')
        
        ax = plt.subplot(3,1,3)
        ax.grid()
        plt.plot(tm[0:i],Q1[0:i],'r-',label=r'$Q_1$')
        plt.plot(tm[0:i],Q2[0:i],'k--',label=r'$Q_2$')
        plt.ylabel('Heaters')
        plt.xlabel('Time (sec)')
        plt.legend(loc='best')
        plt.draw()
        plt.pause(0.05)
    
    #Turn off heater
    a.Q1(0)
    a.Q2(0)
    
    #Turn off LED
    a.LED(0)
    
    #Save plot file
    plt.savefig('testEBmodel.png')
        
    

    
# Allow user to end loop with Ctrl-C           
except KeyboardInterrupt:
    # Disconnect from Arduino
    a.Q1(0)
    a.Q2(0)
    a.LED(0)
    print('Shutting down')
    a.close()
    plt.savefig('testEBmodel.png')

# Make sure serial connection still closes when there's an error
except:           
    # Disconnect from Arduino
    a.Q1(0)
    a.Q2(0)
    a.LED(0)
    print('Error: Shutting down')
    a.close()
    plt.savefig('testEBmodel.png')
    raise