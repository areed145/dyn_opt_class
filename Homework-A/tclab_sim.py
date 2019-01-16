#import tclab
import numpy as np
import time
import matplotlib.pyplot as plt
from gekko import GEKKO

# Set up the heater model in GEKKO
m = GEKKO(remote=False)

# Set up the model time array (10 minutes, 1 second intervals)
m.time = np.linspace(0, 60*10, 60*10+1)

# Constants
t0 = 23.0 + 273.15    # Initial temperature
tAmb = 23.0 + 273.15  # Ambient temperature
area = 12.0e-4   # Area
mass = 4.0e-3
Cp = 500.0
AmCp = 6.0e-4    # A/mCp from notes
sig = 5.67e-8    # Stefsn-Boltzmann constant
eps = 0.9        # Emissivity 
alpha = 0.01         # Heater scaling factor

U = 10         # Initial estimate of HT coefficient

# Parameters:
q0 = np.zeros(m.time.shape)
q0[:60]     = 0.0
q0[60:240]  = 100.0
q0[240:300] = 0.0
q0[300:360] = 100.0
q0[360:]    = 0.0
Q = m.Param(value=q0)

# Variables:
T = m.Var(value=t0)

# Equation
m.Equation(mass*Cp*T.dt()== \
           U*area*(tAmb-T) \
           + eps*sig*area*(tAmb**4-T**4) \
           + alpha*Q)

# Simulation mode
m.options.IMODE = 4

# Solve simulation model
m.solve()

# Show results
plt.figure(1)

plt.subplot(2,1,1)
plt.plot(m.time/60.0,np.array(T.value)-273.15,'r-')
plt.ylabel('Temperature (degC)')
plt.legend(['Step Test (0-100% heater)'])

plt.subplot(2,1,2)
plt.plot(m.time/60.0,np.array(Q.value),'b-')
plt.ylabel('Heater (%)')
plt.xlabel('Time (min)')

plt.show()
