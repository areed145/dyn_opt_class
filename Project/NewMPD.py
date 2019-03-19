#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
import reservoirmodel as rm
from gekko import GEKKO
import matplotlib.pyplot as plt


d = GEKKO()

# Model Constants
Patm = d.Const(1)              # Atmospheric pressure outside choke [bar]    
Md = d.Const(2500)            # Lumped Density per length of Mud in Drill String [kg/m^4 * 1e5]
Ma = d.Const(800)             # Lumped Density per length of Mud in Annulus [kg/m^4 * 1e5]
g = d.Const(9.81)               # Gravitational constant (kg-m/s^2)
r_di = 2.0 * 0.0254             # drillstring inner radius (m) (4" diameter)
r_do = 2.5 * 0.0254             # drillstring outer radius (m) (5" diameter)
r_ci = 4.3125 * 0.0254          # annulus inner radius (m) (8 5/8" diameter)
Ad = d.Const(math.pi*r_di**2)  # drillstring inner area , m^2
Aa = d.Const(math.pi*(r_ci**2 - r_do**2))  # annulus flow area , m^2
Ah = d.Const(math.pi*r_ci**2)  # borehole cross area, m^2

Mdepth = 4200


# Parameters
Kc     = d.Param(0.4)        # Valve Coefficient 
betaD = d.Param(90000)          # Bulk Modulus of Mud in Drill String [bar]
betaA = d.Param(50000)    
Fd    = d.Param(80)             # Friction Factor in the drill string [bar*s^2/m^6]
Fa    = d.Param(330)            # Friction Factor in the Annulus [bar*s^2/m^6]
rhoD   = d.Param(1240)          # Mud Density in the Drill String [kg/m^3]        
rhoA   = d.FV(1290,lb=rhoD)     # Mud Density in the Drill String Annulus [kg/m^3]    
Qpump = d.Param(2.0)               # Flow Rate through Pump [m^3/min] 
Zc = d.Param(30)                # Choke Valve Opening from 0-100 [%]
Qres = d.Param(0)               # reservoir gas influx flow rate [m^3/min]
Qloss = d.Param(0)
Qback = d.Param(0)              # back pressure pump flow rate [m^3/min]
ROP = d.Param(20)                # rate of penetration (m/min)

# Variables
Qbit = d.Var(Qpump - ROP*Ad,lb=0)        # Flow Rate through Bit [m^3/min]

Pp = d.Var(38)                # Pressure at Pump [bar]
Pc = d.Var(2,lb=Patm)          # Pressure at Choke Valve [bar]
TVD = d.Param(3950)             # Total vertical depth of well [m]
Pbit_init = Pc + (rhoA*(Fa/3600)*Mdepth*(Qbit**2) + rhoA*g*TVD)*1e-5
Pbit = d.Var(Pbit_init)             # Bit pressure [bar]
MD = d.Var(value=Mdepth)

Qchoke_init = Kc * Zc * d.sqrt(rhoA*(Pc-Patm)*1e-5) 
Qchoke = d.Var(Qchoke_init,lb=0)        # Flow Rate through Choke [m^3/min]


# Intermediates
M = d.Intermediate(Md+Ma)    # Total Mud Density per length [kg/m^4]
Va = d.Intermediate(Aa*MD)  # Volume of Annulus [m^3]
Vd = d.Intermediate(Ad*MD)  # Volume of Drill String [m^3]
# Bit pressure [bar]
#Pbit_init = d.Intermediate(Pc + (rhoA*(Fa/3600)*MD*(Qbit**2) + rhoA*g*TVD)*1e-5)
#Pbit = d.Var(Pc + (rhoA*(Fa/3600)*MD*(Qbit**2) + rhoA*g*TVD)*1e-5)
# Flow Rate through Choke Valve [m^3/min] 
#Qchoke_init = d.Intermediate(Kc * Zc * d.sqrt(rhoA*(Pc-Patm)*1e-5))

#Equations

#d.Equation(Pbit == Pbit_init)
d.Equation(Pbit == Pc + (rhoA*(Fa/3600)*MD*(Qbit**2) + rhoA*g*TVD)*1e-5)
d.Equation(Qchoke == Kc * Zc * d.sqrt(rhoA*(Pc-Patm)*1e-5))
d.Equation(Qres == rm.reservoir_flow(Pbit, Ah, MD))
# Equation 5.1
d.Equation( Pp.dt() == (betaD/Vd) * (Qpump - Qbit - ROP*Ad) )

# Equation 5.2
#d.Equation( Pc.dt() == (betaA/Va) * (Qres + Qbit + Qback - Qchoke - Qloss \
#                       - ROP*Aa))
d.Equation( Pc.dt() == (betaA/Va) * (Qres + Qbit + Qback - Qchoke - ROP*Aa))

# Equation 5.3
d.Equation( Qbit.dt() == (1e+5/M) * (Pp - Pbit - Fd/3600*(Qbit**2) \
                    + rhoD*g*TVD/1e+5 ) )


# Drilling rate from reservior simulation
d.Equation( MD.dt() == ROP + 0*MD )


# Options
#d.options.solver = 1

# Calculate starting conditions
#d.options.imode = 3
#d.solve()

# Print solution
print('------------------------------------------------')
print('Pressure at Pump [bar]', Pp.value)
print('Pressure at Choke Valve [bar]', Pc.value)
print('Bit pressure [bar]', Pbit.value)
print('Flow Rate through Bit [m^3/min] ', Qbit.value)
print('Flow Rate through Choke [m^3/min]', Qchoke.value)
print('Bit Height [m]', TVD.value)
print('Measured Depth [m] ', MD.value)
print('------------------------------------------------')
print('')

# Simulation time
tf = 1200.0  # final time (sec)
st = 2.0   # simulation time intervals
nt = int(tf/st)+1 # simulation time points
d.time = np.linspace(0,tf,nt)/60.0


# Configure choke valve step
#zc = np.ones(nt)*50.0
#zc[50:] = 25.0  # change to 25% open
#zc[100:] = 20.0 # change to 20% open
#Zc.value = zc

# Configure pump ramp
#q_pump = np.ones(nt)*5.0
## ramp up
#for i in range(5,len(q_pump)):
#    if q_pump[i-1]<=2.95:
#        q_pump[i] = q_pump[i-1]+0.05
#    else:
#       q_pump[i] = 3.0
#Qpump.value = q_pump

d.options.imode = 4 # dynamic simulation
d.solve()

plt.figure(1)

plt.subplot(4,1,1)
plt.plot(d.time,Qpump,'b-',label='Mud Pump Flow')
plt.plot(d.time,Qbit,'g:',label='Bit Mud Flow')
plt.plot(d.time,Qchoke,'r--',label='Choke Mud Flow')
plt.ylabel(r'Flow ($m^3/min$)')
plt.legend(loc='best')

plt.subplot(4,1,2)
plt.plot(d.time,Zc,'k-',label='Choke Opening (%)')
plt.ylabel('Choke (%)')
plt.legend(loc='best')

plt.subplot(4,1,3)
plt.plot(d.time,Pbit,'r-',label='Bit Pressure (bar)')
plt.ylabel('Press (bar)')
plt.legend(loc='best')

plt.subplot(4,1,4)
plt.plot(d.time,Pp,'r:',label='Pump Pressure (bar)')
plt.plot(d.time,Pc,'b--',label='Choke Pressure (bar)')
plt.legend(loc='best')
plt.ylabel('Press (bar)')
plt.xlabel('Time (min)')
plt.show()


