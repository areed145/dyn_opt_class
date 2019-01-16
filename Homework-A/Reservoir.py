#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HW A - Reservoir Simulation
Due on Sat Jan 12
@author: Rob Hawkins

Simulation Basis:

Outflow River Rates (km3/yr) with height in meters
 Vflow_out1 = 0.030 sqrt(h1) 
 Vflow_out2 = 0.015 sqrt(h2) 
 Vflow_out3 = 0.060 sqrt(h3)
 Vflow_out4 = 0

Evaporation Rates (km3/yr)
 Vevap = 0.5e-5 * Area, for salt water (Great Salt Lake)
 Vevap = 1e-5 * Area, for fresh water (all others)

Inflow Rates (km3/yr)
 Vflow_in1 = 0.13 (June-Feb), 0.21 (Mar-May)
 Vflow_in2 = Vflow_out1
 Vflow_in3 = Vflow_out2
 Vflow_in4 = Vflow_out3

Usage Requirements (km3/yr)
 Vuse1 = 0.03
 Vuse2 = 0.05
 Vuse3 = 0.02
 Vuse4 = 0.00

Area of Reservoir / Lake (km2)
 A1 = 13.4
 A2 = 12.0
 A3 = 384.5
 A4 = 4400

Initial Volume of Reservoir / Lake (km3)
 V1 = 0.26
 V2 = 0.18
 V3 = 0.68
 V4 = 22.0
"""

# Import necessary libraries
from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt

# Problem Constants
evap = np.array([1.0e-5, 1.0e-5, 1.0e-5, 0.5e-5])  # Evaporation
areas = np.array([13.4, 12.0, 384.5, 4400.0])      # Reservoir Areas (km2)
vEvap = evap * areas                               # Evaporation rates (km3/yr)

vUse = np.array([0.03, 0.05, 0.02, 0.0])           # Water Usage (km3/yr):

k = np.array([0.03, 0.015, 0.06, 0.0])             # Outflow rate

# Initial conditions
v0 = np.array([0.26, 0.18, 0.68, 22.0])  # Reservoir volumes (km3)
h0 = 1000 * v0 / areas                   # Reservoir heights (m)
vOut0 = k * np.sqrt(h0)                  # Outflow rates (km3/yr)

# Daily timeline
runoff = np.zeros(365)            # Inflow to Jordanelle
runoff[:59]    = 0.13  # Normal inflow (Jan-Feb)
runoff[59:151] = 0.21  # Spring run-off (Mar-May)
runoff[151:]   = 0.13  # Normal inflow (June-Dec)

#Monthly Timeline
#runoff = np.zeros(12)            # Inflow to Jordanelle
#runoff[:3]  = 0.13  # Normal inflow (Jan-Feb)
#runoff[3:6] = 0.21  # Spring run-off (Mar-May)
#runoff[6:]  = 0.13  # Normal inflow (June-Dec)


# Create Simulation model
r = GEKKO(remote=False)
# Time array (daily)
r.time = np.linspace(0, 1, 365)

# Time array (monthly)
#r.time = np.linspace(0, 1, 12)

vIn = [0.0, 0.0, 0.0, 0.0]  #Inlet flows

# Define model parameters
vIn[0] = r.Param(value=runoff)    # Inflow to Jordanelle Reservoir

# Define model Variables
V = [r.Var(value=x) for x in v0]
H = [r.Var(value=x) for x in h0]
vOut = [r.Var(value=x) for x in vOut0]

# Define model Intermediates
vIn[1:4] = [r.Intermediate(vOut[i]) for i in range(3)]

# Equations
r.Equations([1000*V[i] == H[i] * areas[i] for i in range(4)])
r.Equations([vOut[i]**2 == k[i]**2 * H[i] for i in range(4)])
r.Equations([V[i].dt() == vIn[i] - vOut[i] - vEvap[i] - vUse[i] \
             for i in range(4)])
 
# Set Solver options
r.options.IMODE = 4   # Dynamic simulation
r.options.SOLVER = 1  # change solver (1=APOPT,3=IPOPT)

# Simulate the dayamics
r.solve(disp=False)
 
# Display the results

plt.figure(1)

plt.subplot(2,1,1)
plt.plot(r.time, H[0].value, 'k-')
plt.plot(r.time, H[1].value, 'r-')
plt.ylabel('Height (m)')
plt.xlabel('Time (yr)')
plt.legend(['Jordanelle','Deer Creek'])

plt.subplot(2,1,2)
plt.plot(r.time, H[2].value, 'g-')
plt.plot(r.time, H[3].value, 'b-')
plt.ylabel('Height (m)')
plt.xlabel('Day of the year')
plt.legend(['Utah Lake','Great Salt Lake'])

plt.figure(2)

#plt.subplot(1,1,1)
plt.plot(r.time, vIn[0].value, 'k-')
plt.plot(r.time, vOut[0].value, 'r-')
plt.plot(r.time, vOut[1].value, 'g-')
plt.plot(r.time, vOut[2].value, 'b-')
plt.ylabel('Flow (km3/yr)')
plt.xlabel('Time (yr)')
plt.legend(['Supply','Upper Provo','Lower Provo','Jordan'])

plt.show()
