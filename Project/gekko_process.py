# -*- coding: utf-8 -*-
"""
Dynamic Optimization - Group 17 Project
Managed Pressure Drilling Model in GEKKO


"""

#%%Import packages
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

import drillingmodel as dm

depth = 600

print(dm.drillstring(depth))


#%% Specify controller model
tf = 100        # control horizon
npts = 10       # condensation points

c = GEKKO(remote=rmt)
c.time = np.linspace(0, tf, npts+1)  # Controller horizon




#%% Display results
#plt.figure(3)
#plt.subplot(3,1,1)
#plt.plot(m.time, Pbit.VALUE, 'r', label='BHP(Pa)')
#plt.plot(m.time, Pp.VALUE, 'b', label='Pump Pressure(Pa)')
#plt.legend(loc='best')
#plt.subplot(3,1,2)
#plt.plot(m.time, Pc.VALUE, 'g', label ='Choke Pressure(Pa)')
#plt.legend(loc='best')
#plt.subplot(3,1,3)
#plt.plot(m.time, np.array(Qbit.VALUE)*60000, 'k', label ='Qbit (l/min)')
#plt.plot(m.time, np.array(Qchoke.VALUE)*60000, 'r', label ='Qchoke (l/min)')
#plt.legend(loc='best')
#plt.show()

