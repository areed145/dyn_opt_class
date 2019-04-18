#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 15:14:40 2019

@author: areed145
"""

import pandas as pd
from gekko import GEKKO
import math
import numpy as np
import matplotlib.pyplot as plt

rmt = True

def plot_results():
    plt.figure(1)
    
    plt.subplot(6,1,1)
    plt.plot(telapsed_ / 3600, md_, 'r-', label='Measured depth')
    plt.plot(telapsed_ / 3600, tvd_, 'b-', label='True vertical depth')
    plt.ylabel('Depth (m)')
    plt.legend(loc='best')
    
    plt.subplot(6,1,2)
    plt.plot(telapsed_ / 3600, Pmp_, 'k-', label='Mud pump discharge pressure')
    plt.plot(telapsed_ / 3600, Pdh_, 'r-', label='Drill bit pressure')
    plt.plot(telapsed_ / 3600, Pf_, 'g:', label='Formation pressure')
    plt.plot(telapsed_ / 3600, Pchoke_, 'b-', label='Choke valve inlet pressure')
    plt.ylabel('Pressure (bar)')
    plt.legend(loc='best')
    
    plt.subplot(6,1,3)
    plt.plot(telapsed_ / 3600, Qmp_, 'k-', label='Mud pump flowrate')
    plt.plot(telapsed_ / 3600, Qbit_, 'r-', label='Drill bit flowrate')
    plt.plot(telapsed_ / 3600, Qres_, 'g:', label='Formation flowrate')
    plt.plot(telapsed_ / 3600, Qchoke_, 'b-', label='Choke flowrate')
    plt.plot(telapsed_ / 3600, Qbp_, 'y-', label='Back-pressure pump flowrate')
    plt.ylabel('Flow (m3/min)')
    plt.legend(loc='best')
    
    plt.subplot(6,1,4)
    plt.plot(telapsed_ / 3600, rop_, 'k-', label='ROP')
    plt.ylabel('ROP (m/hr)')
    plt.legend(loc='best')
    
    plt.subplot(6,1,5)
    plt.plot(telapsed_ / 3600, hpit_, 'k-', label='Pit Level')
    plt.ylabel('Height (m)')
    plt.legend(loc='best')
    
    plt.subplot(6,1,6)
    plt.plot(telapsed_ / 3600, cv_, 'b-', label='Choke Position')
    plt.ylabel('%')
    plt.legend(loc='best')
    
    plt.xlabel('Time (hr)')
    plt.show()

def load_reservoir():
    '''
    Returns:
        res - reservoir model
    '''
    
    res = pd.read_excel(r'example well pressure and flow function of depth and time.xlsx', skiprows=10)
    res.drop(0, inplace=True)
    res.sort_values(by=['MD Measured Depth (ft)'], inplace=True)

    res['MD Measured Depth (m)']  = res['MD Measured Depth (ft)'] * 0.3048
    res['TVD True Vertical Depth(m)'] = res['TVD True Vertical Depth(ft)']  * 0.3048
    res['ROP (m/sec)'] = res['ROP (m/sec)'] * 10
    res['Permeability_k (m2)'] = res['Permeability_k (m2)'] * 1e5

    return res

def calc_weighted_PIh(md_, k_, mu):
    '''
    '''
    
    PIh = 0
    for idx, i in enumerate(md_):
        try:
            PIh += k_[idx] * 1e3 / mu * (md_[idx+1] - md_[idx])
        except:
            pass
    return PIh

def get_reservoir(res, md):
    '''
    Arguments:
        res   - reservoir model
        md    - Measured depth (m)
    Returns:
        tvd   - true vertical depth (m)
        rop   - rate of penetration (m/s)
        pres    - formation pressure (bar)
        k     - permeability (m2)
        el    - effective length (m)
    '''
    
    try:
        res_rel = res[res['MD Measured Depth (m)'] >= md]
        tvd = res_rel['TVD True Vertical Depth(m)'].values[0]
        rop = res_rel['ROP (m/sec)'].values[0]
        pres = res_rel['Pressure (formation bar)'].values[0]
        k = res_rel['Permeability_k (m2)'].values[0]
        el = 1
    except:
        tvd = res_rel['TVD True Vertical Depth(m)'].values[-1]
        rop = res_rel['ROP (m/sec)'].values[-1]
        pres = res_rel['Pressure (formation bar)'].values[-1]
        k = res_rel['Permeability_k (m2)'].values[-1]
        el = 1
        
    if md > 1000:
        rop = rop / 3 

    return tvd, rop, pres, k, el

res = load_reservoir()

D_max = 2300                                   # m
t_delta = 60.0 * 60.0 * 2                      # seconds
st = 180.0                                    # seconds
mu = 0.042                                    # kinematic viscosity (kg/m*s)

t_delta_ = np.array([t_delta])
st_ = np.array([st])
t_elapsed_ = np.array([0])
md_ = np.array([0.0])

TVD, ROP, PRES, K, EL = get_reservoir(res, 0)

tvd_ = np.array([TVD])
rop_ = np.array([ROP])
P_res_ = np.array([PRES])
k_ = np.array([K])
el_ = np.array([EL])

PIH = calc_weighted_PIh(md_, k_, mu)
PIh_ = np.array([PIH])

cv_ = np.array([20.0])
frac_mud_ = np.array([0.5])
rho_mud_ = np.array([1800.0])
rho_water_ = np.array([1000.0])
rho_pit_ = np.array([1400.0])
rho_res_ = np.array([500.0])

h_pit_ = np.array([4.0])
A_pit_ = np.array([100.0])

Q_mp_ = np.array([1.0])
P_mp_ = np.array([0.0])
Q_bit_ = np.array([0.0])
Q_res_ = np.array([0.0])
Q_bp_ = np.array([0.4])
P_choke_ = np.array([0.0])

rho_C_ = np.array([1240.0])
Q_choke_ = np.array([0.0])

Q_mud_ = np.array([0.0])
Q_water_ = np.array([0.0])

V_mud_ = np.array([100])
V_water_ = np.array([100])

drilling = True

mhe = GEKKO(remote=rmt)
mpc = GEKKO(remote=rmt)

# create 2 models (MHE and MPC) in one loop
for m in [mhe, mpc]:
    
    ######## OTHER CONSTANTS ########
    m.P_atm = m.Const(1)                                   # Atmospheric pressure outside choke [bar]
    m.g = m.Const(9.81)                                   # Gravitational constant (kg-m/s^2)
    m.r_do = 2.5 * 0.0254                                 # drillstring outer radius (m) (5" diameter)
    m.r_ci = 4.3125 * 0.0254                              # annulus inner radius (m) (8 5/8" diameter)
    m.r_di = 2.0 * 0.0254                                 # drillstring inner radius (m) (4" diameter)
    
    ######## DRILLING ########
    m.md  = m.Param(0)                              # Measured depth (m)
    m.tvd = m.Param(tvd_[-1])                             # Total vertical depth of well [m]
    m.rop = m.Param(rop_[-1])                             # rate of penetration (m/sec)
    
    ######## MUDPIT ########
    m.A_pit = m.Const(value=A_pit_[-1])                               # Mud pit area (m2)
    m.V_mud = m.MV(value=V_mud_[-1], lb=0.0, ub=1000.0)
    m.V_water = m.MV(value=V_water_[-1], lb=0.0, ub=1000.0)
    m.rho_mud = m.MV(value=rho_mud_[-1], lb=1600.0, ub=2000.0)
    m.rho_water = m.Const(value=rho_water_[-1])
    
    m.V_pit = m.Var(value=(m.V_mud + m.V_water))
    m.h_pit = m.CV(value=(m.V_pit / m.A_pit))                               # CV - Mud pit level
    m.frac_mud = m.Var(value=(m.V_mud / m.V_pit))
    m.M_pit = m.Var(value=(m.V_mud * m.rho_mud + m.V_water * m.rho_water))
    m.rho_pit = m.Var(value=(m.M_pit / m.V_pit))                             # CV - Mud density
    
    ######## MUDPUMP ########
    m.P_mp = m.Var(value=P_mp_[-1])
    m.Q_mp = m.MV(Q_mp_[-1])                             # Mud pump flow rate[m^3/min]
    m.rho_mp = m.Var(m.rho_pit)
    
    m.Q_mud = m.Var(m.Q_mp * m.frac_mud)
    m.Q_water = m.Var(m.Q_mp - m.Q_mud)
    
    ######## DRILLSTRING ########
    m.A_ds = m.Const(math.pi * m.r_di**2)                     # drillstring inner area , m^2
    m.beta_ds  = m.Param(90000)                             # Bulk Modulus of Mud in Drill String [bar]
    m.F_ds     = m.Param(80)                                # Friction Factor in the drill string [bar*s^2/m^6]
    m.rho_ds   = m.Var(m.rho_pit)                         # Mud Density in the Drill String [kg/m^3]
    m.V_ds = m.Intermediate(m.A_ds * m.md)                        # Volume of Drill String [m^3]
    m.Q_ds = m.Var(m.Q_mp)
    
    ######## BIT ########
    m.P_bit_init = (m.rho_ds * (-m.F_ds / 3600) * m.md * (m.Q_mp**2) + m.rho_ds * m.g * m.tvd) * 1e-5
    m.P_bit = m.Var(m.P_bit_init)                           # Bit pressure [bar]
    m.Q_bit = m.Var(m.Q_mp)
    m.rho_bit = m.Var(m.rho_ds)
    
    ######## RESERVOIR ########
    m.P_res  = m.Param(P_res_[-1])                              # Formation Pressure (bar)
    m.k   = m.Param(k_[-1])                               # Permeability
    m.el  = m.Param(el_[-1])                              # Effective Length (m)
    m.PIh = m.Param(PIh_[-1])
    #m.Q_res_init = m.PIh * (m.P_res - m.P_bit) / math.log(10.0 / m.r_ci)
    m.Q_res = m.Var(0, lb=0) # Reservoir fluid influx flow rate [m^3/min]
    m.rho_res = m.Var(value=rho_res_[-1])
    #m.Q_resm_init = m.PIh * (m.P_bit - m.P_res) / math.log(10.0 / m.r_ci)
    m.Q_resm = m.Var(0, lb=0) # Reservoir fluid influx flow rate [m^3/min]
    
    ######## ANNULUS ########
    m.A_ann = m.Const(math.pi * (m.r_ci**2 - m.r_do**2))         # annulus flow area , m^2
    m.beta_ann  = m.Param(50000)
    m.F_ann     = m.Param(330)                               # Friction Factor in the Annulus [bar*s^2/m^6]
    m.V_ann = m.Intermediate(m.A_ann * m.md)                        # Volume of Annulus [m^3]
    m.rho_ann = m.Var(m.rho_ds)
    
    #Equations
    # Bit pressure [bar]
    m.Equation(m.P_bit == (m.rho_ann * (m.F_ann / 3600) * m.md * (m.Q_mp**2) + m.rho_ann * m.g * m.tvd) * 1e-5)
    #m.Equation(m.P_bit == (m.rho_ds * (-m.F_ds / 3600) * m.md * (m.Q_mp**2) + m.rho_ds * m.g * m.tvd) * 1e-5)

    # Flow Rate from Reservoir
    m.Equation(m.Q_res == m.PIh * (m.P_res - m.P_bit))
    m.Equation(m.Q_resm == m.PIh * (m.P_bit - m.P_res))

    # Mud pump discharge (Equation 5.1)
    m.Equation(m.P_mp.dt() == (m.beta_ds / m.V_ds) * (m.Q_mp - m.Q_mp - 60 * m.rop * m.A_ds) )

    # Flow through drill bit (Equation 5.3)
    m.M = m.Intermediate((m.rho_ds + m.rho_ann) * m.md)                         # Total Mud Density per length [kg/m^4]
    m.Equation(m.Q_bit.dt() == (1e+5 / m.M) * (m.P_mp - m.P_bit - m.F_ds / 3600 * (m.Q_bit**2) + m.rho_ds * m.g * m.tvd / 1e+5))
    
    
    # Mudpit equations
    m.Equation(m.Q_mud == m.Q_mp * m.frac_mud)
    m.Equation(m.Q_water == m.Q_mp - m.Q_mud)    
    #m.Equation(m.V_pit == m.V_water + m.V_mud)
    #m.Equation(m.h_pit == m.V_pit / m.A_pit)
    #m.Equation(m.frac_mud == m.V_mud / m.V_pit)
    #m.Equation(m.M_pit == (m.V_mud * m.rho_mud + m.V_water * m.rho_water))
    #m.Equation(m.rho_pit == (m.V_water * m.rho_water + m.V_mud * m.rho_mud) / m.V_pit)
    
    #m.Equation(m.rho_ann.dt() == (m.rho_bit * (m.Q_bit - m.Q_resm) + m.rho_res * m.Q_res) / 60)                # Mud Density in the Drill String Annulus [kg/m^3]
    m.Equation(m.h_pit.dt() == (1 / m.A_pit) * (m.Q_res - m.Q_resm - m.Q_mp) / 60)
    #m.Equation(m.frac_mud.dt() == (m.Q_mud / m.Q_mp) / 60)
    #m.Equation(m.M_pit.dt() == (m.Q_water * m.rho_water + m.Q_mud * m.rho_mud - m.Q_mp * m.rho_pit) / 60)
    m.Equation(m.V_pit.dt() == (m.Q_res - m.Q_resm - m.Q_mp) / 60)
    #m.Equation(m.rho_pit.dt() == ((m.Q_water * m.rho_water + m.Q_mud * m.rho_mud) / m.Q_mp) / 60)
    
###################################################
# Configure MHE
# 120 second time horizon, steps of 3 sec
ntm = 40
mhe.time = np.linspace(0,ntm*3,ntm+1)

# Measured inputs

mhe.Q_mp.STATUS = 0  # not changed by mhe
mhe.Q_mp.FSTATUS = 1 # measured

mhe.V_mud.STATUS = 0  # not changed by mhe
mhe.V_mud.FSTATUS = 1 # measured

mhe.V_water.STATUS = 0  # not changed by mhe
mhe.V_water.FSTATUS = 1 # measured

mhe.h_pit.FSTATUS = 1    # receive measurement
meas_gap = 2.0
mhe.h_pit.MEAS_GAP = meas_gap # for CV_TYPE = 1

mhe.options.IMODE     = 5 # MHE Mode
mhe.options.EV_TYPE   = 1 # Estimator Objective type
mhe.options.NODES     = 3 # Collocation nodes
mhe.options.ICD_CALC  = 1 # Calculate initial conditions
mhe.options.SOLVER    = 3 # IPOPT
mhe.options.COLDSTART = 0 # COLDSTART on first cycle

mhe.solve()

###################################################
# Configure MPC
# Control horizon, non-uniform time steps
nt = int(t_delta / st) + 1                      # simulation time points
mpc.time = np.linspace(0, t_delta, nt)

# update parameters from mhe
mpc.K1.FSTATUS = 1
mpc.K2.FSTATUS = 1
mpc.K3.FSTATUS = 1
mpc.tau12.FSTATUS = 1
mpc.tau3.FSTATUS = 1

# Measured inputs
mpc.Q_mp.STATUS = 1  # manipulated
mpc.Q_mp.FSTATUS = 0 # not measured
mpc.Q_mp.DMAX = 20.0
mpc.Q_mp.DCOST = 0.1

mpc.V_mud.STATUS = 1  # manipulated
mpc.V_mud.FSTATUS = 0 # not measured
mpc.V_mud.DMAX = 30.0
mpc.V_mud.DCOST = 0.1

mpc.V_water.STATUS = 1  # manipulated
mpc.V_water.FSTATUS = 0 # not measured
mpc.V_water.DMAX = 30.0
mpc.V_water.DCOST = 0.1

mpc.h_pit.STATUS = 1     # drive to setpoint
mpc.h_pit.FSTATUS = 1    # receive measurement
mpc.h_pit.TAU = 40       # response speed (time constant)
mpc.h_pit.TR_INIT = 1    # reference trajectory
mpc.h_pit.TR_OPEN = 0

# Global Options
mpc.options.IMODE   = 6 # MPC Mode
mpc.options.CV_TYPE = 1 # Controller Objective type
mpc.options.NODES   = 3 # Collocation nodes
mpc.options.SOLVER  = 3 # IPOPT
mpc.options.COLDSTART = 1 # COLDSTART on first cycle

print('-[Process Initialization]-----------------------')
print('Time elapsed            =', t_elapsed_[-1], 's')
print('Measured depth          =', md_[-1], 'm')
print('Mud pit level           =', h_pit_[-1], 'm')
#print('Mud density             =', rho_P_[-1], 'kg/m3')
print('Pressure at choke valve =', P_choke_[-1], 'bar')
print('Flow rate through choke =', Q_choke_[-1], 'm3/min')
print('Downhole pressure       =', P_res_[-1], 'bar')

# main loop
i = 0
while drilling:

    
    PIh = calc_weighted_PIh(md_, k_, mu)
    
    if md_[-1] > 1000:
        cv_[-1] = 10.0
        
    if i==10:
            mhe.K1.STATUS = 1
            mhe.K2.STATUS = 1
            
    # Get the relevant reservoir properties
    TVD, ROP, PRES, K, EL = get_reservoir(res, md_[-1])
    tvd_ = np.append(tvd_, TVD)
    rop_ = np.append(rop_, ROP)
    P_res_ = np.append(P_res_, PRES)
    k_ = np.append(k_, K)
    el_ = np.append(el_, EL)
    
    PIH = calc_weighted_PIh(md_, k_, mu)
    PIh_ = np.append(PIh_, PIH)
    
    # Insert measurements to MHE
    mhe.tvd.MEAS = TVD
    mhe.rop.MEAS = ROP
    mhe.P_res.MEAS = PRES
    mhe.k.MEAS = K
    mhe.PIh.MEAS = PIH
    mhe.el.MEAS = EL

    # Update model parameters with MHE
    try:
        mhe.solve(disp=False)

        # Insert updated parameters to MPC
        mpc.K1.MEAS    = mhe.K1.NEWVAL   
        mpc.K2.MEAS    = mhe.K2.NEWVAL   
        mpc.K3.MEAS    = mhe.K3.NEWVAL   
        mpc.tau12.MEAS = mhe.tau12.NEWVAL
        mpc.tau3.MEAS  = mhe.tau3.NEWVAL
    except:
        print('MHE solution failed, using prior values')

    # Insert temperature measurement for MPC
    mpc.TC1.MEAS   = mhe.TC1.MODEL  # or T1m[i] 
    mpc.TC2.MEAS   = mhe.TC2.MODEL  # or T2m[i]

    # Adjust setpoints
    db1 = 1.0 # dead-band
    mpc.TC1.SPHI = T1sp[i] + db1
    mpc.TC1.SPLO = T1sp[i] - db1

    db2 = 0.5
    mpc.TC2.SPHI = T2sp[i] + db2
    mpc.TC2.SPLO = T2sp[i] - db2

    # Adjust Heaters with MPC
    try:
        # Solve MPC
        mpc.solve(disp=False)

        # Retrieve new values
        Q1s[i+1]  = mpc.Q1.NEWVAL
        Q2s[i+1]  = mpc.Q2.NEWVAL
        # get additional solution information
        with open(m.path+'//results.json') as f:
            results = json.load(f)
    except:
        print('MPC solution failed, turn off heaters')
        Q1s[i+1]  = 0.0
        Q2s[i+1]  = 0.0

    # Write new heater values (0-100)
    a.Q1(Q1s[i])
    a.Q2(Q2s[i])

    print('-[Process]--------------------------------------')
    print('Time elapsed            =', t_elapsed_[-1] / 3600, 'hr')
    print('Measured depth          =', md_[-1], 'm')
    print('Mud pit level           =', h_pit_[-1], 'm')
    print('Mud density             =', rho_P_[-1], 'kg/m3')
    print('Pressure at choke valve =', P_choke_[-1], 'bar')
    print('Flow rate through choke =', Q_choke_[-1], 'm3/min')
    print('Downhole pressure       =', P_res_[-1], 'bar')

    print('-[Reservoir]------------------------------------')
    print('TVD                     =', tvd_[-1], 'm')
    print('Rate of penetration     =', rop_[-1] * 3600, 'm/hr')
    print('Pressure of formation   =', P_res_[-1], 'bar')
    print('Permeability            =', k_[-1], 'm/s')
    print('Effective length        =', el_[-1], 'm')

    print('-[MHE]------------------------------------------')
    print('Pit area                =', A_pit_[-1], 'm2')

    print('-[MPC]------------------------------------------')
    print('Choke valve position    =', cv_[-1], '%')
    print('Back-pressure pump flow =', Q_bp_[-1], 'm3/min')


    # Go back to the top
    print('------------------------------------------------\n')

    if md > D_max:
        drilling = False
        
    i +=1

plot_results()
