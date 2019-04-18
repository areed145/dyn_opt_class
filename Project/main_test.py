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
        pf    - formation pressure (bar)
        k     - permeability (m2)
        el    - effective length (m)
    '''
    
    try:
        res_rel = res[res['MD Measured Depth (m)'] >= md]
        tvd = res_rel['TVD True Vertical Depth(m)'].values[0]
        rop = res_rel['ROP (m/sec)'].values[0]
        Pf = res_rel['Pressure (formation bar)'].values[0]
        k = res_rel['Permeability_k (m2)'].values[0]
        el = 1
    except:
        tvd = res_rel['TVD True Vertical Depth(m)'].values[-1]
        rop = res_rel['ROP (m/sec)'].values[-1]
        Pf = res_rel['Pressure (formation bar)'].values[-1]
        k = res_rel['Permeability_k (m2)'].values[-1]
        el = 1
        
    if md > 1000:
        rop = rop / 3 

    return tvd, rop, Pf, k, el

res = load_reservoir()

D_max = 2300                                   # m
t_delta = 60.0 * 60.0 * 2                      # seconds
st = 180.0                                    # seconds
mu = 0.042                                    # kinematic viscosity (kg/m*s)

t_delta_ = np.array([t_delta])
st_ = np.array([st])
t_elapsed_ = np.array([0])
md_ = np.array([0.0])

TVD, ROP, PF, K, EL = get_reservoir(res, 0)

tvd_ = np.array([TVD])
rop_ = np.array([ROP])
P_res_ = np.array([PF])
k_ = np.array([K])
el_ = np.array([EL])

PIH = calc_weighted_PIh(md_, k_, mu)
PIh_ = np.array([PIH])

cv_ = np.array([20.0])
frac_mud_ = np.array([0.5])
rho_mud_ = np.array([1500.0])
rho_water_ = np.array([1000.0])
rho_pit_ = np.array([1250.0])
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

drilling = True

mhe = GEKKO(remote=rmt)
mpc = GEKKO(remote=rmt)

# create 2 models (MHE and MPC) in one loop
for m in [mhe, mpc]:
    
    ######## DRILLSTRING ########
    
    # Model Constants
    m.P_atm = m.Const(1)                                   # Atmospheric pressure outside choke [bar]
    m.M_ds = m.Const(2500)                                  # Lumped Density per length of Mud in Drill String [kg/m^4 * 1e5]
    m.M_ann = m.Const(800)                                   # Lumped Density per length of Mud in Annulus [kg/m^4 * 1e5]
    m.g = m.Const(9.81)                                   # Gravitational constant (kg-m/s^2)
    m.r_di = 2.0 * 0.0254                                 # drillstring inner radius (m) (4" diameter)
    m.r_do = 2.5 * 0.0254                                 # drillstring outer radius (m) (5" diameter)
    m.r_ci = 4.3125 * 0.0254                              # annulus inner radius (m) (8 5/8" diameter)
    m.Kc     = m.Param(0.4)                               # Valve Coefficient
    m.A_ds = m.Const(math.pi * m.r_di**2)                     # drillstring inner area , m^2
    m.beta_ds  = m.Param(90000)                             # Bulk Modulus of Mud in Drill String [bar]
    m.F_ds     = m.Param(80)                                # Friction Factor in the drill string [bar*s^2/m^6]
    m.rho_ds   = m.Param(rho_pit_[-1])                         # Mud Density in the Drill String [kg/m^3]
    m.A_ann = m.Const(math.pi * (m.r_ci**2 - m.r_do**2))         # annulus flow area , m^2
    m.beta_ann  = m.Param(50000)
    m.F_ann     = m.Param(330)                               # Friction Factor in the Annulus [bar*s^2/m^6]
    m.rho_ann   = m.Param(rho_pit_[-1])                   # Mud Density in the Drill String Annulus [kg/m^3]

    # Reservoir parameters
    m.md  = m.Param(md_[-1])                              # Measured depth (m)
    m.rop = m.Param(rop_[-1])                             # rate of penetration (m/sec)
    m.tvd = m.Param(tvd_[-1])                             # Total vertical depth of well [m]
    m.P_res  = m.Param(P_res_[-1])                              # Formation Pressure (bar)
    m.k   = m.Param(k_[-1])                               # Permeability
    m.el  = m.Param(el_[-1])                              # Effective Length (m)
    m.PIh = m.Param(PIh_[-1])

    # Manipulated Variables
    m.cv = m.MV(value=cv_[-1])                               # Choke Valve Opening from 0-100 [%]
    m.frac_mud = m.MV(value=frac_mud_[-1])
    m.rho_mud = m.MV(value=rho_mud_[-1])
    m.Q_mp = m.MV(Q_mp_[-1])                             # Mud pump flow rate[m^3/min]
    
    # Initialized Variables
    m.P_mp = m.SV(value=P_mp_[-1])
    m.h_pit = m.SV(value=h_pit_[-1])                               # CV - Mud pit level
    m.rho_pit = m.Var(value=rho_pit_[-1])                               # CV - Mud density
    m.Q_bp = m.Param(Q_bp_[-1])                             # Back-pressure pump flow rate [m^3/min]
    
    # Fixed Variables
    m.rho_water = m.FV(value=rho_water_[-1])
    m.rho_res = m.FV(value=rho_res_[-1])
    
    # Calculated Variables
    m.Q_mud = m.Var(m.Q_mp * m.frac_mud)                        # MV - Mud concentrate makeup flow
    m.Q_water = m.Var(m.Q_mp * (1 - m.frac_mud))                    # MV - Water makeup flow
    
    
    # Variables
    m.P_choke = m.Var(2,lb=m.P_atm)                               # Pressure at Choke Valve [bar]
    m.Q_bit = m.Var(m.Q_mp - 60 * m.rop * m.A_ds, lb=0)             # Flow Rate through Bit [m^3/min]
    m.P_bit_init = m.P_choke + (m.rho_ann * (m.F_ann / 3600) * m.md * (m.Q_bit**2) + m.rho_ann * m.g * m.tvd) * 1e-5
    m.P_bit = m.Var(m.P_bit_init)                             # Bit pressure [bar]

    m.Q_res_init = m.PIh * (m.P_res - m.P_bit) / math.log(10.0 / m.r_ci)
    m.Q_res = m.Var(m.Q_res_init) # Reservoir fluid influx flow rate [m^3/min]

    m.Q_choke_init = m.Kc * m.cv * m.sqrt(m.rho_ann * (m.P_choke - m.P_atm) * 1e-5)
    m.Q_choke = m.Var(m.Q_choke_init, lb=0)                   # Flow Rate through Choke [m^3/min]

    # Intermediates
    m.M = m.Intermediate(m.M_ds + m.M_ann)                         # Total Mud Density per length [kg/m^4]
    m.V_ann = m.Intermediate(m.A_ann * m.md)                        # Volume of Annulus [m^3]
    m.V_ds = m.Intermediate(m.A_ds * m.md)                        # Volume of Drill String [m^3]

    #Equations
    # Bit pressure [bar]
    m.Equation(m.P_bit == m.P_choke + (m.rho_ann * (m.F_ann / 3600) * m.md * (m.Q_bit**2) + m.rho_ann * m.g * m.tvd) * 1e-5)

    # Flow Rate through Choke Valve [m^3/min] based on valve characteristics
    m.Equation(m.Q_choke == m.Kc * m.cv * m.sqrt(m.rho_ann * (m.P_choke - m.P_atm) * 1e-5))
    m.Equation(m.Q_res == m.PIh * (m.P_res - m.P_bit))

    # Pressure/flow equations up the annulus
    # Pa_1 == (betaA/Ah) * (Qbit + Qres_bit - Q_1)
    # Pa_i == (betaA/Ah) * (Q_(i-1)  + Qres_i - Q_i)
    # Q_i == Q_(i-1) + Qres_i
    # Pa_N == (betaA/Ah) * (Q_(N-1) - Q_out)
    # Qres_i == K_i * Ah * (PF_i - Pa_i)/EL

    # Mud pump discharge (Equation 5.1)
    m.Equation(m.P_mp.dt() == (m.beta_ds / m.V_ds) * (m.Q_mp - m.Q_bit - 60 * m.rop * m.A_ds) )

    # Choke valve pressure (Equation 5.2)
    m.Equation(m.P_choke.dt() == (m.beta_ann / m.V_ann) * (m.Q_res + m.Q_bit + m.Q_bp - m.Q_choke - 60 * m.rop * m.A_ann))

    # Flow through drill bit (Equation 5.3)
    m.Equation(m.Q_bit.dt() == (1e+5 / m.M) * (m.P_mp - m.P_bit - m.F_ds / 3600 * (m.Q_bit**2) + m.rho_ds * m.g * m.tvd / 1e+5))
    
    #cv = m.MV(value=cv_m)                                   # MV - Choke valve
    #Qbp = m.MV(value=Qbp_m)                                 # MV - Back pressure pump flow
    
    ######## MUDPIT ########
    
    # Model Constants and Parameters
    m.rho_M = m.Const(value=1800.0)                            # makeup mud density (kg/m3)
    m.rho_W = m.Const(value=1000.0)                            # water density (kg/m3)
    #Vpit = m.Const(value=10.0)                             # Mud pit volume (m3)
    m.A_pit = m.Const(value=5.0)                               # Mud pit area (m2)

    # Model Parameters

    # Model Variables
    m.M_pit = m.Var(value=m.h_pit * m.A_pit * m.rho_pit)           # Mud mass in pit (kg)
    m.V_pit = m.Var(value=m.h_pit * m.A_pit)                     # Volume of pit (m3)

    # Model Equations
    m.Equation(m.h_pit.dt() == (1 / m.A_pit) * (m.Q_choke + m.Q_mud + m.Q_water - m.Q_mp) / 60)
    m.Equation(m.M_pit.dt() == (m.rho_pit * m.Q_choke + m.rho_mud * m.Q_mud + m.rho_water * m.Q_water - m.rho_pit * m.Q_mp) / 60)
    m.Equation(m.V_pit.dt() == (m.Q_choke + m.Q_mud + m.Q_water - m.Q_mp) / 60)
    m.Equation(m.rho_pit == m.M_pit / m.V_pit)

###################################################
# Configure MHE
# 120 second time horizon, steps of 3 sec
ntm = 40
mhe.time = np.linspace(0,ntm*3,ntm+1)

# Measured inputs
mhe.Q1.STATUS = 0  # not changed by mhe
mhe.Q1.FSTATUS = 1 # measured

mhe.Q2.STATUS = 0  # not changed by mhe
mhe.Q2.FSTATUS = 1 # measured

mhe.TC1.FSTATUS = 1    # receive measurement
mhe.TC2.FSTATUS = 1    # receive measurement
meas_gap = 2.0
mhe.TC1.MEAS_GAP = meas_gap # for CV_TYPE = 1
mhe.TC2.MEAS_GAP = meas_gap # for CV_TYPE = 1

mhe.options.IMODE     = 5 # MHE Mode
mhe.options.EV_TYPE   = 1 # Estimator Objective type
mhe.options.NODES     = 3 # Collocation nodes
mhe.options.ICD_CALC  = 1 # Calculate initial conditions
mhe.options.SOLVER    = 3 # IPOPT
mhe.options.COLDSTART = 0 # COLDSTART on first cycle

###################################################
# Configure MPC
# Control horizon, non-uniform time steps
mpc.time = [0,3,6,10,14,18,22,27,32,38,45,55,65,75,90,110,130,150]

# update parameters from mhe
mpc.K1.FSTATUS = 1
mpc.K2.FSTATUS = 1
mpc.K3.FSTATUS = 1
mpc.tau12.FSTATUS = 1
mpc.tau3.FSTATUS = 1

# Measured inputs
mpc.Q1.STATUS = 1  # manipulated
mpc.Q1.FSTATUS = 0 # not measured
mpc.Q1.DMAX = 20.0
mpc.Q1.DCOST = 0.1

mpc.Q2.STATUS = 1  # manipulated
mpc.Q2.FSTATUS = 0 # not measured
mpc.Q2.DMAX = 30.0
mpc.Q2.DCOST = 0.1

mpc.TC1.STATUS = 1     # drive to setpoint
mpc.TC1.FSTATUS = 1    # receive measurement
mpc.TC1.TAU = 40       # response speed (time constant)
mpc.TC1.TR_INIT = 1    # reference trajectory
mpc.TC1.TR_OPEN = 0

mpc.TC2.STATUS = 1     # drive to setpoint
mpc.TC2.FSTATUS = 1    # receive measurement
mpc.TC2.TAU = 0        # response speed (time constant)
mpc.TC2.TR_INIT = 0    # dead-band
mpc.TC2.TR_OPEN = 1

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
print('Mud density             =', rho_P_[-1], 'kg/m3')
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
    tvd, rop, Pf, k, el = get_reservoir(res, md_[-1])
    tvd_ = np.append(tvd_, tvd)
    rop_ = np.append(rop_, rop)
    P_res_ = np.append(P_res_, P_res)
    k_ = np.append(k_, k)
    PIh_ = np.append(PIh_, PIh)
    el_ = np.append(el_, el)    
    
    # Insert measurements to MHE
    tvd.MEAS = tvd
    rop.MEAS = rop
    P_res.MEAS = P_res
    k.MEAS = k
    PIh.MEAS = PIh
    el.MEAS = el

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
    

    t_delta_ = np.append(t_delta_, t_delta)
    st_ = np.append(st_, st)
    t_elapsed_ = np.append(t_elapsed_, t_elapsed_[-1] + t_delta_[-1])
    md_ = np.append(md_, md)

    
    Q_mp_ = np.append(Q_mp_, Q_mp_[-1])
    Q_bp_ = np.append(Q_bp_, Q_bp)
    rho_P_ = np.append(rho_P_, rho_P)
    cv_ = np.append(cv_, cv)
    P_mp_ = np.append(P_mp_, P_mp)
    P_choke_ = np.append(P_choke_, P_choke)
    Q_bit_ = np.append(Q_bit_, Q_bit)
    Q_res_ = np.append(Q_res_, Q_res)

    h_pit_ = np.append(h_pit_, h_pit)
    rho_C_ = np.append(rho_C_, rho_P)
    Q_choke_ = np.append(Q_choke_, Q_choke)
    Q_mud_ = np.append(Q_mud_, Q_mud)
    Q_water_ = np.append(Q_water_, Q_water)

    A_pit_ = np.append(A_pit_, A_pit)

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
