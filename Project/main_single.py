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
        rop = rop * 3 

    return tvd, rop, Pf, k, el

def drillstring(md, tvd, rop, Pf, k, el, PIh, Qmp, Qbp, rhoP, cv, tdelta, st, rmt):
    '''
    Arguments:
        pumpQ   - Mud pump flowrate for this interval (m3/min)
        backQ   - Back-pressure pump flowrate for this interval (m3/min)
        chokeVP - Choke Valve Opening from 0-100 [%]
        rhoM    - Mud density [kg/m3]
        depth   - Measured depth at the start of this interval (m)
        dTime   - Length of the interval (seconds)
        st       - simulation time interval (seconds)
    Returns:
        Pp     - Mud pump discharge pressure (bar)
        Pc     - Choke valve inlet pressure (bar)
        Qbit   - Flow through the drill bit (m3/min)
        Pbit   - Pressure at the drill bit (bar)
        Qchoke - Flow through the choke (m3/min)
        MD     - Measured depth at the end of the interval (m)
    '''

    d = GEKKO(remote=rmt)

    # Set up the timeline for the requested interval
    nt = int(tdelta / st) + 1                      # simulation time points
    d.time = np.linspace(0, tdelta, nt)

    # Model Constants
    Patm = d.Const(1)                         # Atmospheric pressure outside choke [bar]
    Md = d.Const(2500)                        # Lumped Density per length of Mud in Drill String [kg/m^4 * 1e5]
    Ma = d.Const(800)                         # Lumped Density per length of Mud in Annulus [kg/m^4 * 1e5]
    g = d.Const(9.81)                         # Gravitational constant (kg-m/s^2)
    r_di = 2.0 * 0.0254                       # drillstring inner radius (m) (4" diameter)
    r_do = 2.5 * 0.0254                       # drillstring outer radius (m) (5" diameter)
    r_ci = 4.3125 * 0.0254                    # annulus inner radius (m) (8 5/8" diameter)
    Ad = d.Const(math.pi * r_di**2)             # drillstring inner area , m^2
    Aa = d.Const(math.pi * (r_ci**2 - r_do**2)) # annulus flow area , m^2
    #Ah = d.Const(math.pi*r_ci**2)             # borehole cross area, m^2

    # Calling Arguments
    Qmp = d.Param(Qmp)                    # Mud pump flow rate[m^3/min]
    Qbp = d.Param(Qbp)                    # Back-pressure pump flow rate [m^3/min]
    cv = d.Param(cv)                     # Choke Valve Opening from 0-100 [%]

    # Parameters from Reservoir Model
    MD_init = np.zeros(len(d.time))
    ROP_init = np.zeros(len(d.time))
    TVD_init = np.zeros(len(d.time))
    PF_init = np.zeros(len(d.time))
    K_init = np.zeros(len(d.time))
    EL_init = np.zeros(len(d.time))
    PIh_init = np.zeros(len(d.time))

    MD_init[0] = md
    TVD_init[0] = tvd
    ROP_init[0] = rop
    PF_init[0] = Pf
    K_init[0] = k
    EL_init[0] = el
    PIh_init[0] = PIh

    for s in range(1, len(d.time)):
        MD_init[s] = MD_init[s-1] + ROP_init[s-1] * (d.time[s] - d.time[s-1])
        TVD_init[s], ROP_init[s], PF_init[s], K_init[s], EL_init[s] = get_reservoir(res, MD_init[s])
        PIh_init[s] = PIh

    MD  = d.Param(MD_init)                    # Measured depth (m)
    ROP = d.Param(ROP_init)                   # rate of penetration (m/sec)
    TVD = d.Param(TVD_init)                   # Total vertical depth of well [m]
    PF  = d.Param(PF_init)                    # Formation Pressure (bar)

    #K   = d.Param(K_init)                     # Permeability
    #EL  = d.Param(EL_init)                    # Effective Length (m)
    PIH = d.Param(PIh_init)

    # Other Model Parameters
    Kc     = d.Param(0.4)                     # Valve Coefficient
    betaD  = d.Param(90000)                   # Bulk Modulus of Mud in Drill String [bar]
    betaA  = d.Param(50000)
    Fd     = d.Param(80)                      # Friction Factor in the drill string [bar*s^2/m^6]
    Fa     = d.Param(330)                     # Friction Factor in the Annulus [bar*s^2/m^6]
    #rhoD   = d.Param(1240)                    # Mud Density in the Drill String [kg/m^3]
    rhoD   = d.Param(rhoP)                    # Mud Density in the Drill String [kg/m^3]
    #rhoA   = d.FV(1290,lb=rhoD)               # Mud Density in the Drill String Annulus [kg/m^3]
    rhoA   = d.FV(rhoP, lb=rhoD)               # Mud Density in the Drill String Annulus [kg/m^3]

    # Variables
    Pp = d.Var(38)                            # Pressure at Pump [bar]
    Pc = d.Var(2,lb=Patm)                     # Pressure at Choke Valve [bar]

    Qbit = d.Var(Qmp - 60 * ROP * Ad, lb=0)      # Flow Rate through Bit [m^3/min]

    Pbit_init = Pc + (rhoA * (Fa / 3600) * md * (Qbit**2) + rhoA * g * TVD_init[0]) * 1e-5
    Pbit = d.Var(Pbit_init)                   # Bit pressure [bar]

    # Reservoir gas influx flow rate [m^3/min]
    #Qres_init = rm.reservoir_flow(200, 1, 4000)
    #Qres_init = K * Ah * (PF - Pbit)/EL
    Qres_init = PIH * (PF - Pbit) / math.log(10.0 / r_ci)
    #Qres_init = rm.reservoir_flow(Pbit.value, math.pi*r_ci**2, depth)
    Qres = d.Var(Qres_init)

    # Flow Rate through Choke Valve [m^3/min]
    Qchoke_init = Kc * cv * d.sqrt(rhoA * (Pc - Patm) * 1e-5)
    Qchoke = d.Var(Qchoke_init, lb=0)          # Flow Rate through Choke [m^3/min]

    # Intermediates
    M = d.Intermediate(Md + Ma)                 # Total Mud Density per length [kg/m^4]
    Va = d.Intermediate(Aa * MD)                # Volume of Annulus [m^3]
    Vd = d.Intermediate(Ad * MD)                # Volume of Drill String [m^3]

    #Equations

    # Bit pressure [bar]
    d.Equation(Pbit == Pc + (rhoA * (Fa / 3600) * MD * (Qbit**2) + rhoA * g * TVD) * 1e-5)

    # Flow Rate through Choke Valve [m^3/min] based on valve characteristics
    d.Equation(Qchoke == Kc * cv * d.sqrt(rhoA * (Pc - Patm) * 1e-5))

    # Flow to/from reservoir based on bit pressure and formation data

    #       K * A * (Pf-Pbit)
    #  Q = -------------------
    #        mu * L

    #        K * (Pf-Pbit) * H      K * H
    #  Q = -------------------- = ---------------- * (Pf-Pbit)
    #       ln(re/rci) * m        ln(re/rci) * mu

    d.Equation(Qres == PIh * (PF - Pbit))

    #r_e = 10  # Effective drainage radius (m)
    #mu = 5.25e-5 * 60 # dynamic viscosity (from reservoir model) [m2/min]
    #mu = 0.042 # kinematic viscosity (kg/ms)
    #d.Equation(Qres == 60 * K * (PF - Pbit)*1e5 / (math.log(r_e/r_ci) * mu) )
    #d.Equation(Qres == 60 * K * Ah * (PF - Pbit)*1e5 / EL)
    #Qres_fixed = 1.0
    #d.Equation(Qres == Qres_fixed)

    # Pressure/flow equations up the annulus
    # Pa_1 == (betaA/Ah) * (Qbit + Qres_bit - Q_1)
    # Pa_i == (betaA/Ah) * (Q_(i-1)  + Qres_i - Q_i)
    # Q_i == Q_(i-1) + Qres_i
    # Pa_N == (betaA/Ah) * (Q_(N-1) - Q_out)
    # Qres_i == K_i * Ah * (PF_i - Pa_i)/EL

    # Change in total vertical depth from formation information
    #d.Equation(TVD.dt() == rm.reservoir_dTVD(MD))

    # Mud pump discharge (Equation 5.1)
    d.Equation(Pp.dt() == (betaD / Vd) * (Qmp - Qbit - 60 * ROP * Ad) )

    # Choke valve pressure (Equation 5.2)
    d.Equation(Pc.dt() == (betaA / Va) * (Qres + Qbit + Qbp - Qchoke - 60 * ROP * Aa))

    # Flow through drill bit (Equation 5.3)
    d.Equation(Qbit.dt() == (1e+5 / M) * (Pp - Pbit - Fd / 3600 * (Qbit**2) + rhoD * g * TVD / 1e+5))

    # Drilling rate from reservior simulation
    #d.Equation(MD.dt() == ROP )

    # Options
    d.options.solver = 3
    d.options.imode = 5                       # dynamic simulation
    d.solve(disp=False)

    # Print solution
    #print("PI =",PI.value[-1])

    return Pp.VALUE[-1], Pc.VALUE[-1], Qbit.VALUE[-1], Pbit.VALUE[-1], Qchoke.VALUE[-1], Qres.VALUE[-1], MD.VALUE[-1]

def mudpit(level_st, rho_st, rho_in, inflow, outflow, mud, water, dTime, st, rmt):

    '''
    Arguments:
        level_st - Mud pit level at start (m)
        rho_st   - Mud density in pit at start (kg/m3)
        rho_in   - Density of inflow from choke (kg/m3)
        inflow   - Flow in from Choke (m3/min)
        outflow  - Flow out to mud pumps (m3/min)
        mudflow  - Mud concentrate makeup (m3/min)
        water    - Fresh water makeup (m3/min)
        dTime    - Length of the interval (seconds)
        st       - simulation time interval (seconds)
    Returns:
        level    - Mud pit level at end (m)
        rho      - Mud density in pit at end (kg/m3)
    '''

    # Specify model
    m = GEKKO(remote=rmt)
    nt = int(dTime / st) + 1 # simulation time points
    m.time = np.linspace(0, dTime, nt)

    # Model Constants and Parameters
    rhoM = m.Const(value=900.0)              # makeup mud density (kg/m3)
    rhoW = m.Const(value=1000.0)              # water density (kg/m3)
    #pitV = m.Const(value=10.0)                # Mud pit volume (m3)
    pitA = m.Const(value=5.0)                 # Mud pit area (m2)

    # Model Parameters

    # Model Variables
    level = m.Var(value=level_st)             # Pit level (m)
    rhoP = m.Var(value=rho_st)                 # Mud density (kg/m3)
    massP = m.Var(value=level_st * pitA * rho_st) # Mud mass in pit (kg)
    volP = m.Var(value=level_st * pitA)       # Volume of pit (m3)

    # Model Equations
    m.Equation(level.dt() == (1 / pitA) * (inflow + mud + water - outflow) / 60)
    m.Equation(massP.dt() == (rho_in * inflow + rhoM * mud + rhoW * water - rhoP * outflow) / 60)
    m.Equation(volP.dt() == (inflow + mud + water - outflow) / 60)
    m.Equation(rhoP == massP / volP)

    # Solver options
    m.options.imode = 5                       # Dynamic Simulation

    # Solve the model
    m.solve(disp=False)

    # Return final state
    return level.VALUE[-1], rhoP.VALUE[-1]

def mhe(dTime, st, rmt):
    '''
    Arguments:
        dTime       - Length of the interval (seconds)
    Returns:
        pitA     - Mud pit area estimate (m2)
                 - Choke valve constant
    '''

    # Specify model
    #m = GEKKO(remote=rmt)
    #nt = int(dTime / st) + 1 # simulation time points
    #m.time = np.linspace(0, dTime, nt)

    # Model Constants and Parameters
    #c1 = 1.0       # constant #1

    # Model Parameters
    #p1 = m.Param(value=0.0)     # parameter #1

    # Model Variables
    #y = m.Var(value=-1.0)       # general variable

    #u = m.MV(value=1.0)         # MV

    #x = m.CV(value=1.0)         # CV

    # Objective
    #term = m.Param(value=np.array([int(t>=tmax) for t in m.time]))
    #m.Obj(term*y*y)
    #m.Obj(term*x*x)
    #m.Obj(term*u*u)

    # Model Equations
    #m.Equation( y.dt() == -y + u )
    #m.Equation( 5.0*x.dt() == -x + u )

    # Tuning

    # MV tuning parameters
    #u.STATUS = 1        # turn MV ON
    #u.DCOST  = 0.01     # move penalty
    #u.DMAX   = 100.0    # maximum move

    # CV tuning parameters
    #x.STATUS = 1        # turn CV ON
    #x.SP   = 0.0        # setpoint for L2 norm
    #x.SPLO = -1.0       # low setpoint for L1 norm
    #x.SPHI = 1.0       # high setpoint for L1 norm
    #x.TR_INIT = 1       # initial equal to the current value on coldstart
    #x.TAU     = 2.0     # speed of SP response

    # Solver options
    #m.options.IMODE = 6     # Dynamic Optimization (Control)
    #m.options.CV_TYPE = 2   # L1 or L2 Norm

    # Solve the model
    #m.solve(disp=False)

    pitA = 5.0  # Placeholder.

    return pitA

def mpc(cv_m, cv_l, cv_h, 
        Qbp_m, Qbp_l, Qbp_h,
        Qmud_m, Qmud_l, Qmud_h,
        Qwater_m, Qwater_l, Qwater_h,
        hpit_m, hpit_l, hpit_h,
        rhoP_m, rhoP_l, rhoP_h, 
        Pchoke_m, Pchoke_l, Pchoke_h, 
        Pbit_m, Pbit_l, Pbit_h, 
        tdelta, st, rmt):
    '''
    Arguments:

        Manipulated Measurements:
        chokeVP_m     - Choke Valve Opening from 0-100 [%]
        backQ_m       - Back-pressure pump flowrate for this interval (m3/min)
        ---> Later
        mudQ_m        - Mud concentrate makeup (m3/min)
        waterQ_m      - Fresh water makeup (m3/min)

        Feed Forward Measurements:
        pumpQ_m       - Mud pump flowrate for this interval (m3/min)

        Dependents:
        level_m, level_l, level_h   - Mud pit level (m)
        Pc_m, Pc_l, Pc_h            - Choke pressure (bar)
        Pbit_m, etc.                - Downhole pressure (bar)
        ---> Later
        rho_m, rho_l, rho_h         - Mud density in pit (kg/m3)

        Other:
        dTime       - Length of the interval (seconds)

    Returns (MVs):
        chokeVP     - Choke Valve Opening from 0-100 [%]
        backQ       - Back-pressure pump flowrate for this interval (m3/min)
        ---> Later
        mudQ        - Mud concentrate makeup (m3/min)
        waterQ      - Fresh water makeup (m3/min)
    '''

    # # Specify controller 
    # mpc = GEKKO(remote=rmt)
    # nt = int(tdelta / st) + 1 # simulation time points
    # mpc.time = np.linspace(0, tdelta, nt)

    # # Model Constants and Parameters

    # # Model Parameters
    # #p1 = mpc.Param(value=0.0)     # parameter #1

    # # Model Variables
    # #y = mpc.Var(value=-1.0)       # general variable

    # cv = mpc.MV(value=cv_m)   # MV - Choke valve
    # Qbp = mpc.MV(value=Qbp_m)       # MV - Back pressure pump flow
    # Qmud = mpc.MV(value=Qmud_m)         # MV - Mud concentrate makeup flow
    # Qwater = mpc.MV(value=Qwater_m)     # MV - Water makeup flow

    # hpit = mpc.CV(value=hpit_m)       # CV - Mud pit level
    # Pchoke = mpc.CV(value=Pchoke_m)             # CV - Pressure at choke
    # rhoP = mpc.CV(value=rhoP_m)           # CV - Mud density

    # # Objective
    # #term = m.Param(value=np.array([int(t>=tmax) for t in m.time]))
    # #mpc.Obj(term*y*y)
    # #mpc.Obj(term*x*x)
    # #mpc.Obj(term*u*u)

    # # Model Equations
    # #mpc.Equation( y.dt() == -y + u )
    # #mpc.Equation( 5.0*x.dt() == -x + u )

    # # Tuning

    # # MV tuning parameters
    # #u.STATUS = 1        # turn MV ON
    # #u.DCOST  = 0.01     # move penalty
    # #u.DMAX   = 100.0    # maximum move

    # # CV tuning parameters
    # hpit.STATUS = 1        # turn CV ON
    # hpit.SP   = 0.0        # setpoint for L2 norm
    # hpit.SPLO = hpit_l    # low setpoint for L1 norm
    # hpit.SPHI = hpit_h    # high setpoint for L1 norm
    # hpit.TR_INIT = 1       # initial equal to the current value on coldstart
    # hpit.TAU     = 2.0     # speed of SP response

    # # Solver options
    # mpc.options.IMODE = 6     # Dynamic Optimization (Control)
    # mpc.options.CV_TYPE = 2   # L1 or L2 Norm

    # # Solve the model
    # mpc.solve(disp=False)

    # # Make up the outputs so everything works
    # #cv = 20.0
    # #Qbp   = 0.4
    # #Qmud = 0.0
    # #Qwater = 0.0   

    # return cv.VALUE[-1], Qbp.VALUE[-1], Qmud.VALUE[-1], Qwater.VALUE[-1]
    
    if hpit_m < hpit_l:
        cv_new = min(cv_m + 10, 100)
    elif hpit_m > hpit_h:
        cv_new = max(cv_m - 10, 0)
    else:
        cv_new = cv_m
        
    if cv_new < 20:
        Qbp_new = min(Qbp_m + 1, 10)
    elif cv_new < 80:
        Qbp_new = max(Qbp_m - 1, 0)
    else:
        Qbp_new = Qbp_m
    
    return cv_new, Qbp_new, 0.0, 0.0 

# Initialize things
Dmin = 2000.0
Dmax = 4000.0                               # m
tdelta = 60.0 * 60.0 * 0.1                      # seconds
st = 180.0                                    # seconds

tdelta_ = np.array([tdelta])
st_ = np.array([st])
telapsed_ = np.array([0])
md_ = np.array([Dmin])

tvd_ = np.array([0.0])
rop_ = np.array([0.0])
Pf_ = np.array([0.0])
k_ = np.array([0.0])
el_ = np.array([0.0])
PIh_ = np.array([0.0])
mu = 0.042                                    # kinematic viscosity (kg/m*s)

Qmp_ = np.array([2.0])
Qbp_ = np.array([0.4])
rhoP_ = np.array([900.0])
cv_ = np.array([20.0])
Pmp_ = np.array([0.0])
Pchoke_ = np.array([0.0])
Qbit_ = np.array([0.0])
Qres_ = np.array([0.0])
Pdh_ = np.array([0.0])

hpit_ = np.array([10])
rhoC_ = np.array([900.0])
Qchoke_ = np.array([0.0])
Qmud_ = np.array([0.0])
Qwater_ = np.array([0.0])

Apit_ = np.array([0.0])

res = load_reservoir()

drilling = True
rmt = True

print('-[Process Initialization]-----------------------')
print('Time elapsed            =', telapsed_[-1], 's')
print('Measured depth          =', md_[-1], 'm')
print('Mud pit level           =', hpit_[-1], 'm')
print('Mud density             =', rhoP_[-1], 'kg/m3')
print('Pressure at choke valve =', Pchoke_[-1], 'bar')
print('Flow rate through choke =', Qchoke_[-1], 'm3/min')
print('Downhole pressure       =', Pdh_[-1], 'bar')

# main loop
while drilling:

    # Get the relevant reservoir properties
    tvd, rop, Pf, k, el = get_reservoir(res, md_[-1])
    
    PIh = calc_weighted_PIh(md_, k_, mu)

    # Call the drill string model
    Pmp, Pchoke, Qbit, Pdh, Qchoke, Qres, md = drillstring(md_[-1],
                              tvd, rop, Pf, k, el, PIh, 
                              Qmp_[-1], Qbp_[-1],
                              rhoP_[-1], cv_[-1],
                              tdelta_[-1], st_[-1], rmt)

    # Call the mudpit model
    hpit, rhoP = mudpit(hpit_[-1], rhoP_[-1], rhoC_[-1],
                             Qchoke, Qmp_[-1] + Qbp_[-1],
                             Qmud_[-1], Qwater_[-1], 
                             tdelta_[-1], st_[-1], rmt)
    
    # Call the MHE
    Apit = mhe(tdelta_[-1], st_[-1], rmt)

    # Call the MPC
    cv, Qbp, Qmud, Qwater = mpc(cv_[-1], 1.0, 99.0,
                            Qbp_[-1], 0.0, 4.0,
                            Qmud_[-1], 0.0, 4.0,
                            Qwater_[-1], 0.0, 4.0,
                            hpit_[-1], 9.5, 10.5,
                            rhoP_[-1], 800.0, 1500.0,
                            Pchoke_[-1], 0.0, 5.0,
                            Pdh_[-1], 100.0, 1000.0,
                            tdelta_[-1], st_[-1], rmt)

    r_ci = 4.3125 * 0.0254            # annulus inner radius (m) (8 5/8" diameter)
    Aa = math.pi * (r_ci**2) # annulus flow area , m^2
    Qmud = rop * Aa * 60

    tdelta_ = np.append(tdelta_, tdelta)
    st_ = np.append(st_, st)
    telapsed_ = np.append(telapsed_, telapsed_[-1] + tdelta_[-1])
    md_ = np.append(md_, md)

    tvd_ = np.append(tvd_, tvd)
    rop_ = np.append(rop_, rop)
    Pf_ = np.append(Pf_, Pf)
    k_ = np.append(k_, k)
    PIh_ = np.append(PIh_, PIh)
    el_ = np.append(el_, el)

    Qmp_ = np.append(Qmp_, Qmp_[-1])
    Qbp_ = np.append(Qbp_, Qbp)
    rhoP_ = np.append(rhoP_, rhoP)
    cv_ = np.append(cv_, cv)
    Pmp_ = np.append(Pmp_, Pmp)
    Pchoke_ = np.append(Pchoke_, Pchoke)
    Qbit_ = np.append(Qbit_, Qbit)
    Qres_ = np.append(Qres_, Qres)
    Pdh_ = np.append(Pdh_, Pdh)

    hpit_ = np.append(hpit_, hpit)
    rhoC_ = np.append(rhoC_, rhoP)
    Qchoke_ = np.append(Qchoke_, Qchoke)
    Qmud_ = np.append(Qmud_, Qmud)
    Qwater_ = np.append(Qwater_, Qwater)

    Apit_ = np.append(Apit_, Apit)

    print('-[Process]--------------------------------------')
    print('Time elapsed            =', telapsed_[-1] / 3600, 'hr')
    print('Measured depth          =', md_[-1], 'm')
    print('Mud pit level           =', hpit_[-1], 'm')
    print('Mud density             =', rhoP_[-1], 'kg/m3')
    print('Pressure at choke valve =', Pchoke_[-1], 'bar')
    print('Flow rate through choke =', Qchoke_[-1], 'm3/min')
    print('Downhole pressure       =', Pdh_[-1], 'bar')

    print('-[Reservoir]------------------------------------')
    print('TVD                     =', tvd_[-1], 'm')
    print('Rate of penetration     =', rop_[-1] * 3600, 'm/hr')
    print('Pressure of formation   =', Pf_[-1], 'bar')
    print('Permeability            =', k_[-1], 'm/s')
    print('Effective length        =', el_[-1], 'm')
    print('Reservoir flow          =', Qres_[-1], 'm3/min')

    print('-[MHE]------------------------------------------')
    print('Pit area                =', Apit_[-1], 'm2')

    print('-[MPC]------------------------------------------')
    print('Choke valve position    =', cv_[-1], '%')
    print('Back-pressure pump flow =', Qbp_[-1], 'm3/min')
    print('Qmud                    =', Qmud_[-1], 'm3/min')


    # Go back to the top
    print('------------------------------------------------\n')

    if md > Dmax:
        drilling = False

plt.figure(1)

plt.subplot(6,1,1)
plt.plot(telapsed_ / 3600, md_, 'r-', label='Measured depth')
plt.plot(telapsed_ / 3600, tvd_, 'b-', label='True vertical depth')
plt.ylabel('Depth (m)')
plt.legend(loc='upper left')

plt.subplot(6,1,2)
plt.plot(telapsed_ / 3600, Pmp_, 'k-', label='Mud pump discharge pressure')
plt.plot(telapsed_ / 3600, Pdh_, 'r-', label='Drill bit pressure')
plt.plot(telapsed_ / 3600, Pf_, 'g:', label='Formation pressure')
plt.plot(telapsed_ / 3600, Pchoke_, 'b-', label='Choke valve inlet pressure')
plt.ylabel('Pressure (bar)')
plt.legend(loc='upper left')

plt.subplot(6,1,3)
plt.plot(telapsed_ / 3600, Qmp_, 'k-', label='Mud pump flowrate')
plt.plot(telapsed_ / 3600, Qbit_, 'r-', label='Drill bit flowrate')
plt.plot(telapsed_ / 3600, Qres_, 'g:', label='Formation flowrate')
plt.plot(telapsed_ / 3600, Qchoke_, 'b-', label='Choke flowrate')
plt.plot(telapsed_ / 3600, Qbp_, 'y-', label='Back-pressure pump flowrate')
plt.ylabel('Flow (m3/min)')
plt.legend(loc='upper left')

plt.subplot(6,1,4)
plt.plot(telapsed_ / 3600, rop_, 'k-', label='ROP')
plt.ylabel('ROP (m/hr)')
plt.legend(loc='upper left')

plt.subplot(6,1,5)
plt.plot(telapsed_ / 3600, hpit_, 'k-', label='Pit Level')
plt.ylabel('Height (m)')
plt.legend(loc='upper left')

plt.subplot(6,1,6)
plt.plot(telapsed_ / 3600, cv_, 'b-', label='Choke Position')
plt.ylabel('%')
plt.legend(loc='upper left')

plt.xlabel('Time (hr)')

plt.show()
