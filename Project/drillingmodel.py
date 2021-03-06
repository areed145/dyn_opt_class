#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# reservoirmodel.py - A simple table-lookup model of an oil reservoir
#
__version__   = "1.0.0"
__date__      = "2019.03.18"

import math
import numpy as np
import reservoirmodel as rm
from gekko import GEKKO
import matplotlib.pyplot as plt


#%% Drillstring Model
def drillstring(pumpQ, backQ, rhoM, chokeVP, depth, dTime):
    """
    Arguments:
        pumpQ   - Mud pump flowrate for this interval (m3/min)
        backQ   - Back-pressure pump flowrate for this interval (m3/min)
        chokeVP - Choke Valve Opening from 0-100 [%]
        rhoM    - Mud density [kg/m3]
        depth   - Measured depth at the start of this interval (m)
        dTime   - Length of the interval (seconds)
    Returns:
        Pp     - Mud pump discharge pressure (bar)
        Pc     - Choke valve inlet pressure (bar)
        Qbit   - Flow through the drill bit (m3/min)
        Pbit   - Pressure at the drill bit (bar)
        Qchoke - Flow through the choke (m3/min)
        MD     - Measured depth at the end of the interval (m)
    """
    rmt = True  # Solve local or remote

    d = GEKKO(remote=rmt)

    # Set up the timeline for the requested interval
    st = 10.0   # simulation time interval (seconds)
    nt = int(dTime/st)+1 # simulation time points
    d.time = np.linspace(0,dTime,nt)
    
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
    #Ah = d.Const(math.pi*r_ci**2)  # borehole cross area, m^2
    
    # Calling Arguments
    Qpump = d.Param(pumpQ)          # Mud pump flow rate[m^3/min] 
    Qback = d.Param(backQ)          # Back-pressure pump flow rate [m^3/min]
    Zc = d.Param(chokeVP)           # Choke Valve Opening from 0-100 [%]

    # Parameters from Reservoir Model
    MD_init = np.zeros(len(d.time))
    ROP_init = np.zeros(len(d.time))
    TVD_init = np.zeros(len(d.time))
    PF_init = np.zeros(len(d.time))
    K_init = np.zeros(len(d.time))
    EL_init = np.zeros(len(d.time))
    #PI_init = np.zeros(len(d.time))

    MD_init[0] = depth
    TVD_init[0], ROP_init[0], PF_init[0], K_init[0], EL_init[0] \
                            = rm.reservoir(depth)
    
    for s in range(1,len(d.time)):
        MD_init[s] = MD_init[s-1] + ROP_init[s-1] * ( d.time[s] - d.time[s-1] )
        TVD_init[s], ROP_init[s], PF_init[s], K_init[s], EL_init[s] \
                            = rm.reservoir(MD_init[s])

    mu = 0.042  # kinematic viscosiry (kg/m*s)
    PI_init = K_init * 1e3 * 1e5 / (math.log(10.0/r_ci)*mu)

        
    MD  = d.Param(MD_init)      # Measured depth (m)
    ROP = d.Param(ROP_init)     # rate of penetration (m/sec)
    TVD = d.Param(TVD_init)     # Total vertical depth of well [m]
    PF  = d.Param(PF_init)      # Formation Pressure (bar)

    #K   = d.Param(K_init)       # Permeability
    #EL  = d.Param(EL_init)      # Effective Length (m)
    PI = d.Param(PI_init)
    
    # Other Model Parameters
    Kc     = d.Param(0.4)       # Valve Coefficient 
    betaD  = d.Param(90000)      # Bulk Modulus of Mud in Drill String [bar]
    betaA  = d.Param(50000)    
    Fd     = d.Param(80)         # Friction Factor in the drill string [bar*s^2/m^6]
    Fa     = d.Param(330)        # Friction Factor in the Annulus [bar*s^2/m^6]
    #rhoD   = d.Param(1240)       # Mud Density in the Drill String [kg/m^3]        
    rhoD   = d.Param(rhoM)       # Mud Density in the Drill String [kg/m^3]        
    #rhoA   = d.FV(1290,lb=rhoD)  # Mud Density in the Drill String Annulus [kg/m^3]    
    rhoA   = d.FV(rhoM,lb=rhoD)  # Mud Density in the Drill String Annulus [kg/m^3]    
    

    # Variables
    Pp = d.Var(38)                # Pressure at Pump [bar]
    Pc = d.Var(2,lb=Patm)         # Pressure at Choke Valve [bar]

    Qbit = d.Var(Qpump - 60*ROP*Ad,lb=0)    # Flow Rate through Bit [m^3/min]
    
    Pbit_init = Pc + (rhoA*(Fa/3600)*depth*(Qbit**2) + rhoA*g*TVD_init[0])*1e-5
    Pbit = d.Var(Pbit_init) # Bit pressure [bar]

    # Reservoir gas influx flow rate [m^3/min]
    #Qres_init = rm.reservoir_flow(200, 1, 4000)
    #Qres_init = K * Ah * (PF - Pbit)/EL
    Qres_init = PI * (PF - Pbit)
    #Qres_init = rm.reservoir_flow(Pbit.value, math.pi*r_ci**2, depth)
    Qres = d.Var(Qres_init)
    
    # Flow Rate through Choke Valve [m^3/min]
    Qchoke_init = Kc * Zc * d.sqrt(rhoA*(Pc-Patm)*1e-5) 
    Qchoke = d.Var(Qchoke_init,lb=0)        # Flow Rate through Choke [m^3/min]
    
    # Intermediates
    M = d.Intermediate(Md+Ma)   # Total Mud Density per length [kg/m^4]
    Va = d.Intermediate(Aa*MD)  # Volume of Annulus [m^3]
    Vd = d.Intermediate(Ad*MD)  # Volume of Drill String [m^3]
    
    #Equations
    
    # Bit pressure [bar]
    d.Equation(Pbit == Pc + (rhoA*(Fa/3600)*MD*(Qbit**2) + rhoA*g*TVD)*1e-5)

    # Flow Rate through Choke Valve [m^3/min] based on valve characteristics
    d.Equation(Qchoke == Kc * Zc * d.sqrt(rhoA*(Pc-Patm)*1e-5))

    # Flow to/from reservoir based on bit pressure and formation data
    
    #       K * A * (Pf-Pbit)
    #  Q = -------------------
    #        mu * L
    
    #        K * (Pf-Pbit) * H      K * H
    #  Q = -------------------- = ---------------- * (Pf-Pbit)
    #       ln(re/rci) * m        ln(re/rci) * mu
    
    d.Equation(Qres == PI * (PF - Pbit))
    
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
    d.Equation( Pp.dt() == (betaD/Vd) * (Qpump - Qbit - 60*ROP*Ad) )
    
    # Choke valve pressure (Equation 5.2)
    d.Equation( Pc.dt() == (betaA/Va) * (Qres + Qbit + Qback - Qchoke - 60*ROP*Aa))
    
    # Flow through drill bit (Equation 5.3)
    d.Equation( Qbit.dt() == (1e+5/M) * (Pp - Pbit - Fd/3600*(Qbit**2) \
                        + rhoD*g*TVD/1e+5 ) )
    
    # Drilling rate from reservior simulation
    #d.Equation( MD.dt() == ROP )
    
    
    # Options
    #d.options.solver = 1
    #d.options.imode = 3    # Calculate starting conditions
    d.options.imode = 4     # dynamic simulation
    d.solve(disp=False)
    
    # Print solution
    #print("PI =",PI.value[-1])

    # plt.figure(1)
    
    # plt.subplot(6,1,1)
    # plt.plot(d.time[2:],Qpump[2:],'b-',label='Mud Pump Flow')
    # plt.plot(d.time[2:],Qbit[2:],'g:',label='Bit Mud Flow')
    # plt.plot(d.time[2:],Qchoke[2:],'r--',label='Choke Mud Flow')
    # plt.ylabel(r'Flow ($m^3/min$)')
    # plt.legend(loc='best')
    
    # plt.subplot(6,1,2)
    # plt.plot(d.time[2:],Zc[2:],'k-',label='Choke Opening (%)')
    # plt.ylabel('Choke (%)')
    # plt.legend(loc='best')
    
    # plt.subplot(6,1,3)
    # plt.plot(d.time[2:],Pbit[2:],'r-',label='Bit Pressure (bar)')
    # plt.ylabel('Press (bar)')
    # plt.legend(loc='best')
    
    # plt.subplot(6,1,4)
    # plt.plot(d.time[2:],Pp[2:],'r:',label='Pump Pressure (bar)')
    # plt.plot(d.time[2:],Pc[2:],'b--',label='Choke Pressure (bar)')
    # plt.ylabel('Press (bar)')
    # plt.legend(loc='best')


    # plt.subplot(6,1,5)
    # plt.plot(d.time[2:],Qres[2:],'k-',label='Reservoir Flow (m3/min)')
    # plt.ylabel('Qres (m3/min)')
    # plt.legend(loc='best')
    
    # plt.subplot(6,1,6)
    # plt.plot(d.time[2:],MD[2:],'k-',label='Measured Depth (m)')
    # plt.ylabel('Measured Depth (m)')
    # plt.legend(loc='best')
    
    # plt.xlabel('Time (sec)')
    # plt.show()

    return (Pp.VALUE[-1], \
            Pc.VALUE[-1], \
            Qbit.VALUE[-1], \
            Pbit.VALUE[-1], \
            Qchoke.VALUE[-1], \
            Qres.VALUE[-1], \
            MD.VALUE[-1] )

#%% Utility functions
def maxdepth():
    return rm.reservoir_maxdepth()

#%% Run   

def test():
    main()
    return    

def main():
    mud_pump_flow = 1.0 #2004.0*(1e-3/60)
    bp_pump_flow = 0.2 #804.0*(1e-3/60)
    mud_density = 900
    choke_valve = 50.0
    #meas_depth = 6000.0*.3048  # Normal
    #meas_depth = 8000.0*.3048  # loss of circ
    meas_depth = 15199.5*.3048  # Kick
    time_interval = 5.0

    print('-[Inputs]---------------------------------------')
    print('Flow Rate through Mud Pump [m^3/min] =', mud_pump_flow)
    print('Flow Rate through Back-pressure Pump [m^3/min] =', bp_pump_flow)
    print('Mud density [kg/m3] =', mud_density)
    print('Choke Valve Position [%] =', choke_valve)
    print('------------------------------------------------')
    print('Time interval [minutes] =', time_interval)
    print('Initial measured depth =', meas_depth,'m',meas_depth/0.3048,'ft')
    print('------------------------------------------------\n')
    
    print('Calling reservoir.reservoir()...')
    TVD, ROP, PF, K, EL = rm.reservoir(meas_depth)
    
    print('-[Results]--------------------------------------')
    print('True Vertical Depth [m] =', TVD)
    print('Rate of Penetration [m/s] =', ROP)
    print('Formation Pressure [bar] =', PF)
    print('Permeability [m2] =', K)
    print('Effective Length [m] =', EL)
    print('------------------------------------------------\n')
    
    print('Calling drillstring()...')
    Pp, Pc, Qb, Pb, Qc, Qr, md = drillstring(mud_pump_flow, bp_pump_flow, \
                                         mud_density, choke_valve, meas_depth, \
                                         time_interval*60 )
    
    print('-[Drilling Results]--------------------------------------')
    print('Pressure at Pump [bar] =', Pp)
    print('Pressure at Choke Valve [bar] =', Pc)
    print('Bit pressure [bar] =', Pb)
    print('Flow Rate through Bit [m^3/min] =', Qb)
    print('Flow Rate through Choke [m^3/min] =', Qc)
    print('Flow Rate from reservoir [m^3/min] =', Qr)
    print('Flow Rate from reservoir [bbl/h] =', Qr*60/0.159)
    print('Measured depth =', md,'m',md/0.3048,'ft')
    print()
    print("delta Pressure [bar ... psi]", Pb-PF, "...", (Pb-PF)*14.7)
    print()
    print("RoC of Mud Pit Volume [m3/min]", \
              (Qc - mud_pump_flow - bp_pump_flow))
    print("RoC of Mud Pit Volume [m3/day]", \
              (Qc - mud_pump_flow - bp_pump_flow) * 1440)
    print("RoC of Mud Pit Volume [bbl/h]", \
              (Qc - mud_pump_flow - bp_pump_flow) * 60 / 0.159)
    print('------------------------------------------------\n')
    
    return     

#%% This is only run whne script is executed as a standalone program
if __name__ == '__main__':
    #import sys
    import os, traceback, argparse
    import time
    #import re
    #from pexpect import run, spawn

    try:
        start_time = time.time()
        parser = argparse.ArgumentParser(description='A simple table-lookup reservoir model')
        parser.add_argument('--version', action='version', version='%(prog)s v'+__version__)
        parser.add_argument ('-v', '--verbose', action='store_true', help='produce verbose output')
        parser.add_argument ('-t', '--test', action='store_true', help='run test suite')
        args = parser.parse_args()
        if args.verbose: print (time.asctime())
        if args.test: 
            test()
        else:
            main()
        if args.verbose: print (time.asctime())
        if args.verbose: print ('Elapsed time:', \
                                (time.time() - start_time), 'seconds')

    except KeyboardInterrupt as e: # Ctrl-C
        raise e

    except SystemExit as e: # sys.exit()
        raise e

    except Exception as e:
        print ('ERROR, UNEXPECTED EXCEPTION')
        print (str(e))
        traceback.print_exc()
        os._exit(1)
