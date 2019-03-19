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


#%% Drillstring Model
def drillstring(pumpQ, backQ, chokeVP, depth, dTime):
    """
    Arguments:
        pumpQ   - Mud pump flowrate for this interval (m3/min)
        backQ   - Back-pressure pump flowrate for this interval (m3/min)
        chokeVP - Choke Valve Opening from 0-100 [%]
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
    rmt = False  # Solve local or remote

    d = GEKKO(remote=rmt)

    # Set up the timeline for the requested interval
    st = 2.0   # simulation time interval (seconds)
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
    Ah = d.Const(math.pi*r_ci**2)  # borehole cross area, m^2
    
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

    MD_init[0] = depth
    TVD_init[0], ROP_init[0], PF_init[0], K_init[0], EL_init[0] = rm.reservoir(depth)
    
    for s in range(1,len(d.time)):
        MD_init[s] = MD_init[s-1] + ROP_init[s-1] * ( d.time[s] - d.time[s-1] )
        TVD_init[s], ROP_init[s], PF_init[s], K_init[s], EL_init[s] = rm.reservoir(MD_init[s])
        
    MD  = d.Param(MD_init)      # Measured depth (m)
    ROP = d.Param(ROP_init)     # rate of penetration (m/min)
    TVD = d.Param(TVD_init)     # Total vertical depth of well [m]
    PF  = d.Param(PF_init)      # Formation Pressure (bar)
    K   = d.Param(K_init)       # Permeability
    EL  = d.Param(EL_init)      # Effective Length (m)
    
    # Other Model Parameters
    Kc     = d.Param(0.4)           # Valve Coefficient 
    betaD  = d.Param(90000)          # Bulk Modulus of Mud in Drill String [bar]
    betaA  = d.Param(50000)    
    Fd     = d.Param(80)             # Friction Factor in the drill string [bar*s^2/m^6]
    Fa     = d.Param(330)            # Friction Factor in the Annulus [bar*s^2/m^6]
    rhoD   = d.Param(1240)          # Mud Density in the Drill String [kg/m^3]        
    rhoA   = d.FV(1290,lb=rhoD)     # Mud Density in the Drill String Annulus [kg/m^3]    
    

    # Variables
    Pp = d.Var(38)                # Pressure at Pump [bar]
    Pc = d.Var(2,lb=Patm)         # Pressure at Choke Valve [bar]

    Qbit = d.Var(Qpump - ROP*Ad,lb=0)    # Flow Rate through Bit [m^3/min]
    
    Pbit_init = Pc + (rhoA*(Fa/3600)*depth*(Qbit**2) + rhoA*g*TVD_init[0])*1e-5
    Pbit = d.Var(Pbit_init) # Bit pressure [bar]

    # Reservoir gas influx flow rate [m^3/min]
    #Qres_init = rm.reservoir_flow(200, 1, 4000)
    Qres_init = K * Ah * (PF - Pbit)/EL
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

    # Flow to/from reservior based on bit pressure and formation data
    d.Equation(Qres == K * Ah * (PF - Pbit)/EL)

    # Change in total vertical depth from formation information
    #d.Equation(TVD.dt() == rm.reservoir_dTVD(MD))

    # Mud pump discharge (Equation 5.1)
    d.Equation( Pp.dt() == (betaD/Vd) * (Qpump - Qbit - ROP*Ad) )
    
    # Choke valve pressure (Equation 5.2)
    d.Equation( Pc.dt() == (betaA/Va) * (Qres + Qbit + Qback - Qchoke - ROP*Aa))
    
    # Flow through drill bit (Equation 5.3)
    d.Equation( Qbit.dt() == (1e+5/M) * (Pp - Pbit - Fd/3600*(Qbit**2) \
                        + rhoD*g*TVD/1e+5 ) )
    
    # Drilling rate from reservior simulation
    #d.Equation( MD.dt() == ROP )
    
    
    # Options
    #d.options.solver = 1
    #d.options.imode = 3    # Calculate starting conditions
    d.options.imode = 4     # dynamic simulation
    d.solve()
    
    # Print solution
   

    return (Pp.VALUE[-1], Pc.VALUE[-1], \
            Qbit.VALUE[-1], Pbit.VALUE[-1], \
            Qchoke.VALUE[-1], MD.VALUE[-1] )
    

#%% Run   

def test():
    main()
    return    

def main():
    mud_pump_flow = 2004.0*(1e-3/60)
    bp_pump_flow = 804.0*(1e-3/60)
    choke_valve = 30.0
    meas_depth = 4000.0
    time_interval = 20.0

    print('-[Inputs]---------------------------------------')
    print('Flow Rate through Mud Pump [m^3/min] =', mud_pump_flow)
    print('Flow Rate through Back-pressure Pump [m^3/min] =', bp_pump_flow)
    print('Choke Valve Position [%] =', choke_valve)
    print('Initial measured depth [m] =', meas_depth)
    print('Time interval [minutes] =', time_interval)
    print('------------------------------------------------\n')
    
    print('Calling drillstring()...')
    Pp, Pc, Qb, Pb, Qc, md = drillstring(mud_pump_flow, bp_pump_flow, \
                                         choke_valve, meas_depth, \
                                         time_interval*60 )
    
    print('-[Results]--------------------------------------')
    print('Pressure at Pump [bar] =', Pp)
    print('Pressure at Choke Valve [bar] =', Pc)
    print('Bit pressure [bar] =', Pb)
    print('Flow Rate through Bit [m^3/min] =', Qb)
    print('Flow Rate through Choke [m^3/min] =', Qc)
    print('Measured Depth [m] =', md)
    print('------------------------------------------------\n')
    
    return     

#%% This is only run whne script is executed as a standalone program
if __name__ == '__main__':
    import sys, os, traceback, argparse
    import time
    #import matplotlib.pyplot as plt
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
        #sys.exit(0)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print ('ERROR, UNEXPECTED EXCEPTION')
        print (str(e))
        traceback.print_exc()
        os._exit(1)
