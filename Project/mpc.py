#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# mpc.py MPC for drilling problem
#
__version__   = "1.0.1"
__date__      = "2019.03.17"

"""
SYNOPSIS

    mpc [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    This is the model predictive controller (MPC) for the 
    managed-pressure drilling problem

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

LICENSE

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# Version   Date        Notes
# -------   ----------  -------------------------------------------------------
# 1.0.0     2013.12.01  Starting script template
# 1.0.1     2019-03-17  Adapted as a GEKKO template
# 1.0.2     2019-04-02  First version of MPC function
#

import numpy as np
from gekko import GEKKO
#import matplotlib.pyplot as plt

def mpc(chokeVP_m, chokeVP_l, chokeVP_h, \
        backQ_m, backQ_l, backQ_h, \
        level_m, level_l, level_h, \
        Pc_m, Pc_l, Pc_h, \
        Pbit_m, Pbit_l, Pbit_h, \
        dTime):
#def mpc(dTime):
#        rho_m, rho_l, rho_h, \
    """
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
        rho_m, rho_l,rho_h          - Mud density in pit (kg/m3)

        Other:
        dTime       - Length of the interval (seconds)
 
   Returns (MVs):
        chokeVP     - Choke Valve Opening from 0-100 [%]
        backQ       - Back-pressure pump flowrate for this interval (m3/min)
        ---> Later
        mudQ        - Mud concentrate makeup (m3/min)
        waterQ      - Fresh water makeup (m3/min)

    """

    #%% Non-model parameters
    rmt = True  # Solve local or remote
    npts = 100    # time steps

    #%% Specify controller 
    mpc = GEKKO(remote=rmt)
    mpc.time = np.linspace(0, dTime, npts+1)  # Model Timeline

    # Model Constants and Parameters

    # Model Parameters
    #p1 = mpc.Param(value=0.0)     # parameter #1

    # Model Variables
    #y = mpc.Var(value=-1.0)       # general variable

    #chokeVP = mpc.MV(value=chokeVP_m)   # MV - Choke valve
    #backQ = mpc.MV(value=backQ_m)       # MV - Back pressure pump flow
    #mudQ = mpc.MV(value=mudQ_m)         # MV - Mud concentrate makeup flow
    #waterQ = mpc.MV(value=waterQ_m)     # MV - Water makeup flow

    #level = mpc.CV(value=level_m)       # CV - Mud pit level
    #Pc = mpc.CV(value=Pc_m)             # CV - Pressure at choke
    #rho = mpc.CV(value=rho_m)           # CV - Mud density

    # Objective
    #term = m.Param(value=np.array([int(t>=tmax) for t in m.time]))
    #mpc.Obj(term*y*y)
    #mpc.Obj(term*x*x)
    #mpc.Obj(term*u*u)

    # Model Equations
    #mpc.Equation( y.dt() == -y + u )
    #mpc.Equation( 5.0*x.dt() == -x + u )

    # Tuning

    # MV tuning parameters
    #u.STATUS = 1        # turn MV ON
    #u.DCOST  = 0.01     # move penalty
    #u.DMAX   = 100.0    # maximum move

    # CV tuning parameters
    #level.STATUS = 1        # turn CV ON
    #level.SP   = 0.0        # setpoint for L2 norm
    #level.SPLO = level_l    # low setpoint for L1 norm
    #level.SPHI = level_h    # high setpoint for L1 norm
    #level.TR_INIT = 1       # initial equal to the current value on coldstart
    #level.TAU     = 2.0     # speed of SP response

    # Solver options
    mpc.options.IMODE = 6     # Dynamic Optimization (Control)
    mpc.options.CV_TYPE = 2   # L1 or L2 Norm

    # Solve the model
    #mpc.solve(disp=False)

    # Make up the outputs so everything works
    chokeVP = 30.0
    backQ   = 0.4
    #mudQ = 0.0
    #waterQ = 0.0

    return (chokeVP, backQ)  #, mudQ, waterQ)

#%% The main function in this module
def main ():

    global options, args

"""
    #%% Non-model parameters
    rmt = False  # Solve local or remote
    tf = 10.0     # Final time
    npts = 100    # time steps
    tmax = 8.0    # end point for "terminal" kinds of limits 

    #%% Specify model
    m = GEKKO(remote=rmt)
    m.time = np.linspace(0, tf, npts+1)  # Model Timeline

    # Model Constants and Parameters
    c1 = 1.0       # constant #1

    # Model Parameters
    p1 = m.Param(value=0.0)     # parameter #1

    # Model Variables
    y = m.Var(value=-1.0)       # general variable

    u = m.MV(value=1.0)         # MV

    x = m.CV(value=1.0)         # CV

    # Objective
    term = m.Param(value=np.array([int(t>=tmax) for t in m.time]))
    m.Obj(term*y*y)
    m.Obj(term*x*x)
    m.Obj(term*u*u)

    # Model Equations
    m.Equation( y.dt() == -y + u )
    m.Equation( 5.0*x.dt() == -x + u )

    # Tuning

    # MV tuning parameters
    u.STATUS = 1        # turn MV ON
    u.DCOST  = 0.01     # move penalty
    u.DMAX   = 100.0    # maximum move

    # CV tuning parameters
    x.STATUS = 1        # turn CV ON
    x.SP   = 0.0        # setpoint for L2 norm
    x.SPLO = -1.0       # low setpoint for L1 norm
    x.SPHI = 1.0       # high setpoint for L1 norm
    x.TR_INIT = 1       # initial equal to the current value on coldstart
    x.TAU     = 2.0     # speed of SP response

    # Solver options
    m.options.IMODE = 6     # Dynamic Optimization (Control)
    m.options.CV_TYPE = 2   # L1 or L2 Norm

    # Solve the model
    m.solve(disp=False)

    #%% Display the results
    plt.figure()

    plt.subplot(2,1,1)
    plt.plot(m.time,u.value,'k-',label=r'$u$')
    plt.legend(loc='best')
    plt.ylabel('MV')
    plt.subplot(2,1,2)
    plt.plot(m.time,y.value,'r--',label=r'$y$')
    plt.plot(m.time,x.value,'g--',label=r'$x$')
    plt.legend(loc='best')
    plt.ylabel('CV')
    plt.xlabel('time')
    plt.show()
"""

#%% Define a standard test suite for this function
def test ():
    global options, args
    # TODO: Do something more interesting here...
    print ('Hello from the test() function!')

#%% This is only run whne script is executed as a standalone program
if __name__ == '__main__':
    import sys, os, traceback, argparse
    import time
    import re
    #from pexpect import run, spawn

    try:
        start_time = time.time()
        #parser = argparse.ArgumentParser(description="This is the program description",  usage=globals()['__doc__'])
        parser = argparse.ArgumentParser(description='This is the program description')
        parser.add_argument('--version', action='version', version='%(prog)s v'+__version__)
        parser.add_argument ('-v', '--verbose', action='store_true', help='produce verbose output')
        parser.add_argument ('-t', '--test', action='store_true', help='run test suite')
        args = parser.parse_args()
        #if len(args) < 1:
        #    parser.error ('missing argument')
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
