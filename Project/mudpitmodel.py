#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# TODO prog_base.py - A starting template for Python scripts
#
"""
SYNOPSIS

    TODO prog_base [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    TODO This describes how to use this script. This docstring
    will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

LICENSE

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""

# Version   Date        Notes
# -------   ----------  -------------------------------------------------------
# 1.0.0     2013.12.01  Starting script template
#
__version__   = "1.0.0"
__date__      = "2013.12.01"

#%%Import packages
import numpy as np
from gekko import GEKKO
#import matplotlib.pyplot as plt

def mudpit(level_st, rho_st, rho_in, inflow, outflow, mud, water, dTime):
    """
    Arguments:
        level_st - Mud pit level at start (m)
        rho_st   - Mud density in pit at start (kg/m3)
        rho_in   - Density of inflow from choke (kg/m3)
        inflow   - Flow in from Choke (m3/min)
        outflow  - Flow out to mud pumps (m3/min)
        mudflow  - Mud concentrate makeup (m3/min)
        water    - Fresh water makeup (m3/min)
        dTime    - Length of the interval (seconds)
    Returns:
        level    - Mud pit level at end (m)
        rho      - Mud density in pit at end (kg/m3)
    """
    #%% Non-model parameters
    rmt = False  # Solve local or remote

    #%% Specify model
    m = GEKKO(remote=rmt)

    st = 10.0   # simulation time interval (seconds)
    nt = int(dTime/st)+1 # simulation time points
    m.time = np.linspace(0,dTime,nt)

    # Model Constants and Parameters
    rhoM = m.Const(value=1800.0)    # makeup mud density (kg/m3)
    rhoW = m.Const(value=1000.0)    # water density (kg/m3)
    #pitV = m.Const(value=10.0)     # Mud pit volume (m3)
    pitA = m.Const(value=5.0)      # Mud pit area (m2)

    # Model Parameters

    # Model Variables
    level = m.Var(value=level_st)                 # Pit level (m)
    rho   = m.Var(value=rho_st)                   # Mud density (kg/m3)
    massP = m.Var(value=level_st * pitA * rho_st) # Mud mass in pit (kg)
    volP = m.Var(value=level_st * pitA)           # Volume of pit (m3)

    # Model Equations
    m.Equation( level.dt() == (1/pitA)*(inflow + mud + water - outflow)/60 )
    m.Equation( massP.dt() == (rho_in*inflow + rhoM*mud + rhoW*water - rho*outflow)/60)
    m.Equation( volP.dt() == (inflow + mud + water - outflow)/60)
    m.Equation( rho == massP / volP)

    # Solver options
    m.options.IMODE = 4     # Dynamic Simulation

    # Solve the model
    m.solve(disp=False)

    # Return final state
    return (level.VALUE[-1], rho.VALUE[-1])



def test ():

    global options, args
    # TODO: Do something more interesting here...
    print ('Hello from the test() function!')

def main ():

    global options, args
    mud_pump_flow = 1.0 #2004.0*(1e-3/60)
    bp_pump_flow = 0.2 #804.0*(1e-3/60)
    choke_flow = 1.19601
    mud = 0.0       # Makeup mud flow
    water = 0.0     # Makeup water
    level = 10      # mud pit level (m)
    rhoC = 1240      # choke flow density
    rhoP = 1240      # mud pit density

    time_interval = 5.0    # time step (minutes)

    outflow = mud_pump_flow + bp_pump_flow

    print('-[Inputs]---------------------------------------')
    print('Mud Pump Flow           =', mud_pump_flow, 'm3/min')
    print('Back-pressure Pump Flow =', bp_pump_flow, 'm3/min')
    print('Choke Flow              =', choke_flow, 'm3/min')
    print('Makeup Mud Flow         =', mud, 'm3/min')
    print('Makeup Water Flow       =', water, 'm3/min')
    print('------------------------------------------------')
    print('Time interval           =', time_interval, 'minutes')
    print('Initial mud pit level   =', level, 'm')
    print('Initial mud density     =', rhoP, 'kg/m3')
    print('------------------------------------------------')
    print('Expected mud pit level   =', level+time_interval*(choke_flow-outflow+mud+water)/5, 'm')
    print('------------------------------------------------\n')

    print ('calling mudpit()')
    level, rhoP = mudpit(level, rhoP, rhoC, choke_flow, outflow, mud, water, time_interval*60)

    print('-[Results]--------------------------------------')
    print('Final mud pit level   =', level, 'm')
    print('Final mud density     =', rhoP, 'kg/m3')
    print('------------------------------------------------\n')


if __name__ == '__main__':
    #import sys
    import os, traceback, argparse
    import time
    #import re
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
        if args.verbose: print ('TOTAL TIME IN MINUTES:',)
        if args.verbose: print ((time.time() - start_time) / 60.0)

    except KeyboardInterrupt as e: # Ctrl-C
        raise e

    except SystemExit as e: # sys.exit()
        raise e

    except Exception as e:
        print ('ERROR, UNEXPECTED EXCEPTION')
        print (str(e))
        traceback.print_exc()
        os._exit(1)
