#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# process.py - A wrapper for the drilling process
__version__   = "1.0.2"
__date__      = "2019.04.01"
__author__    = "Robert B. Hawkins <braynebuddy@gmail.com>"

"""
SYNOPSIS

    TODO process [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    TODO This describes how to use this script. This docstring
    will be printed by the script if there is an error or
    if the user requests help (-h or --help).

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
# 1.0.2     2019-04-01  Drilling process function
#

#%% The main function in this module
def process (level_start, rho_start, depth_start, pumpQ, backQ, chokeVP, \
             mudQ, waterQ, dTime):
    """
    Arguments:
        level_start - Mud pit level at start (m)
        rho_start   - Mud density in pit at start (kg/m3)
        depth_start - Measured depth at the start of this interval (m)
        pumpQ       - Mud pump flowrate for this interval (m3/min)
        backQ       - Back-pressure pump flowrate for this interval (m3/min)
        chokeVP     - Choke Valve Opening from 0-100 [%]
        mudQ        - Mud concentrate makeup (m3/min)
        waterQ      - Fresh water makeup (m3/min)
        dTime       - Length of the interval (seconds)
    Returns:
        level    - Mud pit level at end (m)
        rho      - Mud density in pit at end (kg/m3)
        depth    - Measured depth at end of interval [m]
    """

    #import numpy as np
    #from gekko import GEKKO
    #import matplotlib.pyplot as plt
    
    import drillingmodel as drill
    import mudpitmodel as mud

    global options, args

    print('\tCalling drillstring()...')
    Pp, Pc, Qb, Pb, Qc, Qr, depth = drill.drillstring(pumpQ, backQ, \
                                         chokeVP, depth_start, dTime*60 )

    print('\t-[Drilling Results]-----------------------------')
    print('\tPressure at Pump [bar] =', Pp)
    print('\tPressure at Choke Valve [bar] =', Pc)
    print('\tBit pressure [bar] =', Pb)
    print('\tFlow Rate through Bit [m^3/min] =', Qb)
    print('\tFlow Rate through Choke [m^3/min] =', Qc)
    print('\tFlow Rate from reservoir [m^3/min] =', Qr)
    print('\tFlow Rate from reservoir [bbl/h] =', Qr*60/0.159)
    print('\t------------------------------------------------\n')

    rhoC = rho_start
    
    print ('\tCalling mudpit()\n')
    level, rho = mud.mudpit(level_start, rho_start, rhoC, Qc, pumpQ+backQ, \
                        mudQ, waterQ, dTime*60)

    return (level, rho, depth)

#%% Define a stand-alone call for the function
def main():
    level = 10.0
    rhoP = 1240
    #meas_depth = 6000.0*.3048  # Normal
    #meas_depth = 8000.0*.3048  # loss of circ
    meas_depth = 15400.0*.3048  # Kick

    mud_pump_flow = 1.0 #2004.0*(1e-3/60)
    bp_pump_flow = 0.2 #804.0*(1e-3/60)
    choke_valve = 50.0

    mudQ = 0.0
    waterQ = 0.0
    
    time_interval = 5.0

    print('-[Inputs]---------------------------------------')
    print('Initial mud pit level   =', level, 'm')
    print('Initial mud density     =', rhoP, 'kg/m3')
    print('Initial measured depth  =', meas_depth, 'm')
    print('------------------------------------------------')
    print('Mud Pump Flow           =', mud_pump_flow, 'm3/min')
    print('Back-pressure Pump Flow =', bp_pump_flow, 'm3/min')
    print('Choke Valve Position    =', choke_valve, '%')
    print('------------------------------------------------')
    print('Makeup Mud Flow         =', mudQ, 'm3/min')
    print('Makeup Water Flow       =', waterQ, 'm3/min')
    print('------------------------------------------------')
    print('Time interval           =', time_interval, 'minutes')
    print('------------------------------------------------\n')

    print('Calling process...')
    level, rhoP, meas_depth = process(level, rhoP, meas_depth, \
                                      mud_pump_flow, bp_pump_flow, \
                                      choke_valve, \
                                      mudQ, waterQ, \
                                      time_interval*60 )
    
    print('-[Results]--------------------------------------')
    print('Mud pit level   =', level, 'm')
    print('Mud density     =', rhoP, 'kg/m3')
    print('Measured depth  =', meas_depth, 'm')
    print('------------------------------------------------\n')
    
    return
    
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
