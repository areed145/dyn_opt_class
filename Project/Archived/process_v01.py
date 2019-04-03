#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# TODO prog_base.py - A starting template for Python scripts
#
# Copyright 2013 Robert B. Hawkins
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

AUTHOR

    Rob Hawkins <webwords@txhawkins.net>

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

VERSION

    1.0.0
"""
__author__    = "Rob Hawkins <webwords@txhawkins.net>"
__version__   = "1.0.0"
__date__      = "2013.12.01"

# Version   Date        Notes
# -------   ----------  -------------------------------------------------------
# 1.0.0     2013.12.01  Starting script template
#

#%%Import packages
import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt
import sys, os, traceback, argparse
import time
import re
#from pexpect import run, spawn

def process():

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

    return

def test ():

    global options, args
    # TODO: Do something more interesting here...
    print ('Hello from the test() function!')

def main ():

    global options, args
    # TODO: Do something more interesting here...
    process()

if __name__ == '__main__':
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
        if args.verbose: print ('TOTAL TIME IN MINUTES:',end="")
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
