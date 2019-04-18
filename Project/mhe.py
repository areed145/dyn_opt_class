#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# mhe.py Moving horizon estimator for drilling problem
#
__version__   = "1.0.2"
__date__      = "2019.04.02"

"""
SYNOPSIS

    mhe [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    This is the moving horizon estimator (MHE) for the 
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
# 1.0.2     2019-04-02  First version of MHE function
#

import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

def mhe(dTime):
    """
    Arguments:
        dTime       - Length of the interval (seconds)
    Returns:
        pitA     - Mud pit area estimate (m2)
                 - Choke valve constant
    """

    #%% Non-model parameters
    rmt = True  # Solve local or remote
    npts = 100    # time steps

    #%% Specify model
    #m = GEKKO(remote=rmt)
    #m.time = np.linspace(0, dTime, npts+1)  # Model Timeline

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

    return (pitA)


#%% The main function in this module
def main ():

    global options, args


    #%% Display the results
    #plt.figure()

    #plt.subplot(2,1,1)
    #plt.plot(m.time,u.value,'k-',label=r'$u$')
    #plt.legend(loc='best')
    #plt.ylabel('MV')
    #plt.subplot(2,1,2)
    #plt.plot(m.time,y.value,'r--',label=r'$y$')
    #plt.plot(m.time,x.value,'g--',label=r'$x$')
    #plt.legend(loc='best')
    #plt.ylabel('CV')
    #plt.xlabel('time')
    #plt.show()

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
        parser = argparse.ArgumentParser(description='Drilling problem MHE')
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
