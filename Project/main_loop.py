#
# main_loop.py - The main program loop for the MPD project
"""
SYNOPSIS

    main_loop.py [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    TODO This describes how to use this script. This docstring
    will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    Probably no exit codes from this routine. :)

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
__version__   = "1.0.0"
__date__      = "2019-04-02"

# Version   Date        Notes
# -------   ----------  -------------------------------------------------------
# 1.0.0     2013.12.01  Starting script template
#           2019-04-02  Managed pressure drilling project mainline
#

import sys, os, traceback, argparse
import time
import re
#from pexpect import run, spawn

def main ():

    global options, args
    # TODO: Do something more interesting here...
    print ('Welcome to the MPD Project main loop!')

def test ():

    global options, args
    # TODO: Do something more interesting here...
    print ('Hello from the test() function!')

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
