#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# reservoirmodel.py - A simple table-lookup model of an oil reservoir
#
__version__   = "1.0.0"
__date__      = "2019.03.18"

import pandas as pd
import numpy as np
#%% Reservoir Model
res_data_xls = pd.read_excel(r'example well pressure and flow function of depth and time.xlsx', \
                             skiprows=10)
res_data_xls.drop(0, inplace=True)
res_data_xls.sort_values(by=['MD Measured Depth (ft)'],inplace=True)

res_data_xls_MD  = res_data_xls['MD Measured Depth (ft)'].values  # measured depth
res_data_xls_TVD = res_data_xls['TVD True Vertical Depth(ft)'].values  # total vertical depth
res_data_xls_ROP = res_data_xls['ROP (m/sec)'].values   # rate of penetration
res_data_xls_K   = res_data_xls['Permeability_k (m2)'].values  # permeability
res_data_xls_PF  = res_data_xls['Pressure (formation bar)'].values  # formation pressure
res_data_xls_EL  = res_data_xls['equivilant_length (m)'].values  # effective length

#print ('res_MD =', (res_MD,))
#print ('res_TVD =', (res_TVD,))
#print ('res_ROP =', (res_ROP,))
#print ('res_K =', (res_K,))
#print ('res_PF =', (res_PF,))
#print ('res_EL =', (res_EL,))

def res_data():
    MD = res_data_xls_MD * 0.3048  # convert ft > m
    TVD = res_data_xls_TVD * 0.3048  # convert ft > m
    ROP = res_data_xls_ROP
    K = res_data_xls_K
    PF = res_data_xls_PF
    EL = res_data_xls_EL

    return (MD, TVD, ROP, K, PF, EL)

def reservoir_mindepth():
    return res_data_xls_MD[0]

def reservoir_flow(pressure, area, depth):
    """
    Parameters:
        pressure - drill string pressure
        area     - drill string flowing area
        depth    - measured depth of the drill string
    
    Returns:
        Qres - flow in/our of the reservoir
    """
    res_MD, res_TVD, res_ROP, res_K, res_PF, res_EL = res_data()
    k = 0.0
    formation_pressure = 0.0
    effective_length = 1.0

    idx = 0
    max_idx = len(res_MD)
    for idx in range(max_idx-1):
        if depth < res_MD[idx]:
            idx += 1
        else:
            k = res_ROP[idx]
            formation_pressure = res_PF[idx]
            effective_length = res_EL[idx]
            break

    Qres = k * area * (formation_pressure - pressure)/effective_length
    return Qres

def reservoir(depth):
    """
    Parameters:
        depth - the measured depth of the drill string
    
    Returns:
        tvd - total vertical depth
        rop - rate of penetration of the drill
        pf - pressure in the formation
    """
    res_MD, res_TVD, res_ROP, res_K, res_PF, res_EL = res_data()
    tvd = 0
    rop = 0
    pf = 0
    k = 0
    el = 1
    idx = 0
    while idx < len(res_MD)-1 and depth >= res_MD[idx]:
        idx += 1
    
    rop = res_ROP[idx]
    pf = res_PF[idx]
    k = res_K[idx]
    el = res_EL[idx]

    if (idx+1<len(res_MD)):
        tvd_slope = (res_TVD[idx+1]-res_TVD[idx]) / (res_MD[idx+1]-res_MD[idx])
    else:
        tvd_slope = (res_TVD[idx]-res_TVD[idx-1]) / (res_MD[idx]-res_MD[idx-1])
    tvd =  res_TVD[idx] + (depth-res_MD[idx]) * tvd_slope
        
    # Does the model need to return the tvd and pf derivatives as well?
    return (tvd, rop, pf, k, el)

def reservoir_TVD(depth):
    """
    Calling Argument:
        depth - the measured depth of the drill string
    
    Function Returns:
        tvd - total vertical depth
    """
    res_MD, res_TVD, res_ROP, res_K, res_PF, res_EL = res_data()
    tvd = 0
    idx = 0
    while idx < len(res_MD)-1 and depth >= res_MD[idx]:
        idx += 1
    
    if (idx+1<len(res_MD)):
        tvd_slope = (res_TVD[idx+1]-res_TVD[idx]) / (res_MD[idx+1]-res_MD[idx])
    else:
        tvd_slope = (res_TVD[idx]-res_TVD[idx-1]) / (res_MD[idx]-res_MD[idx-1])
    tvd =  res_TVD[idx] + (depth-res_MD[idx]) * tvd_slope
        
    return tvd

def reservoir_dTVD(depth):
    """
    Calling Argument:
        depth - the measured depth of the drill string
    
    Function Returns:
        tvd - derivative of total vertical depth
    """
    res_MD, res_TVD, res_ROP, res_K, res_PF, res_EL = res_data()
    tvd_slope = 1.0
    idx = 0
    while idx < len(res_MD)-1 and depth >= res_MD[idx]:
        idx += 1
    
    if (idx+1<len(res_MD)):
        tvd_slope = (res_TVD[idx+1]-res_TVD[idx]) / (res_MD[idx+1]-res_MD[idx])
    else:
        tvd_slope = (res_TVD[idx]-res_TVD[idx-1]) / (res_MD[idx]-res_MD[idx-1])
    return tvd_slope

def reservoir_ROP(depth):
    """
    Calling Argument:
        depth - the measured depth of the drill string
    
    Function Returns:
        rop - rate of penetration of the drill
    """
    res_MD, res_TVD, res_ROP, res_K, res_PF, res_EL = res_data()
    idx = 0
    while idx < len(res_MD)-1 and depth >= res_MD[idx]:
        idx += 1
    return res_ROP[idx]

def test(pressure, area, depth):
    print('Reservoir variables:')
    tvd, rop, pf, tvds = reservoir(depth)

    print ('measured depth =', depth)
    print ('total vertical depth =', tvd)
    print ('total vertical depth slope =', tvds)
    print ('rate of penetration =', rop)
    print ('formation pressure =', pf)
    
    print ('flow from reservoir =', reservoir_flow(pressure, area, depth))

def main():
    print('     Depth      TVD  slope(TVD)    ROP     Form Pr  Res Flow')
    for d in range(1000, 30000, 1000):
        tvd, rop, pf, tvds = reservoir(d)
        print("{0:10d}{1:10.1f}{2:10.5f}{3:10.5f}{4:10.1f}{5:10.1f}".format( \
              d, reservoir_TVD(d), reservoir_dTVD(d), reservoir_ROP(d), \
              pf, reservoir_flow(1800, 0.11, d) ) )

#%% This is only run whne script is executed as a standalone program
if __name__ == '__main__':
    import sys, os, traceback, argparse
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
            test(1850, 4.3125*0.0254, 4000)
        else:
            main()
        if args.verbose: print (time.asctime())
        if args.verbose: print ('Elapsed time:', \
                                (time.time() - start_time), 'seconds')
        sys.exit(0)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print ('ERROR, UNEXPECTED EXCEPTION')
        print (str(e))
        traceback.print_exc()
        os._exit(1)
