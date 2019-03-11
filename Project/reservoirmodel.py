#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

#%% Reservoir Model
res_data = pd.read_csv('reservoir_v1.csv')
res_MD  = res_data['MD'].values     # measured depth
res_TVD = res_data['TVD'].values    # total vertical depth
res_ROP = res_data['ROP'].values    # rate of penetration
res_PF  = res_data['Pf'].values     # formation pressure

def reservoir(depth):
    tvd = 0
    rop = 0
    pf = 0
    for i in range(len(res_MD)):
        if (depth >= res_MD[i]):
            if (i+1<len(res_MD)):
                tvd_slope = (res_TVD[i+1]-res_TVD[i]) / (res_MD[i+1]-res_MD[i])
            else:
                tvd_slope = (res_TVD[i]-res_TVD[i-1]) / (res_MD[i]-res_MD[i-1])
        tvd =  res_TVD[i] + (depth-res_MD[i]) * tvd_slope
        rop = res_ROP[i]
        pf = res_PF[i]
        
    # Does the model need to return the tvd and pf derivatives as well?
    return (tvd, rop, pf)