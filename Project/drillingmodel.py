#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import math
import numpy as np
import reservoirmodel as rm
from gekko import GEKKO


rmt = True  # Solve local or remote

#%% Drillstring Model
def drillstring(pumpQ, backQ, depth, dTime):
    
    
    
    # Specify Model
    
    #How to define the number of elements?
    #Temperature change? Depth? 
    N = 10                 #Number of elements along the borehole, min = 3
    
    
    nSteps = 10
    d = GEKKO(remote=rmt)
    d.time = np.linspace(0, dTime, nSteps+1)

    # Reservoir Model
    TVD, ROP, Pf = rm.reservoir(depth)

    # Model Constants
    g = 9.81    # Gravitational constant (kg-m/s^2)
    #h = 1951.0  # Height of wellbore (m)
    r_di = 2.0 * 0.0254    # drillstring inner radius (m) (4" diameter)
    r_do = 2.5 * 0.0254    # drillstring outer radius (m) (5" diameter)
    r_ci = 4.3125 * 0.0254 # annulus inner radius (m) (8 5/8" diameter)
    
    Fa = d.Const(value = 2.0e+9)   # friction coefficient inside the annulus
    Fd = d.Const(value = 5.0e+9)   # friction coefficient inside the drill string
    Fb = d.Const(value = 1.0e+9)   # friction coefficient inside the drill bit
    
    betaD = d.Const(value = 2.0e+9)   # bulk modulus of the drill string
    betaA = d.Const(value = 1.0e+9)   # bulk modulus of the annulus

    #Ma    = d.FV(value = 0.0)      # mass of fluid in the annulus
    #Md    = d.FV(value = 0.0)      # mass of fluid in the drill string
    mass  = d.FV(value = 4.3e+8)   # Total mass of fluid (Ma + Md)

    rhoA  = d.FV(value = 1580.0)   # density of the fluid in the drill string
    rhoD  = d.FV(value = 1580.0)   # density of the fluid in the annulus
    #rhoM  = d.FV(value = 1580.0)   # density of the mud used

    Qres  = d.FV(value = 0.0)      # volume flow rate from the reservoir
    Qloss  = d.FV(value = 0.0)      # volume flow rate loss to the reservoir
    
    # Define Variables
    Qpump = d.FV(value = pumpQ)  # volume flow rate through the pump
    Qback = d.FV(value = backQ)  # volume flow rate through the back pressure pump
    
    Vd_init = math.pi*r_di*r_di * depth
    Vd      = d.Var(value = Vd_init)  # volume of the drill string (m3)
    
    Va_init = math.pi*(r_ci*r_ci - r_do*r_do) * depth    
    Va      = d.Var(value = Va_init)  # volume of the annulus (m3)
    
    Pp    = d.Var(value = 40.0e+5)  # pump pressure inside the drill string
    Pc    = d.Var(value = 10.0e+5)  # choke pressure inside the annulus
    Qbit  = d.Var(value = pumpQ)   # volume flow rate through the drill bit
     # volume flow rate through the choke (Qpump+Qback)
    Qchoke = d.Intermediate( Qbit + Qback )

    MD = d.Var(value = depth)    # measured depth
    

    Pbit  = d.Var(value = 0.0)      # Pressure at the drill bit / BHP
    #dPa   = m.Var(value = 0.0)      # differential pressure across the annulus (?)
        
    Te = 59.5   # surface earth temperature (degF)
    ge = 0.025  # geothermal gradient (degF/ft)
    
    # Model equations
    
    # Drillstring and annulus volumes
    d.Equation ( Vd == math.pi*r_di*r_di * MD )
    d.Equation ( Va == math.pi*(r_ci*r_ci - r_do*r_do) * MD )
    
    # Equation 5.1
    #d.Equation( Pp.dt() == (betaD/Vd) * (Qpump - Qbit - Vd.dt()) )
    #d.Equation( Vd*Pp.dt() == (betaD) * (Qpump - Qbit) )
    
    # Equation 5.2
    #d.Equation( Pc.dt() == (betaA/Va) * (Qres + Qbit + Qback - Qchoke - Va.dt()) )
    #d.Equation( Va*Pc.dt() == (betaA) * (Qres + Qbit + Qback - Qchoke) )
    
    # Equation 5.3
    #Using the equation by replacing dPa by the pressure balance in drillpipe
    d.Equation( Qbit.dt() == (1/mass) * (Pp - Pbit - Fd*(Qbit**2) \
                        + rhoD*g*TVD ) )
    
    # Equation 5.6
    d.Equation( Pbit == Pc + rhoA*g*TVD + Fa*(Qbit**2) )
    
    # Drilling rate from reservior simulation
#    d.Equation( TVD.dt() == ROP )
    
    
    ###########################################################################
    ##Finite elements approach
    ###########################################################################
    
    #Need to complete this. But the flow of reading the model is to imagine
    #the drillstring and the annulus being divided into N sections
    #Top of the well is the start for the drillstring element 1 -> botton is element N
    #Bottom of the well is the start for the annulus element 1 -> wellhead is element N
    #The equations are in the sequential order of flow across each element inside the 
    #drillstring, across the bit and finally up the annulus
    
    
    Qd = [d.Var(value=0) for i in range(N)]     #Flow inside the drillstring, m^3/s
    Qa = [d.Var(value=0) for i in range(N)]     #Flow inside the annulus, m^3/s
    Pd = [d.Var(value=0) for i in range(N)]     #Pressure inside drillstring, Pa
    Pa = [d.Var(value=0) for i in range(N)]     #Pressure inside annulus, Pa
    
    #Pressure inside the drillstring
    d.Equation(Pd[0].dt() == (betaD/Vd)*(Qpump - Qd[0]))    
    for i in range(1,N-1):
        d.Equation(Pd[i].dt() == (betaD/Vd)*(Qd[i-1] - Qd[i]))    
    d.Equation(Pd[N-1].dt() == (betaD/Vd)*(Qd[N-2] - Qbit))
    
    #Flow inside the drillstring    
    for i in range(0,N-1):
        d.Equation(Qd[i].dt() == (1/mass) * (Pd[i] - Pd[i+1] - Fd*(Qd[i]**2) \
                                    + rhoD*g*TVD ))   

         
    #Flow across the drillbit
    d.Equation(Qbit.dt() == (1/mass) * (Pd[N-1] - Pa[0] - Fb*(Qbit**2) \
                        + 0.5*(rhoD - rhoA)*g*TVD) )  
    
    
    #Pressure inside the annulus
    d.Equation(Pa[0].dt() == (betaA/Va)*(Qbit + Qres - Qloss - Qa[0]))
    for i in range(1,N):
        d.Equation(Pa[i].dt() == (betaA/Va)*(Qa[i-1] - Qa[i]))
    d.Equation(Pa[N-1].dt == (betaA/Va)*(Qa[N-2] - Qchoke + Qback))
    
    #Flow inside the annulus
    for i in range(0,N-1):
        d.Equation(Qa[i].dt() == (1/mass) * (Pa[i] - Pa[i+1] - Fa*(Qa[i]**2) \
                        - 0.5*rhoA*g*TVD) ) 
    
    
    
    
    
    
    
    
    
    # @@ Needed to complete simulation...
    
    # Calculate Qres from reservior pressure and Pbit
    # Add temperature effects on mud density from TVD
    # Finite elements equations replacing 5.1-5.6
    
    
    
    
    
    
    
    # Solver options
    d.options.IMODE = 4  # Dynamic Simulation
    
    
    # Solve the model
    d.solve()

    return (Pp.VALUE[-1], Pc.VALUE[-1], \
            Qbit.VALUE[-1], Pbit.VALUE[-1], \
            Qres.VALUE[-1], MD.VALUE[-1] )
    
print(drillstring(86e+5, 10e+5,1900,10))