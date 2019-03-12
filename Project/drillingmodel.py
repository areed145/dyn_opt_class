#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import math
import numpy as np
#import reservoirmodel as rm
from gekko import GEKKO


rmt = True  # Solve local or remote

#%% Drillstring Model
def drillstring(pumpQ, backQ, depth, dTime):
    
    
    
    # Specify Model
    
    #How to define the number of elements?
    #Temperature change? Depth? 
    N = 3                 #Number of elements along the borehole, min = 3
    
    
    nSteps = 10
    d = GEKKO(remote=rmt)
    d.time = np.linspace(0, dTime, nSteps+1)

    # Reservoir Model
    #TVD, ROP, Pf = rm.reservoir(depth)
    TVD = 1900
    ROP = 0
    Pf = 1500

    # Model Constants
    g = 9.81    # Gravitational constant (kg-m/s^2)
    #h = 1951.0  # Height of wellbore (m)
    r_di = 2.0 * 0.0254    # drillstring inner radius (m) (4" diameter)
    r_do = 2.5 * 0.0254    # drillstring outer radius (m) (5" diameter)
    r_ci = 4.3125 * 0.0254 # annulus inner radius (m) (8 5/8" diameter)
    A_d = math.pi*r_di**2  # drillstring inner area , m^2
    A_a = math.pi*(r_ci**2 - r_do**2)  # annulus flow area , m^2
    A_h = math.pi*r_ci**2  # borehole cross area, m^2
    
    Fa = d.Const(value = 2.0e+9)   # friction coefficient inside the annulus
    Fd = d.Const(value = 5.0e+9)   # friction coefficient inside the drill string
    Fb = d.Const(value = 1.0e+9)   # friction coefficient inside the drill bit
    
    betaD = d.Const(value = 2.0e+9, name='betad')   # bulk modulus of the drill string
    betaA = d.Const(value = 1.0e+9, name='betaA')   # bulk modulus of the annulus

    #Ma    = d.FV(value = 0.0)      # mass of fluid in the annulus
    #Md    = d.FV(value = 0.0)      # mass of fluid in the drill string
    mass  = d.FV(value = 4.3e+8, name='mass')   # Total mass of fluid (Ma + Md)

    rhoA  = d.FV(value = 1580.0)   # density of the fluid in the drill string
    rhoD  = d.FV(value = 1580.0)   # density of the fluid in the annulus
    #rhoM  = d.FV(value = 1580.0)   # density of the mud used

    Qres  = d.Param(value = 0.0, name='Qres')      # volume flow rate from the reservoir
    Qloss  = d.Param(value = 0.0, name='Qloss')      # volume flow rate loss to the reservoir
    
    # Define Variables
    Qpump = d.FV(value = pumpQ, name='Qpump')  # volume flow rate through the pump
    Qback = d.FV(value = backQ, name='Qback')  # volume flow rate through the back pressure pump
    
    Vd_init = A_d * depth
    Vd      = d.Var(value = Vd_init, name='Vd')  # volume of the drill string (m3)
    
    Va_init = A_a * depth    
    Va      = d.Var(value = Va_init, name='Va')  # volume of the annulus (m3)
    
    Pp    = d.Var(value = 40.0e+5, name='Pp')  # pump pressure inside the drill string
    Pc    = d.Var(value = 10.0e+5, name='Pc')  # choke pressure inside the annulus
    Qbit  = d.Var(value = pumpQ, name='Qbit')   # volume flow rate through the drill bit
     # volume flow rate through the choke (Qpump+Qback)
    #Qchoke = d.Intermediate( Qbit + Qback + Qres - Qloss)
    Qchoke = d.Var(value = 0, name='Qchoke')
    d.Equation( Qchoke == Qbit + Qback + Qres - Qloss )

    MD = d.Var(value = depth, name='MD')    # measured depth
    

    Pbit  = d.Var(value = 0.0, name='Pbit')      # Pressure at the drill bit / BHP
    #dPa   = m.Var(value = 0.0)      # differential pressure across the annulus (?)
        
    Te = 59.5   # surface earth temperature (degF)
    ge = 0.025  # geothermal gradient (degF/ft)
    
    # Model equations
    
    # Drillstring and annulus volumes
    d.Equation ( Vd.dt() == A_d * ROP )
    d.Equation ( Va.dt() == A_a * ROP )
    
    # Equation 5.1
    d.Equation( Pp.dt() == (betaD/Vd) * (Qpump - Qbit - Vd.dt()) )
    #d.Equation( Vd*Pp.dt() == (betaD) * (Qpump - Qbit) )
    
    # Equation 5.2
    d.Equation( Pc.dt() == (betaA/Va) * (Qres + Qbit + Qback - Qchoke - Va.dt()) )
    #d.Equation( Va*Pc.dt() == (betaA) * (Qres + Qbit + Qback - Qchoke) )
    
    # Equation 5.3
    #Using the equation by replacing dPa by the pressure balance in drillpipe
    d.Equation( Qbit.dt() == (1/mass) * (Pp - Pbit - Fd*(Qbit**2) \
                        + rhoD*g*TVD ) )
    
    # Equation 5.6
    d.Equation( Pbit == Pc + rhoA*g*TVD + Fa*(Qbit**2) )
    
    # Drilling rate from reservior simulation
    d.Equation( MD.dt() == ROP )
    
    
    ###########################################################################
    ##Finite elements approach
    ###########################################################################
    
    #Need to complete this. But the flow of reading the model is to imagine
    #the drillstring and the annulus being divided into N sections
    #Top of the well is the start for the drillstring element 0 -> bottom is element N-1
    #Bottom of the well is the start for the annulus element 0 -> wellhead is element N-1
    #The equations are in the sequential order of flow across each element inside the 
    #drillstring, across the bit and finally up the annulus
    qp = 2000.0*1e-3/60
    pp = 40.0e+5
    
    Qd = [d.Var(value=0) for i in range(N)]     #Flow inside the drillstring, m^3/s
    Qa = [d.Var(value=0) for i in range(N)]     #Flow inside the annulus, m^3/s
    Pd = [d.Var(value=0) for i in range(N)]     #Pressure inside drillstring, Pa
    Pa = [d.Var(value=0) for i in range(N)]     #Pressure inside annulus, Pa
    
    #Pressure inside the drillstring
    d.Equation(Pd[0].dt() == (betaD/(Vd/N))*(Qpump - Qd[0]))    
    for i in range(1,N-1):
        d.Equation(Pd[i].dt() == (betaD/(Vd/N))*(Qd[i-1] - Qd[i]))    
    d.Equation(Pd[N-1].dt() == (betaD/(Vd/N))*(Qd[N-2] - Qbit))
    
    #Flow inside the drillstring    
    for i in range(0,N-1):
        d.Equation(Qd[i].dt() == (1/(mass/N)) * (Pd[i] - Pd[i+1] - (Fd/N)*(Qd[i]**2) \
                                    + rhoD*g*TVD/N ))   

         
    #Flow across the drillbit
    d.Equation(Qbit.dt() == (1/(mass/N)) * (Pd[N-1] - Pa[0] - Fb*(Qbit**2) \
                        + (rhoD - rhoA)*g*TVD/N) )  
    
    
    #Pressure inside the annulus
    d.Equation(Pa[0].dt() == (betaA/(Va/N))*(Qbit + Qres - Qloss - Qa[0]))
    for i in range(1,N):
        d.Equation(Pa[i].dt() == (betaA/(Va/N))*(Qa[i-1] - Qa[i]))
    d.Equation(Pa[N-1].dt() == (betaA/(Va/N))*(Qa[N-2] - Qchoke + Qback))
    
    #Flow inside the annulus
    for i in range(0,N-1):
        d.Equation(Qa[i].dt() == (1/(mass/N)) * (Pa[i] - Pa[i+1] - (Fa/N)*(Qa[i]**2) \
                        - rhoA*g*TVD/N) ) 
    
    
    
    
    
    
    
    
    
    # @@ Needed to complete simulation...
    
    # Calculate Qres from reservior pressure and Pbit
    # Add temperature effects on mud density from TVD
    # Finite elements equations replacing 5.1-5.6
    
    
    
    
    
    
    #d.options.solver = 1
    # Solver options
    d.options.IMODE = 4  # Dynamic Simulation
    d.options.MAX_ITER = 5000
    
    
    # Solve the model
    d.solve()

    return (Pp.VALUE[-1], Pc.VALUE[-1], \
            Qbit.VALUE[-1], Pbit.VALUE[-1], \
            Qres.VALUE[-1], MD.VALUE[-1] )
    
print(drillstring(2000.0*(1e-3/60), 800.0*(1e-3/60), 1900,10))