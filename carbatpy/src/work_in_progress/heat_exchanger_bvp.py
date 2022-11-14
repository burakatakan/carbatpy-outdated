# -*- coding: utf-8 -*-
""" heat exchanger from ODE solving (boundary value problem)

using coolprop for property evaluation
formulated as an energy balance differential equation
T is evaluated for each enthalpy
then the convective heat transfer is calculated
Created on Thu Jan 30 13:12:21 2020

@author: atakan
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
import CoolProp.CoolProp as CP
from coolPropTGlide0 import tp, hps

mdot=np.array((.0029711, .0351)) # kg/s for both fluids
alpha = 500  # heat transfer coefficient through the wall (from total resistance)
r = 13e-2   # radius tube , m
L = 1.    # length tube, m
U = 2 * np.pi * r  # circimference tube
p = [6.545e5, 4.e5]  # pressure for each fluid, Pa
fl =["ISOBUTANE","Water"]   # which fluids?
Tin = [354, 313]
ha_in=tp(Tin[0], p[0], fl[0])[2]  # state of fluid 1 left
hb_in=tp([Tin[1]],p[1],fl[1])[2]  # state of fluide 2 right (at L)
ha_outMax = tp(Tin[1],p[0],fl[0])[2]  # state of fluid 1 left
hb_outMax = tp(Tin[0],p[1],fl[1])[2]  # state of fluide 2 right (at L)
npo = 147                         # nuber of initial points
print(mdot[0]*(ha_in-ha_outMax))


def energy(x,h): 
    """
    energy balance 
    spatial coordinate and both enthalpies as input
    output: both changes in enthalpy in x-direction
    
    depends on global variables: 
        fl (fluid-names)
        alpha: convection coefficient W/m2/K
        mdot: mass flow rates kg/s
        U: circumference of tube m
    function hps returns an array, the first value is temperature
    """
  
    q_konv=(hps(h[1],p[1],fl[1])[0]-hps(h[0],p[0],fl[0])[0])
    # problem = np.where(q_konv>0)
    # q_konv[problem] = -1e-1
    #print(q_konv,h,hps(h[1],p[1],fl[1])[0],hps(h[0],p[0],fl[0])[0])
    dh0 = alpha*U/mdot[0]*q_konv
    dh1 = alpha*U/mdot[1]*(q_konv)
    return np.array([dh0,dh1])


def bc(ha,hb): #boundary conditions
    # return np.array([ha[0]-ha_in,hb[1]-hb_in,])
    return np.array([ha[0]-ha_in,hb[1]-hb_in,])

x=np.linspace(0,L,npo)
y=np.zeros((2,npo))
dh0=ha_in-ha_outMax
y[0,:]=np.linspace(ha_in,ha_outMax,npo)
y[1,:]=np.linspace(hb_in,hb_in,npo)

result=solve_bvp(energy,bc,x,y, tol=5e-3, max_nodes=1000)
if result.success:
    T0=hps(result.y[0],p[0],fl[0])[0]
    T1=hps(result.y[1],p[1],fl[1])[0]
    s0=hps(result.y[0],p[0],fl[0])[4]
    s1=hps(result.y[1],p[1],fl[1])[4]
    plt.figure()
    plt.plot(result.x,T0)
    plt.plot(result.x,T1)
    plt.show()
    ds = (s0[0]-s0[-1]) * mdot[0] + (s1[0]-s1[-1])*mdot[1]
    print("Entropieproduktion:%3.2f" % (ds))
else: 
    print("Fehler, keine Lösung!", result.message)
    T0=hps(result.y[0],p[0],fl[0])[0]
    T1=hps(result.y[1],p[1],fl[1])[0]
    plt.figure()
    plt.plot(result.x,T0)
    plt.plot(result.x,T1)
    plt.show()