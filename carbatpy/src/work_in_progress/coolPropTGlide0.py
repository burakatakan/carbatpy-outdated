# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 14:51:49 2019

@author: atakan
"""
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
fluid   = 'ISOBUTANE'
fluid2  = 'Water'
Tc1     = CP.PropsSI(fluid, "Tcrit")
Tc2     = CP.PropsSI(fluid2, "Tcrit")

pressure_at_critical_point = CP.PropsSI(fluid,'pcrit')
p_limits=np.array((1e5,10e5))
# Massic volume (in m^3/kg) is the inverse of density
# (or volumic mass in kg/m^3). Let's compute the massic volume of liquid
# at 1bar (1e5 Pa) of pressure
Tm          = 385
Tu          = 273.15+15
superheat   = 5
pinch       = 5
psatm       = CP.PropsSI('P','T',Tm,'Q',0,fluid)
vL          = 1/CP.PropsSI('D','T',Tm,'Q',0,fluid)
# Same for saturated vapor
vG          = 1/CP.PropsSI('D','T',Tm,'Q',1,fluid)
T           = CP.PropsSI('T','P',1e5,'Q',1,fluid)
print (fluid,vL,vG,T,psatm)

def tp(T,p,fluid=fluid):
    h = CP.PropsSI('H','P',p,'T',T,fluid)
    v = 1/CP.PropsSI('D','P',p,'T',T,fluid)
    s = CP.PropsSI('S','P',p,'T',T,fluid)
    x= CP.PropsSI('Q','P',p,'T',T,fluid)
    #print(T,p,h,v,s)
    return np.array([T,p,h,v,s,x])

def hps(h,p,fluid=fluid):
    
    """
    uses coolprop to evalauate all properties at a given h and p

    Parameters
    ----------
    h : enthalpy TYPE float
       J/kg.
    p : pressure TYPE float in .
    fluid :  CoolProp FluidTYPE, optional
        DESCRIPTION. The default is fluid.

    Returns
    T,p,h,v,s,x.
    -------
    TYPE numpy.array
       

    """
    
    T = CP.PropsSI('T','P',p,'H',h,fluid)
    v = 1/CP.PropsSI('D','P',p,'H',h,fluid)
    s = CP.PropsSI('S','P',p,'H',h,fluid)
    x= CP.PropsSI('Q','P',p,'H',h,fluid)
    #  print(fluid,h,p,T)
    return np.array([T,p,h,v,s,x])

def xp(x, p, fluid=fluid):
    
    """
    uses coolprop to evalauate all properties at a given h and p

    Parameters
    ----------
    h : enthalpy TYPE float
       J/kg.
    p : pressure TYPE float in .
    fluid :  CoolProp FluidTYPE, optional
        DESCRIPTION. The default is fluid.

    Returns
    T,p,h,v,s,x.
    -------
    TYPE numpy.array
       

    """
    
    T = CP.PropsSI('T','P',p,'Q',x,fluid)
    v = 1/CP.PropsSI('D','P',p,'Q',x,fluid)
    s = CP.PropsSI('S','P',p,'Q',x,fluid)
    h= CP.PropsSI('H','P',p,'Q',x,fluid)
    #  print(fluid,h,p,T)
    return np.array([T,p,h,v,s,x])

def xT(x, T, fluid=fluid):
    
    """
    uses coolprop to evalauate all properties at a given x and T

    Parameters
    ----------
    x : quality TYPE float
       .
    T : temperature TYPE float in K .
    fluid :  CoolProp FluidTYPE, optional
        DESCRIPTION. The default is fluid.

    Returns
    T,p,h,v,s,x.
    -------
    TYPE numpy.array
       

    """
    
    p = CP.PropsSI('P','T',T,'Q',x,fluid)
    v = 1/CP.PropsSI('D','T',T,'Q',x,fluid)
    s = CP.PropsSI('S','T',T,'Q',x,fluid)
    h= CP.PropsSI('H','T',T,'Q',x,fluid)
    #  print(fluid,h,p,T)
    return np.array([T,p,h,v,s,x])

def sp(s,p,fluid=fluid):
    T = CP.PropsSI('T','P',p,'S',s,fluid)
    v = 1/CP.PropsSI('D','P',p,'S',s,fluid)
    h = CP.PropsSI('H','P',p,'S',s,fluid)
    x= CP.PropsSI('Q','P',p,'S',s,fluid)
    return np.array([T,p,h,v,s,x])

def z_uv(u, v, fluid=fluid):
    T = CP.PropsSI('T', 'D', 1/v, 'U', u, fluid)
    p = CP.PropsSI('P', 'D', 1/v, 'U', u, fluid)
    h = CP.PropsSI('H', 'D', 1/v, 'U', u, fluid)
    s = CP.PropsSI('S', 'D', 1/v, 'U', u, fluid)
    return np.array([T, p, v, u, h, s])

def z_ps(p, s, fluid=fluid):
    T = CP.PropsSI('T','P',p,'S',s,fluid)
    u = CP.PropsSI('U','P',p,'S',s,fluid)
    h = CP.PropsSI('H','P',p,'S',s,fluid)
    v = 1 / CP.PropsSI('D','P',p,'S',s,fluid)
    return np.array([T, p, v, u, h, s])


def z_Tp(T, p, fluid=fluid):
    s = CP.PropsSI('S','P',p,'T',T,fluid)
    u = CP.PropsSI('U','P',p,'T',T,fluid)
    h = CP.PropsSI('H','P',p,'T',T,fluid)
    v = 1 / CP.PropsSI('D','P',p,'T',T,fluid)
    return np.array([T, p, v, u, h, s])


def z_Tx(T, x, fluid=fluid):
    s = CP.PropsSI('S','Q',x,'T',T,fluid)
    u = CP.PropsSI('U','Q',x,'T',T,fluid)
    h = CP.PropsSI('H','Q',x,'T',T,fluid)
    v = 1 / CP.PropsSI('D','Q',x,'T',T,fluid)
    p = CP.PropsSI('P','Q',x,'T',T,fluid)
    return np.array([T, p, v, u, h, s])

if __name__== "__main__":
    states=np.zeros((6,4))
    states2=np.zeros((6,4))
    #Anfangszustand, KompressoreintrittKompreesoreintritt
    states[:,0]=tp(Tm+superheat,p_limits[0])
    
    #Austritt Kompressor, isentrop
    h2=CP.PropsSI('H','P',p_limits[1],'S',states[4,0],fluid)
    T2=CP.PropsSI('T','P',p_limits[1],'S',states[4,0],fluid)
    states[:,1]=hps(h2,p_limits[1])
    Tw2=T2-pinch
    states2[:,1]=tp(Tw2,p_limits[1],fluid2)
    
    #Kondensator-Austritt
    states[:,2]=tp(Tm+5,p_limits[1])
    states2[:,2]=tp(Tm,p_limits[1],fluid2)
    
    #isenthalpe Drossel
    states[:,3]=hps(states[2,2],p_limits[0])
    states2[:,3]=tp(states[0,3]+5,p_limits[1],fluid2)
    
    dh1=states[2,2]-states[2,1]
    dh2=states2[2,2]-states2[2,1]
    print(dh1,dh2)
    m2=dh1/dh2
    nh=200
    hs1=[]
    hs2=[]
    p=p_limits[1]
    h_=np.linspace(states[2,2],states[2,1],nh)
    h2_=np.linspace(states2[2,2],states2[2,1],nh)
    
    deltah=h_[-1]-h2_[-1]*m2
    for i,h in enumerate(h_):
        hs1.append(hps(h,p,fluid))
        hs2.append(hps(h2_[i],p,fluid2))
    hs2=np.array(hs2)
    hs1=np.array(hs1)
    h_=h_-np.min(h_)    # enthalpieänderung der beiden fluide,spezif.
    h2_=h2_-np.min(h2_)
    hs1[:,4]=hs1[:,4]-np.min(hs1[:,4]) # entropieänderung der beiden fluide,spezif.
    hs2[:,4]=hs2[:,4]-np.min(hs2[:,4])
    plt.figure(1)
    plt.plot(h_,hs1[:,4],"b") 
    plt.plot(h2_*m2,hs2[:,4]*m2,".k") 
    plt.xlabel("h/(J/kg)")  
    plt.ylabel("s/(J/kg/K)")    
    plt.figure(2)
    plt.plot(h_,hs1[:,0],"b") 
    plt.plot(h2_*m2,hs2[:,0],".k")
    plt.xlabel("h/(J/kg)")  
    plt.ylabel("T/(K)")    
    
    
    
    
    