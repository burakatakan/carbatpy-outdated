# TREND 5.0 Python Interface
# Example Call
# @authors: David Celný , Sven Pohl
# Bochum, 16.05.2019

import fluid_properties_rp as fprop

# import bitzer_poly_read as br
import compressor_performance as comper
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

_props = "REFPROP"

fluid_s = "Propane"
comp = [1]
d = 0.038
h = 0.033
f = 50
V_h = np.pi * d ** 2 * h * 0.25
#secondary_fluid = CP.AbstractState("TTSE&HEOS", fluid_s)
# interesting, when using "BICUBIC&HEOS" the exergy of the ambient state is 0.15!

#statedata = fprop.T_prop_sat(Tevap, fluid_s, composition = comp, option=1)

def eta_calc(Tev, Tcond, Tin, fluid_s =fluid_s):
    """
    Calculate the isentropic work for different condenser and evaporator 
    pressures (at given T =Input) and compare it with the specific work
    as derived from a polinomial from Bitzer for their compressors.
    (From Bitzer we ger P_el and mdor , thus, h = P_el/mdot).
    The compressor inlet Temp.  is fixed at Tin
    BA 22.11.2022
    
    Parameters.
    
    ----------
    Tev : TYPE
      in K   DESCRIPTION.
    Tcond : TYPE
        DESCRIPTION.
    Tin : TYPE, optional
        DESCRIPTION. The default is Tin.
    fluid_s : TYPE, optional
        DESCRIPTION. The default is fluid_s.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    comp = [1]
    evap = fprop.T_prop_sat(Tev, fluid_s, composition=comp, option=1)[1][1] # evaporator pressure
    comp_in = fprop.tp(Tin, evap, fluid_s, composition=comp, option=1)  # compressor entrance at 20C
    cond = fprop.T_prop_sat(Tcond, fluid_s, composition = comp, option=1)[1][1] # condenser pressure
    # isentropic state after compressor
    comp_s = fprop.sp(comp_in[4], cond, fluid_s, composition=comp)
    #print(f"comp_s {comp_s}")
    dhs = comp_s[2] - comp_in[2]
    
    to = Tev-273.15
    tc = Tcond-273.15

    coeff_Q, coeff_P, coeff_m = comper.read_file("HESP-2P_20temsuction.xlsx")
    Q, power_el, mdot = comper.cal_range(to, tc, coeff_Q, coeff_P, coeff_m)
    #power_el = br.bitzer_pol(to, tc, br.cP)
    #mdot = br.bitzer_pol(to, tc, br.cm) / 3600
    dh = power_el / mdot
    comp_out = fprop.hp(comp_in[2]+dh, cond, fluid_s, composition=comp)
    # fld.ALLPROP('HP', comp_in["H"] + dh, cond["P"])  # real stae after compressor
    #   d = fld.ALLPROP('TP', Tout, cond["P"])  # real stae after compressor (Bitzer)
    v_ratio = comp_out[3] / comp_in[3]
    p_ratio = comp_out[1] / comp_in[1]
    etas = dhs / dh
    # degree of delivery
    eta_vol = mdot / (V_h * f * 1/comp_in[3])
    
    return np.array([etas, comp_out[0], comp_out[1], evap, v_ratio,
                     p_ratio, power_el, mdot, eta_vol])

def get_data(name_compressor):
    file_dir = os.path.dirname(os.path.realpath('__file__'))
    file_path = os.path.join(file_dir, 'data', 'R290 Verdichter', 'R290 Verdichter', name_compressor)
    print(file_path)
    read_file = pd.read_excel(file_path, header=33, usecols=[3, 4, 5])
    Te_min = float(read_file.iloc[0, 0][:-2])
    Te_max = float(read_file.iloc[0, 2][:-2])
    Tc_min = float(read_file.iloc[1, 0][:-2])
    Tc_max = float(read_file.iloc[1, 2][:-2])
    read_file2 = pd.read_excel(file_path, header=42, usecols=[2, 3, 4])
    d = read_file2.iloc[0, 0]
    h = read_file2.iloc[0, 1]
    f = read_file2.iloc[0, 2]
    return Te_min, Te_max, Tc_min, Tc_max, d, h, f

if __name__ == "__main__":
    name_compressor = "2CESP-3P.xlsx"
    Te_min, Te_max, Tc_min, Tc_max, d, h, f = get_data(name_compressor)
    fluid = "Propane"
    comp = [1.0]
    points = 10
    Te = np.linspace(Te_min, Te_max, points)+273.15
    Tc = np.linspace(Tc_min, Tc_max, points)+273.15
    V_h = 0.25 * np.pi * d ** 2 * h
    T_in = 20 + 273.15
    fi, ax = plt.subplots(2, 2)
    for i, te in enumerate(Te):
        result = []
        for tc in Tc:
            result.append(eta_calc(te, tc, T_in))
        result = np.array(result)
        ax[0, 0].plot(Tc, result[:, 0], label="%3i" % (te))
        ax[0, 1].plot(Tc, result[:, 1], label="%3i" % (te))
        ax[1, 0].plot(Tc, result[:, -2], label="%3i" % (te))
        ax[1, 1].plot(Tc, result[:, -1], label="%3i" % (te))
    for i in ax.flat:
        i.set_xlabel("tc, condensing temperature in °C")
    ax[0,0].set_ylabel("eta_s")
    ax[0,1].set_ylabel("T_out")
    ax[1,0].set_ylabel("m_dot")
    ax[1,1].set_ylabel("eta_vol")
    ax[1, 1].legend()
    plt.show()


overview_picture = "no"
if overview_picture == "yes":
    fi, ax = plt.subplots(3, 3)
    Te = np.linspace(-25, 0 , 36) + 273.15
    Tc = np.linspace(30, 60, 31) + 273.15
    result = np.zeros([31, 36, 9])
    grid_x = np.linspace(0, 36, 37)
    grid_y = np.linspace(0, 31, 32)
    for i, te in enumerate(Te):
        Tin = te + 10
        for j, tc in enumerate(Tc):
            result[j, i] = eta_calc(te, tc, Tin)
    t = np.linspace(0, 8, 9)
    name = ["isentropic efficiency", "outlet temperature", "outlet pressure", "evaporator pressure", "V_out / V_in", "pressure ratio",
            "electric power", "massflow", "volumetric efficiency"]
    for step in t:
        if any(step==[0,1,2]):
            x = 0
            y = int(step)
        if any(step==[3,4,5]):
            x = 1
            y = int(step-3)
        if any(step==[6,7,8]):
            x = 2
            y = int(step - 6)
        plot1 = ax[x, y].pcolormesh(result[:, :, int(step)])

        ax[x, y].set_xlabel("te, evaporating temperature")
        ax[x, y].set_ylabel("tc, condensing temperature")
        ax[x, y].set_title(name[int(step)])
        plt.colorbar(plot1)

    plt.subplots_adjust(left=0.125,
                       bottom=0.057,
                       right=0.9,
                       top=0.964,
                       wspace=0.2,
                       hspace=0.36)


    #ax.title("eta")
    plt.show()


    #T,p,h,v,s,x







backup = "no"
if backup == "yes":
    fi, ax = plt.subplots(2, 2)
    n_no = 5
    col = ["b.", "ro-", "k", "k.", "bv-"]
    Te = np.linspace(-25, 0, n_no)+273.15
    Tc = np.linspace(30, 60, 10)+273.15
    for i, te in enumerate(Te):
        result = []
        for tc in Tc:
            result.append(eta_calc(te, tc))
        result = np.array(result)
        ax[0, 0].plot(Tc, result[:, 0], col[i], label="%3i" % (te))
        ax[0, 1].plot(Tc, result[:, 1], col[i], label="%3i" % (te))
        ax[1, 0].plot(Tc, result[:, -2], col[i], label="%3i" % (te))
        ax[1, 1].plot(Tc, result[:, -1], col[i], label="%3i" % (te))
    for i in ax.flat:
        i.set_xlabel("tc, condensing temperature in °C")
    ax[0,0].set_ylabel("eta_s")
    ax[0,1].set_ylabel("T_out")
    ax[1,0].set_ylabel("m_dot")
    ax[1,1].set_ylabel("eta_vol")
    ax[1, 1].legend()
    fi.savefig("bitzer2EESP-05PProbe.png")
    #print(eta_calc(Tevap, Tcond))
    plt.show()
