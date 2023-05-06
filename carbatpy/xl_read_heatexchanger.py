# -*- coding: utf-8 -*-
"""
Einlesen einer Exel-Datei mit pandas
-Demo-
Created on Wed Jan 24 12:27:07 2018

@author: atakan
"""

# from pylab import *
import pandas as pd
import numpy as np

# directory = r"C:\Users\atakan\sciebo\Lehre\ThermoUni\ThermoD_19_20"
fn = "heat_exchanger_input0.xlsx"          # Dateiname


def read_hex_file(fn, out):
    # global d, dfs, v_sp, compo_all
    druck = False
    xl = pd.ExcelFile(fn, engine='openpyxl')   # Einlesen
    if druck:
        print("Blattnamen in Datei: %12s" % (xl.sheet_names))

    data = []
    fluids = ["", ""]
    compo_all = []
    dfs = []

    for sheet in xl.sheet_names:
        blatt0 = sheet    # erstes Blatt in der Datei
        df1 = xl.parse(blatt0)
        dfs.append(df1)
        d = df1.values
        data.append(d)
        # print(sheet, sheet == 'Geometry', d[:,0])
        # print("Alles im Blatt: \n%10s" % (df1))
        n_app = ""
        v_sp = np.where(df1.columns == "value")

        if sheet[:3] == 'Flu':
            fln = int(sheet[-1])-1
            n_app = str(fln)
            fl_names = []
            compo = []
            for j, i in enumerate(df1["name_fluid"]):
                if pd.notna(i):
                    fl_names.append(i)
                    n_zeile = np.where(d[:, 0] == 'fl'+str(j-1))[0]
                    compo.append(*d[n_zeile, v_sp])
                    if druck:
                        print(j, i, compo, n_zeile, v_sp, 'fl'+str(j-1))

            fluid_no = len(fl_names)
            if fluid_no > 1:
                fluids[fln] = "*".join(fl_names)

            else:
                fluids[fln] = fl_names[0]
                compo = [1.]
            compo_all.append(np.hstack(compo))

        for i, name in enumerate(d[:, 0]):
            if pd.isna(d[i, 2]):
                globals()[name+n_app] = d[i, 1]
            else:
                globals()[name+n_app] = d[i, 1:3]

    xl.close()
    if out == "HEXsimple":
        # all what is needed for heat_exchanger.counterflow_hex :
        return (fluids, [m_dot0[0], m_dot1[0]], [p_in0[0], p_in1[0]], length,
                d_in, U, tubes, [props0, props1], compo_all)
    else: print(f"{out} is not implemented yet!")


if __name__ == "__main__":
    print(read_hex_file(fn, "HEXsimple"))

    # print("Ein Wert: %10s"%(df1.wert1[2]))
    # a=df1.Wert2[0]
    # # Neue Datei mit dem Dataframe als Inhalt schreiben
    # df1.to_excel("output.xlsx")

    # # Speichern eines Arrays (auf die Schnelle)
    # x = np.zeros((2,5))
    # x0=np.linspace(0,1,5)
    # x[0,:]=x0
    # x[1,:]=2*x0**2 -0.5
    # fb = open("testOut.dat", "bw")
    # x.tofile(fb)
    # fb.close()
    # #Auslesen des Arrays
    # fz=open("testOut.dat", "r")
    # z=np.fromfile(fz)
    # print(z)
    # fz.close()
