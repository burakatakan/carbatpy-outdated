# -*- coding: utf-8 -*-
"""
Part of carbatpy
Reading an excel file with sheets for the definition of an heat exchanger 
(Geometry)
the two fluids (Fluid1 is the working fluid, when used in a cycle and
                Fluid 2 is the secondary fluid) and a sheet with the 
Problem_description, it is not used yet, but will be used for optimizations etc.


@author: atakan
"""

# from pylab import *
import pandas as pd
import numpy as np

# directory = r"C:\Users\atakan\sciebo\Lehre\ThermoUni\ThermoD_19_20"
fn = "heat_exchanger_input1.xlsx"          # Dateiname


def read_hex_file(fn, out, all_out = False):
    """
    Reads an Excel File with information about heat exchangers, fluids etc.
    is used within carbatpy and the module heat_exchanger (see there for the
     specific output)

    Parameters
    ----------
    fn : string
        name of an excel-file with at least 4 sheets as mentioned above, which 
        list the variable names, values and further parameters and comments.
    out : string 
        only "HEXSimple" implemented so far, for the simple heat exchanger 
        calculations with actual fluid properties in a counterflow 
        configuration.
    all_out : boolean
        if True a list with 4 Instances with variable names and values are 
        returned together with all local variables (dictionaries) for usage
        in other modules.

    Returns
    -------
    list 
        if HEXSimple: a list with the required values for setting up the problem
        if (internally).
        if all_out is True, see above.

    """
    
   
    druck = False  # printing?
    xl = pd.ExcelFile(fn, engine='openpyxl')   # Einlesen
    if druck:
        print(f"sheet names in {fn}: {xl.sheet_names}")

    data = []
    dfs = []
    outputs = []

    for sheet in xl.sheet_names:
        blatt0 = sheet    # erstes Blatt in der Datei
        df1 = xl.parse(blatt0)
        dfs.append(df1)
        # d = df1.to_numpy()
        # data.append(d)
        n_app = ""
        c_names = df1.columns
        dict_list = locals()[blatt0] = [v.dropna().to_dict()
                                        for k, v in df1.iterrows()]
        bi = locals()[blatt0+"_Instance"] = type(blatt0, (object,), {})
        if blatt0 == "Geometry":
            Geometry_Instance = bi
        elif blatt0 == "Fluid_1":
            Fluid1_Instance = bi
        elif blatt0 == "Fluid_2":
            Fluid2_Instance = bi
        elif blatt0 == "Problem_description":
            Problem_description = bi
        else:
            print(f"Warning: sheet {blatt0} is not used!")

        fluids = []
        compo_all = []

        for di in dict_list:
            # name a class using a string
            # add an attribute name from a strng
            setattr(bi, di["variable_name"], di)

        if sheet[:3] == 'Flu':

            fl_names = []
            compo = []
            fluid_no = bi.number_compounds["value"]
            for ii in range(fluid_no):
                fname = getattr(bi, "fl"+str(ii+1))
                fl_names.append(fname['name_fluid'])

                compo.append(fname["value"])
                if druck:
                    print(ii, compo, 'fl'+str(ii+1), fl_names)

            if fluid_no > 1:
                print(fl_names, "---")
                fluids = "*".join(fl_names)

            else:
                fluids = fl_names[0]

            compo_all.append(np.hstack(compo))
            setattr(bi, "composition_all", compo_all)
            setattr(bi, "fluidNamesREFPROP", fluids)
            dict_list.append({"composition_all": compo_all,
                             "fluidNamesREFPROP": fluids})

        outputs.append(bi)

    xl.close()
    if all_out:
        return outputs, locals()

    elif out == "HEXsimple":
        outputHex = [[Fluid1_Instance.fluidNamesREFPROP,
                      Fluid2_Instance.fluidNamesREFPROP],
                     [Fluid1_Instance.m_dot["value"],
                     Fluid2_Instance.m_dot["value"]],
                     [Fluid1_Instance.p_in["value"],
                     Fluid2_Instance.p_in["value"]],
                     Geometry_Instance.length["value_1"],
                     [Geometry_Instance.d_in["value_1"],
                     Geometry_Instance.d_in["value_2"]],
                     Geometry_Instance.U["value_1"],
                     Geometry_Instance.tubes["value_1"],
                     [Fluid1_Instance.props["value"],
                     Fluid2_Instance.props["value"]],
                     [Fluid1_Instance.composition_all,
                     Fluid2_Instance.composition_all],

                     ]
            # all what is needed for heat_exchanger.counterflow_hex :
        return outputHex
            

    else:
            print(f"{out} is not implemented yet!")


if __name__ == "__main__":
    a = read_hex_file(fn, "HEXsimple")
    print(a)
    a = read_hex_file(fn, "",True)
    locals().update(a[1])  # set all variables from excel file

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
