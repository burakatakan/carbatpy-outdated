# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 15:46:10 2022

@author: atakan
"""

import ctypes as ct
import os



class ReVal(ct.Structure):
    _fields_ = [("T", ct.c_double),
                ("D", ct.c_double),
                ("P", ct.c_double),
                ("U", ct.c_double),
                ("H", ct.c_double),
                ("S", ct.c_double),
                ("G", ct.c_double),
                ("A", ct.c_double),
                ("CP", ct.c_double),
                ("CV", ct.c_double),
                ("WS", ct.c_double),
                ("B", ct.c_double),
                ("C", ct.c_double),
                ("CP0", ct.c_double),
                ("Q", ct.c_double)
                ]

if __name__ == "__main__":
    dp = 'C:\\Users\\atakan\\sciebo\\Progs\\TREND_50\TREND_x64.dll'
    # os.chdir(dp)
    t_lib = ct.windll.LoadLibrary(dp)
    print(ct.cdll.TREND_x64)

    t_lib.ALLPROP_STDCALL.restype = ReVal
