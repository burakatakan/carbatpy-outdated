from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary


os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])
MASS_BASE_SI = RP.GETENUMdll(0, "MASS BASE SI").iEnum
SI = RP.GETENUMdll(0, "SI").iEnum

def fun_rp(p_wf_0, p_sf_0, di, ds, dh, da, area_wf, area_sf, lam_rohr, m_wf, m_sf, z, T_sf_0):

    # Stoffeigenschaften der Fluide: Dichte, Wärmekapazität, kinematische Viskosität, Thermische Leitfähigkeit, Prandtl-Zahl, Temperatur
    Zustand_wf_x = RP.REFPROPdll('Isobutane', "PH", "D;CP;VIS;TCX;PRANDTL;T", MASS_BASE_SI, 0, 0, 10e5, h_wf, [0]).Output[0:6]