# -*- coding: utf-8 -*-
"""
Compressor model of DEnnis Roskosch
changed to work with the actual Refprop module 10.12 2022 (B. Atakan)
It uses extensivly global variables. This should be changed and rewritten 
to use numba and an ODE-solver

Created on Thu Jan 31 17:46:25 2019

@author: roskosch
"""
import numpy as np
import matplotlib.pyplot as plt
from fl_props_compressor import z_uv, z_ps, z_Tp, z_Tx, z_mm, z_px



Rm = 8.3145  # Gaskonstante J/mol/K
Tu = 25.+273.15  # Umgebungstemperatur
dTK = 273.15  # Umrechnung °C / K
Ver0 = [34e-3, 34e-3, 2., .04]  # Fit-Verdichter: D, H, Zylinder, Außenfläche


###########################################
#  Verdichterspezifische Parameter
# pV:    0 - D, Kolbendurchmesser, m
#        1 - H, Hub, m
#        2   Verhältnis aus Kurbelwellenradius und Pleuellänge
#        3   A_a, Oberfläche Verdichter, m²
#        4 - c1, Totvolumen relativ zum Hubvolumen
#        5 - preib, Reibdruck, kPa
#        6 - f, elektrische Verdichterfrequenz Hz
#        7 - Verdichterdrehzahl Hz
#        8 - Zylinderanzahl 



# pZ: 
# Index S: Saugseite
# Index D: Druckseite
#       0 - TS, °C  # BA: seems to be K
#       1 - pS, kPa
#       2 - vS, m³/kg
#       3 - uS, kJ/kg
#       4 - hS, kJ/kg
#       5 - sS, J/kg/K
#       6 - pD, kPa

#pZyk: Größen, die während eines Zyklusses konstant bleiben bzw. über die
#      Zyklen iteriert werden
#       0 - Aeff_e, effektiver Ströumgsquerschnitt Eintritt, m²
#       1 - Aeff_a, effektiver Strömungsquerschnitt Austritt, m²


# z_it: Nimmt für jeden Iterationsschritt folgende Werte auf
#       0 - tet, Kurbelwinkel
#       1 - x, Kolbenposition
#       2 - V, Zylindervolumen
#       3 - Fläche Wärmeübertragung Zylinderwand 
#       4 - Zyklusschritt, 0=Verdichten, 1=Ausschieben, 2=Expandieren, 3=Ansaugen
#       5 - T, Temperatur im Zylinder °C  # BA: seems to be K
#       6 - p, Druck im Zylinder kPa
#       7 - v, spezifisches Volumen im Zylinder m³/kg
#       8 - u, innere Energie im Zylinder kJ/kg
#       9 - h, Enthalpie im Zylinder kJ/kg
#      10 - s, Entropie im Zylinder J/kg/K
#      11 - m, Masse im Zylinder kg
#      12 - T_th, Temperatur der thermischen Masse
#      13 - alpha, Wärmeübergangskoeffizient innen
#      14 - dm, Eingeströmte bzw. ausgeströmte Masse im Iterationsschritt kg
#      15 - Übertragene Wärme






IS = 360  # Anzahl der differentiellen Schritte für einen Zyklus
IS0 = IS
pV = np.zeros(8, float)
pZ = np.zeros(7, float)
pZyk = np.zeros(2, float)
z_it = np.zeros([IS, 16])
fluid = []
comp = [1.0]  # must be checked BA


def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res

def getalp():
    global z_it
    '''
    Berechnet den Wärmeübergangskoeffizient Gas/Zylinderwand
    Woschni Korrelation
    '''
    if z_it[i,  4] == 0 or z_it[i, 4] == 2:  # Ventile geschlossen
        k = 2.28
    else:  # Ventile offen, Ansaugen oder Ausstoßen
        k = 5.18
    v_p = np.abs(z_it[i, 1] - z_it[i-1, 1]) / ((z_it[i, 0] - z_it[i-1, 0]) /
            (2. * np.pi * pV[7]))  # dX/dt  
    alp =127.93 * pV[0]**(-.2) * (z_it[i-1, 6]*1e-2)**.8 * \
            (z_it[i-1, 5] )**(-.55) * (k * v_p)**.8
    z_it[i, 13] = alp
    return alp


def Zustand_th_Masse(Q):
    global z_it
    '''Berechnet die Temperaturänderung der thermischen Masse als Funktion des
    Wärmeübergang Innen (Q) und des Wärmeübergangs zur Umgebung (Q_u)
    '''
    ### Masse und cv der thermischen Masse sind im stationären Betrieb nicht 
    ### entscheidend, die Parameter sind so gewählt, dass das Modell schnell konvergiert
    ### aber keine Schwingungen auftreten
    m = .0001  # kg
    cv = .502  # kJ/kg/K
    alp_a = 6.  # Wärmeübergangskoeffizient zur Umgebung
    A = Ver0[3] * pV[8] / Ver0[2] * pV[0] / Ver0[0] * pV[1] / Ver0[1]  # Außenfläche Zylinder abgeschätzt über Geometrie bezogen auf Fitting-Verdichter
    Q_u = alp_a * A * (Tu -z_it[i-1, 12]) * ((z_it[i, 0] - z_it[i-1, 0]) / 
          (2. * np.pi * pV[7])) * 1e-3  # kJ
    z_it[i, 12] = (Q + Q_u) / cv / m + z_it[i-1, 12]
    return


def Kompression():
    global z_it, fluid, comp
    Schritt = 0
    W = -z_it[i-1, 6] * (z_it[i, 2] - z_it[i-1, 2])  # Kompressionsarbeit, kJ
    Wr = -pV[5] * (z_it[i, 2] - z_it[i-1, 2])  # Reibarbeit, kJ
    getalp()
    Q = z_it[i, 13] * z_it[i, 3] * (z_it[i-1, 12] - z_it[i-1, 5]) * \
        ((z_it[i, 0] - z_it[i-1, 0]) / (2. * np.pi * pV[7])) * 1e-3  # kJ
    Zustand_th_Masse(-Q)
    dm = 0.  # Keine Ein- oder Ausströmende Masse
    mi = z_it[i-1, 11]  # Masse im Zylinder, kg
    ui = (Q + W + Wr) / mi + z_it[i-1, 8]  # kJ/kg
    vi = z_it[i, 2] / mi  # Spezifisches Volumen im Zylinder, m³/kg
    # Abfrage Stoffdaten Zustand i mit u und v, zurück: T, p, v, u, h, s
    zi = z_uv(ui,vi,fluid, comp) # fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)

    # Setzen der Daten in z_it
    z_it[i, 4] = Schritt
    z_it[i, 5:11] = zi
    z_it[i, 11] = mi
    z_it[i, 14:16] = dm, Q
    return z_it


def Ausschieben():
    global z_it, fluid, comp
    Schritt = 1
    W=-z_it[i-1,6]*(z_it[i,2]-z_it[i-1,2])#Kompressionsarbeit, kJ
    Wr=-pV[5]*(z_it[i,2]-z_it[i-1,2])#Reibarbeit, kJ
    getalp()
    Q = z_it[i, 13] * z_it[i, 3] * (z_it[i-1, 12] - z_it[i-1, 5]) * \
        ((z_it[i, 0] - z_it[i-1, 0]) / (2. * np.pi*pV[7]))*1e-3  # kJ
    Zustand_th_Masse(-Q)
    m_dot = pZyk[1] / z_it[i-1, 7] * np.sqrt(2. * (z_it[i-1, 6] - pZ[6]) * \
            1000. * z_it[i-1, 7])  # Massenstrom der den Zylinder verlässt, kg/s
    dm = m_dot * ((z_it[i, 0] - z_it[i-1, 0]) / (2. * np.pi * pV[7]))
    # Ausgeschobene Masse, kg
    mi = z_it[i-1, 11] - dm  # Resultierende Masse im Zylinder, kg
    vi = z_it[i, 2] / mi  # Resultierendes spezifisches Volumen im Zylinder, m³/kg
    ui = (Q + W + Wr - dm * z_it[i-1, 9] + z_it[i-1, 11] * z_it[i-1, 8]) / mi
    # Energiebilanz, resultierende Innere Energie im Zylinder, kj/kg
    # Abfrage Stoffdaten Zustand i mit u und v
    zi = z_uv(ui, vi, fluid, comp)  # fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)

    # Setzen der Daten in z_it
    z_it[i, 4] = Schritt
    z_it[i, 5:11] = zi
    z_it[i, 11] = mi
    z_it[i, 14:16] = dm, Q


def Expansion():
    global z_it, fluid, comp
    Schritt = 2            
    W = -z_it[i-1, 6] * (z_it[i, 2] - z_it[i-1, 2])  # Kompressionsarbeit, kJ
    Wr = pV[5] * (z_it[i, 2] - z_it[i-1, 2])  # Reibarbeit, kJ
    getalp()
    Q = z_it[i, 13] * z_it[i, 3] * (z_it[i-1, 12] - z_it[i-1, 5]) * \
        ((z_it[i, 0] - z_it[i-1, 0]) / (2. * np.pi * pV[7])) * 1e-3  # kJ

    Zustand_th_Masse(-Q)
    dm = 0.  # Keine Ein- oder Ausströmende Masse
    mi = z_it[i-1, 11]  # Masse im Zylinder, kg
    ui = (Q + W + Wr) / z_it[i-1, 11] + z_it[i-1, 8]  # Energiebilanz
    vi = z_it[i, 2] / mi  # Resultierendes spezifisches Volumen im Zylinder, m³/kg
    # Abfrage Stoffdaten Zustand i mit u und v
    zi = z_uv(ui, vi, fluid, comp)  # zi=fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)

    # Setzen der Daten in z_it
    z_it[i, 4] = Schritt
    z_it[i, 5:11] = zi
    z_it[i, 11] = mi
    z_it[i, 14:16] = dm, Q


def Ansaugen():
    global z_it, fluid, comp
    Schritt = 3
    W = -z_it[i-1, 6] * (z_it[i, 2] - z_it[i-1, 2])  # Kompressionsarbeit, kJ
    Wr = pV[5] * (z_it[i, 2] - z_it[i-1, 2])  # Reibarbeit, kJ
    getalp()
    Q = z_it[i, 13] * z_it[i, 3] * (z_it[i-1, 12] - z_it[i-1, 5]) * \
        ((z_it[i, 0] - z_it[i-1, 0]) / (2. * np.pi * pV[7])) * 1e-3  # kJ
    Zustand_th_Masse(-Q)
    
    m_dot = pZyk[0] / pZ[2] * np.sqrt(2. * (pZ[1] - z_it[i-1, 6]) * 1000 * pZ[2])  # Massenstrom der den Zylinder verlässt, kg
    dm = m_dot * ((z_it[i ,0] - z_it[i-1, 0]) / (2. * np.pi * pV[7]))  # Ausgeschobene Masse, kg
    mi = z_it[i-1, 11] + dm  # Resultierende Masse im Zylinder, kg
    vi = z_it[i, 2] / mi  # Resultierendes spezifisches Volumen im Zylinder, m³/kg
    ui = (Q + W + Wr + dm * pZ[4] + z_it[i-1, 11] * z_it[i-1, 8]) / mi  # isenthalpe Drosselung auf Zylinderdruck
    # Abfrage Stoffdaten Zustand i mit u und v
    zi = z_uv(ui, vi, fluid, comp)   # zi=fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)
    # Setzen der Daten in z_it
    z_it[i, 4] = Schritt
    z_it[i,5: 11] = zi
    z_it[i,11] = mi
    z_it[i,14:16] = dm, Q
    
    

def ProzessIter():
    global pZyk, z_it, i, IS, IS0, fluid, comp
    
    #Setzen von Aeff_e, explizite Funktion
    M = z_mm(300, 100.,fluid, comp)[-1]  # CP.PropsSI("M",fluid) # Molmasse kg/mol
    pZyk[0] = 2.0415e-3 * (Rm / M)**(-.9826) * pV[0]**2. / Ver0[0]**2.  # Effektiver Strömungsquerschnitt Eintritt, m²
    
    #Setzen von Aeff_a, implizite Funktion bzgl der mittleren Massenstromdichte über Ventil
    #Bei der 1. Iteration ist Massenstromdichte unbekannt, typischer Wert wird gesetzt
    pZyk[1] = 1.5e-5 * pV[0]**2. / Ver0[0]**2.
    # print(pZyk)
    zaehl=0
    f, ax = plt.subplots(1,1)

    while 1:
        zaehl += 1
        
        if zaehl>1:
            ax.plot(z_it[:, 0], z_it[:, 5], z_it[:, 12])
        for i in range(1, IS):
            if z_it[i,0] <= np.pi:
                if z_it[i-1, 6] <= pZ[6]:
                    Kompression()
                else:
                    Ausschieben()
            else:
                if z_it[i-1,6] >= pZ[1]:
                    Expansion()
                else:
                    Ansaugen()

        # Fehlerquadratsumme T, p, T_th_mittel
        error = np.sqrt((z_it[-1, 5] - z_it[0, 5])**2.) + np.sqrt((z_it[-1, 6]
                   - z_it[0, 6])**2.) + np.sqrt((np.average(z_it[-1, 12])
                   - np.average(z_it[0, 12]))**2.)
        # print(IS,zaehl)
        print("Iteration:", zaehl, error, IS)


        if error < .01: #Zweistufig erst viele Rechnungen mit geringer Aufloesung (BA)
            if IS == IS0:
                IS = 5 * IS0  # dann erhoehte Aufloesung
                z0_ = z_it[:,:]
                z_it = np.zeros([IS, 16])
                z_it[:IS0,:] = z0_
                z_it[IS0:,:] = z0_[-1,:]
                geometrie()  # die Winkel etc. richtig berechnen
            else:
                break
        else:
            Zellen_ausschieben =  find(z_it[:, 4] == 1)
            m_aus = np.sum(z_it[Zellen_ausschieben, 14])  # Insgesamt ausgeschobene Masse
            t_aus = (z_it[Zellen_ausschieben[-1], 0] - 
                   z_it[Zellen_ausschieben[0], 0]) / (2. * np.pi * pV[7])  # Zeit des Ausschiebens
            m_dichte = m_aus / t_aus / pZyk[1]  # Massenstromdichte kg/s/m²
            #print("BA:",m_dichte,pV[0], error)
            pZyk[1] = 5.1109e-4 * (m_dichte)**(-.486) * pV[0]**2. / Ver0[0]**2.  # Aeff_a neu
            z_it[0,5:14] = z_it[-1, 5:14]  # Endwerte letzter Zyklus = Anfangswerte nächster Zyklus
    
    
    # Auswertung Wirkungsgrade    
    Zellen_ausschieben = find(z_it[:, 4] == 1)  # BA find()
    m_aus = np.sum(z_it[Zellen_ausschieben, 14])  # Insgesamt ausgeschobene Masse
    m0 = np.pi * pV[0]**2. * pV[1] / pZ[2] / 4.  # Angesaugte Masse idealer Verdichter
    Liefer = m_aus / m0 #Liefergrad
    
    h_aus = np.sum(z_it[Zellen_ausschieben, 9] * z_it[Zellen_ausschieben, 14])\
        / m_aus  # Mittlere Enthalpie Ausschieben
    h_aus_s = z_ps(pZ[6], pZ[5], fluid, comp)[4]  # fl.zs_kg(['p','s'],[pZ[6],pZ[5]],['h'],fluid)[0]  # Austrittsenthalpie isentrop
    Isentr = (h_aus_s - pZ[4]) / (h_aus - pZ[4])  # Isentroper Wirkungsgrad

    return Isentr, Liefer


def getETA(T_e, p_e, p_a, fluid_in, comp):
    global pV, pz, z_it
    fluid = fluid_in
    comp = comp

###############################     Parametersatz spezifisch für Verdichter   ###################################
  
    pV = [34e-3, 34e-3, 3.5, .04, .06071, 48.916, 50., 50. / 2., 2.]  # Parameter siehe oben

#################################################################################################################
    pZ[0:6] = z_Tp(T_e,p_e,fluid, comp)  # fl.zs_kg(['T','p'],[T_e,p_e],['T','p','v','u','h','s'],fluid) #Zusatnd Saugleitung
    pZ[6] = p_a #Druck in Druckleitung
    print(pZ)
    ############### Setze Geometrie ##################################
    z_it[:, 0] = np.linspace(0., 2 * np.pi, IS)
    z_it[:, 1] = -(pV[1] / 2. * (1. -np.cos(z_it[:, 0]) + pV[2] * 
        (1.-np.sqrt(1. - (1. / pV[2] * np.sin(z_it[:,0]))**2.)))) + \
        pV[4] * pV[1] + pV[1]  # Kolbenpositionen, x=0 bei UT
    A_kopf = np.pi / 4. * pV[0]**2.
    z_it[:,2] = A_kopf * z_it[:,1]  # Zylindervolumina
    z_it[:, 3] = np.pi * pV[0] * z_it[:,1] + 2. * A_kopf  # Wärmeübertragungsflächen
    
    ########## Setze Startbedingungen im Zylinder ##################
    z_it[0, 5:11] = pZ[0:6]    
    z_it[0, 11] = z_it[0,2] / z_it[0, 7]  # V/v, Zylinder vollständig gefüllt mit Sauggas
    ##### Starttemp thermische Masse, pro Zyklusdurchlauf gibt es stehts nur eine Temperatur
    ##### Startwert frei wählbar, beeinflusst maßgeblich Iterationszeit.
    z_it[:, 12] = 42.+273
    Isentr,Liefer = ProzessIter()
    return np.array((Isentr, Liefer))

def geometrie():
    global pV, pz, z_it, IS
    z_it[:, 0] = np.linspace(0., 2 * np.pi, IS)
    z_it[:, 1] = -(pV[1] / 2. * (1. -np.cos(z_it[:, 0]) + pV[2] * 
        (1.-np.sqrt(1. - (1. / pV[2] * np.sin(z_it[:,0]))**2.)))) + \
        pV[4] * pV[1] + pV[1]  # Kolbenpositionen, x=0 bei UT
    A_kopf = np.pi / 4. * pV[0]**2.
    z_it[:,2] = A_kopf * z_it[:,1]  # Zylindervolumina
    z_it[:, 3] = np.pi * pV[0] * z_it[:,1] + 2. * A_kopf  # Wärmeübertragungsflächen




#Beispiel #################################################
fluid = 'Propane * Hexane *Butane'
comp = [1.0, 0.]
fluid = "Isobutane *Propane"  #"Dimethylether"  # "Butane"
x0=.25
x1 = 1-x0 
comp = [x0, x1, 1-x0-x1]
pe0 = z_Tx(273.15+17.2, 0, fluid, comp)[1]  # fl.zs_kg(['T','q'],[0.,0.],['p'],fluid)[0]
pa0 =  z_Tx(53+273.15, 0, fluid, comp)[1] # fl.zs_kg(['T','q'],[35.,0.],['p'],fluid)[0]
pe = pe0 #190.
pa = 1885.  # 8 * pe

T0 =z_px(pe,1, fluid, comp)[0]+20
T0,pe, pa=(278., 178., 1050.0)
fo = open("Daten.txtx","w")
print("Drücke %2.2f kPa %2.2f kPa" % (pe, pa))
dt_all=np.array([0.15])
out=[]
for dt in dt_all:
    o1 = getETA(dt + T0,pe,pa,fluid, comp)
    #o1.append((np.max(z_it[:,11]) - np.min(z_it[:,11]) * pV[7]))  # Massenstrom
    out.append(o1)
    print(dt, o1)
out = np.array(out)
plt.plot(dt_all, out)
fo.write(str(out))
fo.close()







