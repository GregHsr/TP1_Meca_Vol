# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 10:25:35 2023

@author: Christophe Airiau

Mécanique du vol,

Module pour le TP1 : atmosphère, pression et vitesse


"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import fsolve
from flight.general import *
from flight.atmosphere import *
from flight.aerodynamics import *


# *** VELOCITIES ***

    
def Pi(q, p):
    """
    Pi function, used in TAS and EAS 
    """
    return pow(1 + q / p, (gamma-1) / gamma) -1


def EAS2TAS(EAS, sigma=1):
    """
    EAS = TAS / sqrt(sigma)
    """
    return EAS / np.sqrt(sigma)


def TAS2EAS(TAS, sigma=1):
    """
    TAS = EAS    sqrt(sigma)
    """
    return TAS * np.sqrt(sigma)
    

def CAS2EAS(q, p, p_0, CAS):
    """
    solve CAS from EAS
    """
    return CAS * np.sqrt(p * Pi(q, p) / (p_0 * Pi(q, p_0))) 


def EAS2CAS(q, p, p_0, EAS):
    """
    solve EAS from CAS
    """
    c = np.sqrt(p_0 * Pi(q, p_0) / (p * Pi(q, p)))
    print("CAS/EAS :  ", c)
    return EAS * c


def CAS2TAS(q, p, p_0, Theta, CAS):
    """
    solve TAS from CAS
    """
    return CAS * np.sqrt(Theta * Pi(q, p) / Pi(q, p_0)) 


def display_velocites(KIAS, KCAS, KEAS, KTAS, KGS):
    """
    in kt and in km/h
    """
    n = 70
    print()
    print("=" * n)
    print("   KIAS" + " "*8 + "KCAS" + " "*8 + "KEAS" + " "*8 + "KTAS" + " "*8
          + "KGS     [kt]")
    print("-" * n)
    form = "  %7.2f   " * 5
    print(form % (KIAS, KCAS, KEAS, KTAS, KGS))
    print("=" * n)
    print("    IAS" + " "*9 + "CAS" + " "*9 + "EAS" + " "*9 + "TAS" + " "*9
          + "GS     [km / h]")
    print("-" * n)
    print(form % (kt2ms(KIAS) * 3.6, kt2ms(KCAS) * 3.6, kt2ms(KEAS) * 3.6, 
                  kt2ms(KTAS) * 3.6, kt2ms(KGS) * 3.6))
    print("=" * n)
    
    
def display_values(h, KIAS, KCAS, KEAS, KTAS, KGS,
                   p, T, rho, a, q, Mach, qc):
    """ 
    display velocities, atmosphere quantities and dynamic pressures
    """
    display_state("Altitude  : %.2f m" % h, p, rho, T, a)
    display_velocites(KIAS, KCAS, KEAS, KTAS, KGS)
    display_dynamic_pressure(Mach, q, qc)
    

def TAS_from_CAS(h_ft=25000, KCAS=303, atm0={'p': 101325, 'T': 288.15}):
    """ 
    h : altitude en feets
    """
    CAS = kt2ms(KCAS)
    atm0["rho"] = atm0["p"] / (r_ref * atm0["T"])
    atm0["a"] = sound_velocity(atm0["T"])
    Theta, delta, sigma = atmos(ft2m(h_ft))
    # rho = sigma * atm0["rho"]
    T = Theta * atm0["T"]
    p = delta * atm0["p"]
    a = sound_velocity(T)
   
    qc = dynamic_pressure(p=atm0["p"], M=CAS/atm0["a"], pi=True)   
    EAS = CAS2EAS(qc, p, atm0["p"], CAS)
    KEAS = ms2kt(EAS)
    TAS = EAS / np.sqrt(sigma)
    # TAS = CAS2TAS(qc, p, atm0['p'], Theta, CAS)
    KTAS = ms2kt(TAS)
 
    Mach = TAS / a
    q = dynamic_pressure(p=p, M=Mach)
    return KEAS-KCAS, KTAS - KCAS, Mach, q

def velocities_from_KEAS(h_ft=25000, KEAS=292, Kwind=20, Delta_error=3,
                        atm0={'p': 101325, 'T': 288.15}):
    """ 
    h_ft : altitude en feets
    """
    print("*" * 50)
    print("Velocities from KEAS")
    print("*" * 50)
    atm0["rho"] = atm0["p"] / (r_ref * atm0["T"])
    atm0["a"] = sound_velocity(atm0["T"])
    Theta, delta, sigma = atmos(ft2m(h_ft))
    rho = sigma * atm0["rho"]
    T = Theta * atm0["T"]
    p = delta * atm0["p"]
    a = sound_velocity(T)
    #         
    EAS = kt2ms(KEAS)
    #
    TAS = EAS2TAS(EAS, sigma)
    KTAS = ms2kt(TAS)
    KGS = KTAS - Kwind
    M = TAS / a                                       
    q = q_from_Mach(M, p)          # c'est la vrai pression dynamique ici
    qc = qc_from_Mach(M, p) 
    CAS = EAS2CAS(qc, p, atm0["p"], EAS) 
    KCAS = ms2kt(CAS)
    KIAS = KCAS - Delta_error
    display_values(ft2m(h_ft), KIAS, KCAS, KEAS, KTAS, KGS,
                   p, T, rho, a, q, M, qc)
    
def velocities_from_KTAS(h_ft=25000, KTAS=436, Kwind=20, Delta_error=3,
                        atm0={'p': 101325, 'T': 288.15}):
    """ 
    h_ft : altitude en feets
    """
    print("*" * 50)
    print("Velocities from KTAS")
    print("*" * 50)
    atm0["rho"] = atm0["p"] / (r_ref * atm0["T"])
   
    atm0["a"] = sound_velocity(atm0["T"])
    Theta, delta, sigma = atmos(ft2m(h_ft))
    rho = sigma * atm0["rho"]
    T = Theta * atm0["T"]
    p = delta * atm0["p"]
    a = sound_velocity(T)
    #
    TAS = kt2ms(KTAS)         
    KEAS = TAS2EAS(KTAS, sigma)
    KGS = KTAS - Kwind
    M = TAS / a                                       
    q = q_from_Mach(M, p)          # c'est la vrai pression dynamique ici
    qc = qc_from_Mach(M, p) 
    KCAS = EAS2CAS(qc, p, atm0["p"], KEAS) 
    KIAS = KCAS - Delta_error
    display_values(ft2m(h_ft), KIAS, KCAS, KEAS, KTAS, KGS,
                   p, T, rho, a, q, M, qc)
    
    
def velocities_from_KIAS(h_ft=25000, KIAS=300, Kwind=20, Delta_error=3,
                        atm0={'p': 101325, 'T': 288.15}):
    """ 
    h_ft : altitude en feets
    """
    print("*" * 50)
    print("Velocities from KIAS")
    print("*" * 50)
    atm0["rho"] = atm0["p"] / (r_ref * atm0["T"])
   
    atm0["a"] = sound_velocity(atm0["T"])
    Theta, delta, sigma = atmos(ft2m(h_ft))
    rho = sigma * atm0["rho"]
    T = Theta * atm0["T"]
    p = delta * atm0["p"]
    a = sound_velocity(T)
    #
    KCAS = KIAS + Delta_error
    CAS = kt2ms(KCAS)
    qc = qc_from_Mach(CAS/atm0["a"], atm0["p"])
    KEAS = CAS2EAS(qc, p, atm0["p"], KCAS)
    KTAS = EAS2TAS(KEAS, sigma)
    KGS = KTAS - Kwind
    M = kt2ms(KTAS) / a                                       
    q = q_from_Mach(M, p)             
    display_values(ft2m(h_ft), KIAS, KCAS, KEAS, KTAS, KGS,
                   p, T, rho, a, q, M, qc)

     
def solve(h_ft, KCAS_table=np.linspace(100, 650, 501)):
    """
    h_ft        : table of altitude
    KCAS_table  : table of KCAS 
    """
    data = []
    for h in h_ft:
        KCAS, dKEAS, dKTAS, M, qt = [], [], [], [], [] 
        for kcas in KCAS_table:
            Delta_V, Delta_Vt, Mach, qd = TAS_from_CAS(h_ft=h, KCAS=kcas)
            if Mach <= 1:
                KCAS.append(kcas)
                dKEAS.append(Delta_V)
                dKTAS.append(Delta_Vt)
                M.append(Mach)
                qt.append(qd)
            else:
                break
        data.append([h/1000, KCAS, dKEAS, dKTAS, M, qt])
    return data


def plot_diagram(data, kc=2):
    """
    component is given by cas
    kc = 2 : KEAS - KCAS
       = 3 : KTAS - KCAS
       = 4 : Mach
       = 5 : q 
    """
    title = ("", "KCAS", r"$\Delta= KEAS - KCAS$", r"$\Delta=KTAS - KCAS$", "Mach", "q")
    sub = "    and h / 1000 in ft"
    locator =([], [50, 10], [2, 1], [50, 10],[0.1, 0.02], [10000, 2000])
    fig = plt.figure(figsize=(9, 7))
    plt.title(title[kc] + sub)
    plt.xlabel("KCAS")
    plt.ylabel(title[kc])
    for k in range(len(data)):
        cas = data[k][1]
        print(data[k][0], len(cas))
        plt.plot(cas, data[k][kc], lw=2)
        if kc == 1:
            plt.text(cas[-1]-1, data[k][kc][-1]-1, "%d" % data[k][0], rotation=-60)
        elif kc == 2:
            plt.text(cas[-1]+10, data[k][kc][-1], "%d" % data[k][0])
        else:
            plt.text(cas[-1]+10, data[k][kc][-1], "%d" % data[k][0])
    axs = fig.get_axes()
    for ax in axs:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(locator[1][0]))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(locator[1][1]))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(locator[kc][0]))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(locator[kc][1]))
    plt.grid(which="major", axis="both", color="black", alpha=0.5)
    plt.grid(which="minor", axis="both", color="grey", alpha=0.5)


def get_diagrams_CAS_2_TAS():
    """ 
    """
    h_max = m2ft(11000)
    print("hmax : %.1f ft" % h_max)

    KCAS_table = np.linspace(100, 650, 501)
    h_ft = np.linspace(5000, 35000, 7)
    data = solve(h_ft, KCAS_table)
    for k in range(2, 6):
        plot_diagram(data, k)
    
    