# -*- coding: utf-8 -*-
"""
 
@author: Christophe Airiau

Mécanique du vol,

Module pour le TP1 : atmosphère, pression et vitesse

Module about atmosphere

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from flight.general import R_ref, g_ref, T_ref, Delt, Sigm, Thet
 


def define_constantes():
    """
    Physical constantes
    """
    M_air = 28.9645e-3       # Masse molaire de l'air en kg/mol
    r_ref = R_ref / M_air
    k_T = -6.5/1000          # gradient de température par mètre d'altitude (z <= 11 000)
    p_ref = 101325
    print("r de l'air    : %.4f" % r_ref)
    print("p_ref /r      : %.4f" % (p_ref / r_ref))
    return r_ref, k_T, p_ref


def set_exponents(r_ref, k_T):
    """
    exponents in the atmosphere model
    """
    n = -g_ref / (r_ref * k_T)
    m = g_ref / r_ref * 1000
    print("n      : %.5f" % n)
    print("m      : %.5f" % m, "\t Gudmundsson: ",  34.163195)
    return n, m


def h_from_temperature(ratio, k_T=-6.5e-3, t0=15):
    """
    return the altitude in meters
    """
    kappa = k_T / (T_ref + t0)
    return (ratio - 1) / kappa


def atmos(h, k_T=-6.5e-3, t0=15, n=5.255894):
    """
    standard atmosphere ISA 1976, ratios
    h in meters
    """
    kappa = k_T / (T_ref + t0)
    Theta = 1 + kappa * h
    delta = Theta ** n
    sigma = Theta ** (n-1)
    return Theta, delta, sigma


def set_sea_level_state(r_ref, t_0=15, p_ref=101325):
    """
    reference state at sea level
    """
    p_0, T_0 = p_ref, T_ref + t_0
    rho_0 = p_0 / (r_ref * T_0)
    print("rho_0    %.3f kg/m^3" % rho_0)
    return p_0, T_0, rho_0 


def humidity(x):
    """
    rho/rho_0 with humidity, correction factor
    """
    return (1+x) / (1 +1.609*x)


def density2humidity(r):
    """
    humidity x = f(r),  r = rho/rho_0
    """
    return (r - 1) / (1 - r * 1.609)


def plot_density_humidity():
    """
    plot: rho/rho_0 with humidity, correction factor
    """
    x = np.linspace(0, 1)
    plt.figure()
    plt.plot(x, humidity(x), lw=3)
    plt.xlabel("taux d'humidité")
    plt.ylabel(r"$\dfrac{\rho}{\rho_{std}}$")
    plt.grid()
    plt.title("humidity correction factor")


def solve_humidity(r_target=1):
    """
    get humidity from a given correction
    """
    r_min = humidity(1)
    
    if r_target < r_min:
        raise ValueError("Value under the limit of rho / rho_sl = %.3f"
                         % r_min)
    else:
        return density2humidity(r_target)


def display_atmosphere(h, delta, sigma, theta):
    """ 
    p/p_0, sigma / sigma_0, theta / theta_0 
    
    h in feets
    """
    n = 50
    print()
    print("=" * n)
    print("Standard atmosphere")
    title = "  h [ft]    p / p_0    rho / rho_0    T / T_0 "
    subtitle = " " * 15 + Delt + " " * 12 + Sigm + " " * 12 + Thet 
    print("=" * n)
    print(title)
    print(subtitle)
    print("-" * n)
    form = " %6d      %6.4f     %6.4f         %6.4f "
    print(form % (h, delta, sigma, theta))
    print("=" * n)

# MAIN

r_ref, k_T, p_ref = define_constantes()