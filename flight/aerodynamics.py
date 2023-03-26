# -*- coding: utf-8 -*-
"""
 
@author: Christophe Airiau

Mécanique du vol,

Module pour le TP1 : atmosphère, pression et vitesse

Module about aerodynamics

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import fsolve
from flight.general import *
from flight.atmosphere import * 
from flight.aerodynamics import *

def display_state(msg='state', p=0, rho=0, T=0, a=0):
    """
    display the gas state
    """
    n = 60
    print()
    print("=" * n)
    print("# \t " + msg)
    print("=" * n)
    title = "  p [Pa]     rho [kg/m^3]   a [m/s]   T [K]     theta [°]"
    print(title)
    print("-" * n)
    form = " %6.2f     %6.4f       %6.2f    %6.2f     %6.2f "
    print(form % (p, rho, a, T, T-273.15))
    print("=" * n)
      
def display_state_plus(Mach=0, T_i=0, p_i=0, u=0):
    """
    display the gas state, in addition
    """
    n = 60
    print()
    print("=" * n)
    title = "  Mach      T_i [K ]        p_i [Pa]    u [km/h]    u [kt]"
    print(title)
    print("-" * n)
    form = " %6.2f     %6.4f       %6.2f    %6.2f      %6.2f"
    print(form % (Mach, T_i, p_i, u * 3.6, ms2kt(u)))
    print("=" * n)



# *** PRESSURES ***

def pressure_exponents(gamma=1.4):
    """ 
    exponent in various equations
    """
    g_1 = "(" + Gam + "-1)"
    print()
    print("=" * 30)
    print("Pressure exponents")
    print("=" * 30)
    r = (1 / (gamma-1), 2 /(gamma - 1), (gamma-1) / gamma, gamma /(gamma-1))
    e = (
        "1 / "+ g_1 + " : %.2f",
        "2 / "+ g_1 + " : %.2f",
        g_1 + " / " + Gam + " : %.5f",
        Gam + " / " + g_1 + " : %.2f")
    for leg, x in zip(e, r):
        print(leg % x)


def dynamic_pressure(rho=0, M=0, U=0, p=0, pi=False, gamma=1.4):
    """ 
    Le calcul de la pression dynamique dépendra des arguments 
    employés
        rho, U
        p, M
    si pi = True on retourne q_c, la pression dynamique compressible
    """
    if rho * U > 0:
        return 1 / 2 * rho * U**2
    elif p * M > 0:
        if pi:
            return p * (pow(1 + (gamma-1)/2 * M**2, gamma/(gamma-1)) - 1)
        else:
            return gamma * p * M**2 / 2
    else:
        raise ValueError("probleme : Il manque une informatio")


def pi_sur_p(M, gamma=1.4):
    """
    p_i / p = f(Mach)
    """
    return pow(1 + (gamma-1) / 2 * M**2, gamma/(gamma-1))


def q_from_Mach(M, p, gamma=1.4):
    """
    dynamic pressure
    q = f(Mach, p)
    """
    return p * gamma / 2 * M**2


def qc_from_Mach(M, p, gamma=1.4):
    """
    q_c = p_i - p = f(Mach, p)
    compressible dynamic pressure
    """
    return  (pi_sur_p(M, gamma=gamma) - 1) * p


def pt_sur_p(M, gamma=1.4):
    """
    p_t / p = 1 + gamma M^2 / 2
    """
    return 1 + gamma/2 * M**2


def plot_errors():
    """
    error on total pressure and dynamic pressure
    """
    M = np.linspace(0, 0.9, 91)
    epsilon = 1 - pt_sur_p(M) / pi_sur_p(M) 
    epsilon_q = qc_from_Mach(M[1:], 1,
                             gamma=1.4) /  q_from_Mach(M[1:], 1, gamma=1.4) - 1

    fig, ax = plt.subplots(1, 2, figsize=(9, 7))
    fig.suptitle('différences entre les pressions')
    ax[0].plot(M, epsilon * 100, 'ko', label=r"$\varepsilon$")
    ax[1].plot(M[1:], epsilon_q * 100, 'r-', label=r"$\varepsilon_q$")
    ax[1].plot(M, M**2/4* 100, 'b--', label=r"$M^2/4$")
    ax[0].set_title(r"$\dfrac{p_t-p_i}{p_i}$")
    ax[1].set_title(r"$\dfrac{q_c-q}{q}$")
    for i in range(2):
        ax[i].legend(loc="best")
        ax[i].grid()
        ax[i].set_xlabel("Mach")
        ax[i].set_ylabel("erreur en %")
    return M, epsilon, epsilon_q


def display_dynamic_pressure(M=0.5, q=15000, qc=17000, p_0=101325):
    """
    table of values
    """
    n = 65
    print()
    print("=" * n)
    print("  Mach" + " "*8 + "q [Pa]" + " "*8 + "q_c [Pa]" + " "*4 + "q / p_0",
          " "*4 +  "q_c / p_0")
    print("-" * n)
    form = " %6.3f      %8.3f     %8.3f     %5.3f       %5.3f"
    print(form % (M, q, qc, q / p_0, qc / p_0))
    print("=" * n)
    

# *** MACH NUMBER ***


def sound_velocity(T, gamma=1.4):
    """ a from temperature """
    return np.sqrt(gamma * r_ref * T)


def Mach_from_q(x, gamma=1.4):
    """ 
    x = p / q 
    """
    return np.sqrt(2 / gamma / x)


def Mach_from_pt(x, gamma=1.4):
    """ 
    x = p_t / q 
    """
    return np.sqrt(2 / gamma * (x - 1))


def Mach_from_pi(x, gamma=1.4):
    """ x = p_i / q """
    return  np.sqrt(2 / (gamma - 1) * (pow(x, (gamma-1)/gamma) - 1))


def Mach_number_calculus():
    """
    Mach given from three different conditions
    """
    n = 50
    print()
    print("=" * n)
    print("   q / q    M from q   M from p_t  M from p_i ")
    print("-" * n)
    form =  "  %6.2f     %6.2f     %6.2f     %6.2f  " 
    data = []
    for rapport in [1.2, 10]:
        data.append([rapport, Mach_from_q(rapport), Mach_from_pt(rapport),
                     Mach_from_pi(rapport)])
    for line in data:
        print(form % (line[0], line[1], line[2], line[3]))
    print("=" * n)
    return data


# *** NORMAL SHOCK WAVES ***


def P2_P1(mach, gamma=1.4):
    """ 
    pressure ratio across the shock wave, oblique or normal shock wave
    
    Args:
        * mach (real) : normal Mach number 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : downstream p / upstream p
    """
    return 2.0 * gamma / (gamma + 1.0) * mach ** 2 - (gamma - 1.0) / (gamma + 1.0)


def rho2_rho1(mach, gamma=1.4):
    """    
    density ratio across the shock wave, oblique or normal shock wave
    
    Args:
        * mach (real) : Mach number 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : downstream rho / upstream rho   
    """
    return 1.0 / (2.0 / ((gamma + 1.0) * mach ** 2) + (gamma - 1.0) / (gamma + 1.0))


def pi2_pi1(mach, gamma=1.4):
    """
    isentropic pressure ratio across the shock wave, oblique or normal shock wave
    
    Args:
        * mach (real) : normal Mach number 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : downstream p_i / upstream p_i
    """
    return pow(P2_P1(mach, gamma=gamma), -1 / (gamma - 1)) * pow(rho2_rho1(mach, gamma=gamma), gamma / (gamma - 1))


def downstream_Mach(mach, gamma=1.4):
    """       
    downstream normal Mach across a shock
    
    Args:
        * mach (real) : normal Mach number 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : downstream normal Mach number
    """
    return np.sqrt((1.0 + 0.5 * (gamma - 1.0) * mach ** 2) / (gamma * mach ** 2 - 0.5 * (gamma - 1.0)))


def upstream_Mach(mach, gamma=1.4):
    """       
    upstream normal Mach across a shock when downstream Mach is known
    
    Args:
        * mach (real) : normal Mach number 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : upstream normal Mach number
    """
    return np.sqrt((1.0 + 0.5 * (gamma - 1.0) * mach ** 2) / (gamma * mach ** 2 - 0.5 * (gamma - 1.0)))


# ***  ISENTROPIC FLOWS ***

def p_pi(Mach, gamma=gamma):
    """
    ratio pressure/isentropic pressure function of Mach number

    Args:
        * M (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`

    Returns:
        real : p/p_i
    """
    return (T_Ti(Mach, gamma=1.4)) ** (gamma / (gamma - 1))


def rho_rhoi(Mach, gamma=1.4):
    """
    ratio rho/isentropic rho function of Mach number
    
    Args:
        * M (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`

    Returns:
        real : rho/rho_i
    """
    return (T_Ti(Mach, gamma=gamma)) ** (1. / (gamma - 1))


def T_Ti(Mach, gamma=1.4):
    """
    ratio Temperature/isentropic Temperature function of Mach number
    
    Args:
        * M (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`

    Returns:
        real : T/T_i
    """
    return 1. / (1 + (gamma - 1) / 2 * Mach ** 2)


# ***  PITOT TUBE ***


def pitot_rayleigh(Mach, rhs=0, gamma=1.4):
    """
    pi_2  / p_1 used in the Pitot tube for supersonic flow 
    """
    r_1, r_2 = - 1 / (gamma - 1), gamma / (gamma -1)
    cste = pow((gamma - 1) / (1 + gamma), r_1) * pow((gamma+1) / 2, r_2)
    return cste * pow(2 * r_2 * Mach**2 -1, r_1) * pow(Mach, 2 * r_2) - rhs


def plot_pitot_rayleigh(gamma=1.4, opt=0):
    """
    plot of the function for supersonic flow
    opt = 1 : Delta p / p_0 else : p_i2 / p_1
    """
    Mach = np.linspace(1, 3, 51)
    yleg = (r"$p_{i_2} / p_1$", r"$\dfrac{\Delta p}{p_0}$" )
    r = pitot_rayleigh(Mach, gamma=gamma)
    if opt == 1 :
        r -= 1
    plt.figure()
    plt.title("Pitot tube, supersonic flow")
    plt.ylabel(yleg[opt])
    plt.xlabel("Upstream Mach")
    plt.plot(Mach, r)
    plt.grid()
    plt.show()


def Delta_p_Pitot(Mach, gamma=1.4):
    """
    Delta_p / p_0 for supersonic Pitot tube
    """
    return pitot_rayleigh(Mach, gamma=gamma) - 1


def M1_from_Delta_p(r=4, gamma=1.4):
    """
    r = delta_p / p_0
    return upstream Mach number of a supersonic Pitot tube
    """
    return fsolve(pitot_rayleigh, 1.1, args=(1 + r, gamma))
 

def display_pitot_tube(Mach, p_ratio, delta_p, CAS):
    """
    data used or resuls of the supersonic Pitot tube
    """
    if Mach >= 1:
        n = 65
        print()
        print("=" * n)
        print("Pitot tube (true for supersonic only)")
        print("=" * n)
        print("  Mach     p_i2 / p_1    Delta p [Pa]  CAS [km/h]   KCAS [kt]")
        print("-" * n)
        form = " %6.3f      %6.3f       %8.1f        %6.1f      %6.3f "
        print(form % (Mach, p_ratio, delta_p, CAS*3.6, ms2kt(CAS)))
        print("=" * n)    
    else:
        print("Not possible to display, M must be greater than 1")