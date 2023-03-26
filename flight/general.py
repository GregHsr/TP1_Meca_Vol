# -*- coding: utf-8 -*-
"""
General module containing:
    
    constantes
    conversion

Mécanique du vol,

Module pour le TP1 : atmosphère, pression et vitesse

@author: Christophe Airiau 
""" 
 

T_ref = 273.15    # 0 degrés Celcius en Kelvins
g_ref = 9.80665   # Accélération de la gravité moyenne en m/s^2
R_ref = 8.31432   # constante des gaz parfaits
gamma = 1.4

Gam = u"\u03b3"
Rho = u"\u03c1"
Delt = u"\u03b4"
Sigm = u"\u03c3"
Thet = u"\u0398"
 

# ***  UNIT CONVERSION  ***


def kt2ms(v):
    """ knots to m/s"""
    return v * 1.852 / 3.6


def ms2kt(v):
    """ m/s to knots"""
    return v * 3.6 / 1.852


def ft2m(z):
    """ feets to meters"""
    return z * 0.3048


def m2ft(z):
    """ meters to feets"""
    return z / 0.3048