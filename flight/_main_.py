# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 10:25:35 2023

@author: Christophe Airiau

Mécanique du vol,

Module pour le TP1 : atmosphère, pression et vitesse


"""

from flight.flight_physics import *

   
def main(option):
    """
    Main methods for testing
    """
    
    if option == 1:
        Theta, delta, sigma = atmos(8000)
        print("Theta : %.4f, delta : %.4f,  sigma : %.4f" % (Theta, delta, sigma))
    
        n, m = set_exponents(r_ref, k_T)
    
        # reference state
        p_0, T_0, rho_0 = set_sea_level_state(r_ref, t_0=15)
    
        Mach_number_calculus()
    
        display_state("h=8000 m", 1e4, 0.5, 260, 300)
    
        display_atmosphere(8000, delta, sigma, Theta)
    
        pressure_exponents()
        
    elif option == 2:
    
        velocities_from_KEAS(h_ft=25000, KEAS=292, Kwind=20, Delta_error=3,
                                atm0={'p': 101325, 'T': 288.15})
        
        velocities_from_KTAS(h_ft=25000, KTAS=436, Kwind=20, Delta_error=3,
                                atm0={'p': 101325, 'T': 288.15})
        
        velocities_from_KIAS(h_ft=25000, KIAS=300, Kwind=20, Delta_error=3,
                                atm0={'p': 101325, 'T': 288.15})
    elif option == 3:
        
        get_diagrams_CAS_2_TAS()
    
    elif option ==4:
        # Tube de Pitot en supersonique
        p_0, T_0, rho_0 = set_sea_level_state(r_ref, t_0=15)
        a_0 = sound_velocity(T_0)
        # Downstream the shock 
        p_2 = 118960
        q_c = 30150
        T_i2 = T_i1 = 401.7
        
        p_i2 = p_2 + q_c 
        M_2 =  np.sqrt(2 / (gamma - 1) * (pow(p_2/p_i2, (1-gamma)/gamma) - 1))
        T_2 = T_Ti(M_2) * T_i2
        rho_2 = p_2 / (r_ref * T_2)
        a_2 = sound_velocity(T_2)
        v_2 = M_2 * a_2
        
        M_1 = downstream_Mach(M_2)
        T_1 = T_Ti(M_1) * T_i2
        a_1 = sound_velocity(T_1)
        v_1 = M_1 * a_1
        p_i1 = p_i2 / pi2_pi1(M_1)
        p_1 = p_2 / P2_P1(M_1)
        rho_1 = p_1 / (r_ref * T_1)
        h = h_from_temperature(T_1 / (T_ref + 15))
        
        q_1 = q_from_Mach(M_1, p_1) 
        q_2 =  q_from_Mach(M_2, p_2)
        
        print("h       : ", h)
        Theta, delta, sigma = atmos(h)
        
        delta_p = (p_i2 / p_1 -1) * p_0
        CAS_1 = M_1 * a_0
        
        # OUTPUTS: 
        display_atmosphere(m2ft(h), delta, sigma, Theta)
          
        print("\n Downstream the shock wave: \n")
           
        display_state(msg='state 2', p=p_2, rho=rho_2, T=T_2, a=a_2)
        display_state_plus(Mach=M_2, T_i=T_i2, p_i=p_i2, u=v_2)
        display_dynamic_pressure(M=M_2, q=q_2, qc=p_i2 - p_2)
        
        print("\n Upstream the shock wave: \n")
        
        display_state(msg='state 1', p=p_1, rho=rho_1, T=T_1, a=a_1)
        display_state_plus(Mach=M_1, T_i=T_i1, p_i=p_i1, u=v_1)
        display_dynamic_pressure(M=M_1, q=q_1, qc=p_i1 - p_1)
        display_pitot_tube(M_1, p_i2 / p_1, delta_p, CAS_1)
        
    elif option == 5:
        h = 10000
        M_1 = 2
        p_0, T_0, rho_0 = set_sea_level_state(r_ref, t_0=15)
        a_0 = sound_velocity(T_0)
        Theta, delta, sigma = atmos(h)
        T_1 = Theta * T_0
        p_1 = delta * p_0
        rho_1 = sigma * rho_0
        T_i1 = T_1 / T_Ti(M_1)
        p_i1 = p_1 / p_pi(M_1)
        a_1 = sound_velocity(T_1)
        v_1 = M_1 * a_1
        
        T_i2 = T_i1  
        p_2 = p_1 * P2_P1(M_1)
        p_i2 = p_i1 * pi2_pi1(M_1)
        qc2 = p_i2 - p_2
        M_2 = downstream_Mach(M_1)
        T_2 = T_i1 * T_Ti(M_2)
        a_2 = sound_velocity(T_2)
        v_2 = M_2 * a_2     
        rho_2 = p_2 / (r_ref * T_2)
        q_1 = q_from_Mach(M_1, p_1) 
        q_2 =  q_from_Mach(M_2, p_2)
        
        delta_p = (p_i2 / p_1 -1) * p_0
        CAS_1 = M_1 * a_0
        
        # OUTPUTS: 
        display_atmosphere(m2ft(h), delta, sigma, Theta)
   
        print("\n Downstream the shock wave: \n")
       
        display_state(msg='state 2', p=p_2, rho=rho_2, T=T_2, a=a_2)
        display_state_plus(Mach=M_2, T_i=T_i2, p_i=p_i2, u=v_2)
        display_dynamic_pressure(M=M_2, q=q_2, qc=p_i2 - p_2)
        
        print("\n Upstream the shock wave: \n")
        
        display_state(msg='state 1', p=p_1, rho=rho_1, T=T_1, a=a_1)
        display_state_plus(Mach=M_1, T_i=T_i1, p_i=p_i1, u=v_1)
        display_dynamic_pressure(M=M_1, q=q_1, qc=p_i1 - p_1)
        display_pitot_tube(M_1, p_i2 / p_1, delta_p, CAS_1)
        
    elif option == 6:
         
        # Tube de Pitot en supersonique, on a Delta p / p_0
        p_0, T_0, rho_0 = set_sea_level_state(r_ref, t_0=15)
        a_0 = sound_velocity(T_0)
        # Downstream the shock 
        r = 4.64   # Delta_p / p_0
        delta_p = r * p_0
        p_1 = 27000
        T_i2 = T_i1 = 401.7
        
        M_1 = M1_from_Delta_p(r=r, gamma=1.4)
        CAS_1 = M_1 * a_0
        T_1 = T_Ti(M_1) * T_i1
        a_1 = sound_velocity(T_1)
        v_1 = M_1 * a_1
        p_i1 = p_1 / p_pi(M_1)
        rho_1 = p_1 / (r_ref * T_1)
        
        M_2 = downstream_Mach(M_1)
        p_i2 = p_i1 * pi2_pi1(M_1)
        p_2 = p_1 * P2_P1(M_1) 
        T_2 = T_Ti(M_2) * T_i2
        rho_2 = p_2 / (r_ref * T_2)
        a_2 = sound_velocity(T_2)
        v_2 = M_2 * a_2
        
        h = h_from_temperature(T_1 / (T_ref + 15))
        q_1 = q_from_Mach(M_1, p_1) 
        q_2 =  q_from_Mach(M_2, p_2)
        
        print("h       : ", h)
        Theta, delta, sigma = atmos(h)
        
        # OUTPUTS: 
        display_atmosphere(m2ft(h), delta, sigma, Theta)
   
        print("\n Downstream the shock wave: \n")
       
        display_state(msg='state 2', p=p_2, rho=rho_2, T=T_2, a=a_2)
        display_state_plus(Mach=M_2, T_i=T_i2, p_i=p_i2, u=v_2)
        display_dynamic_pressure(M=M_2, q=q_2, qc=p_i2 - p_2)
        
        print("\n Upstream the shock wave: \n")
        
        display_state(msg='state 1', p=p_1, rho=rho_1, T=T_1, a=a_1)
        display_state_plus(Mach=M_1, T_i=T_i1, p_i=p_i1, u=v_1)
        display_dynamic_pressure(M=M_1, q=q_1, qc=p_i1 - p_1)
        display_pitot_tube(M_1, p_i2 / p_1, delta_p, CAS_1)
        
    elif option == 7:
        plot_pitot_rayleigh(opt=1, gamma=1.4)
        print("Delta_p / p_0 (M=1)  %.3f" % Delta_p_Pitot(1))
        delta_p = 4
        print(M1_from_Delta_p(r=delta_p, gamma=1.4))

    
# *******************************************
#  MAIN
# *******************************************


if __name__ == '__main__':
    main(option=1)
    




 