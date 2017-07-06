# -*- coding: utf-8 -*-
"""
Created on Sat May 21 13:57:44 2016

@author: Ronan
"""

from pylab import *

if __name__ == "__main__":
    close('all')

#-------------------------------
# Quartz Material
#-------------------------------

# Density [kg/m3]   
rho_q = 2200.0

# Molar Weight [kg/mole] 
W_q   = 60.1e-3  

# Thermal Conductivity [W/m/K] - Variable

T_tab    = array([300.0, 500.0, 700.0, 900.0, 1100.0, 1300.0])
cond_tab = array([0.3, 0.34, 0.39,0.45, 0.51,0.62])

def func_cond_q(T):
    return 0.3 + 2.5e-4*(T-300.0)

if __name__ == "__main__":
    fig =figure()
    ax = fig.add_subplot(111)
    ax.set_title("Quartz Thermal Conductivity", fontsize=22,fontname='serif')

    ax.set_xlabel("$T\quad[K]$", fontsize=18)
    ax.set_ylabel("$\lambda_q\quad[W/m/K]$", fontsize=18)
    ax.plot(T_tab, cond_tab,"o", label = "Data")
    ax.plot(T_tab, func_cond_q(T_tab), "--", label ="Fitted")
    ax.legend(loc=1)


# Heat Capacity [J/kg/K] - Variable

def func_c_q(T):
    return (-6.076591 + 251.6755*(T/1000.0) -324.7964*(T/1000.0)**2 + 168.5604*(T/1000.0)**3 + 0.002548/(T/1000.0))/W_q

if __name__ == "__main__":
    T_tab = linspace(300,1000,50)
    c_q_tab = zeros(50)
    for i in range(len(T_tab)):
        c_q_tab[i] = func_c_q(T_tab[i])
        
    fig =figure()
    ax = fig.add_subplot(111)
    ax.set_title("Quartz Heat Capacity", fontsize=22,fontname='serif')

    ax.set_xlabel("$T\quad[K]$", fontsize=18)
    ax.set_ylabel("$c_q\quad[J/kg/K]$", fontsize=18)
    ax.plot(T_tab, c_q_tab)

# Thermal Diffusivity [m2/s] - Variable

def func_a_q(T):
    return func_cond_q(T)/(rho_q*func_c_q(T))

if __name__ == "__main__":
    T_tab = linspace(300,1000,50)
    a_q_tab = zeros(50)
    for i in range(len(T_tab)):
        a_q_tab[i] = func_a_q(T_tab[i])
        
    fig =figure()
    ax = fig.add_subplot(111)
    ax.set_title("Quartz Thermal Diffusivity", fontsize=22,fontname='serif')

    ax.set_xlabel("$T\quad[K]$", fontsize=18)
    ax.set_ylabel("$a_q\quad[m^2/s]$", fontsize=18)
    ax.plot(T_tab, a_q_tab)

#-------------------------------
# Air Material
#-------------------------------

# Prandtl Number [-]
Pr_a = 0.72 

# Heat Capacity [J/kg/K] - Variable

def func_cp_a(T):
    return 1.9327E-10*T**4 - 7.9999E-07*T**3 + 1.1407E-03*T**2 - 4.4890E-01*T + 1.0575E+03

if __name__ == "__main__":
    T_tab = linspace(300,1000,50)
    cp_a_tab = zeros(50)
    for i in range(len(T_tab)):
        cp_a_tab[i] = func_cp_a(T_tab[i])
        
    fig =figure()
    ax = fig.add_subplot(111)
    ax.set_title("Air Heat Capacity", fontsize=22,fontname='serif')

    ax.set_xlabel("$T\quad[K]$", fontsize=18)
    ax.set_ylabel("$c_{p_a} \quad[J/kg/K]$", fontsize=18)
    ax.plot(T_tab, cp_a_tab)
    
# Density [kg/m3] - Variable

def func_rho_a(T):
    return 1.225*(273+15)/T
  
if __name__ == "__main__":
    T_tab = linspace(300,1000,50)
    rho_a_tab = zeros(50)
    for i in range(len(T_tab)):
        rho_a_tab[i] = func_rho_a(T_tab[i])
        
    fig =figure()
    ax = fig.add_subplot(111)
    ax.set_title("Air Density", fontsize=22,fontname='serif')

    ax.set_xlabel("$T\quad[K]$", fontsize=18)
    ax.set_ylabel("$\rho_a \quad[kg/m^3]$", fontsize=18)
    ax.plot(T_tab, rho_a_tab)


# Dynamic Viscosity [kg/m/s] - Variable

def func_mu_a(T):
    return 1.716e-5*((T/273.15)**(1.5))*(273.15+110.4)/(T+110.4)

if __name__ == "__main__":
    T_tab = linspace(300,1000,50)
    mu_a_tab = zeros(50)
    for i in range(len(T_tab)):
        mu_a_tab[i] = func_mu_a(T_tab[i])
        
    fig =figure()
    ax = fig.add_subplot(111)
    ax.set_title("Air Dynamic Viscosity", fontsize=22,fontname='serif')

    ax.set_xlabel("$T\quad[K]$", fontsize=18)
    ax.set_ylabel("$\mu_a\quad[kg/m/s]$", fontsize=18)
    ax.plot(T_tab, mu_a_tab)

# Thermal Diffusivity [m2/s] - Variable

def func_a_a(T):
    return func_mu_a(T)/func_rho_a(T)/Pr_a

if __name__ == "__main__":
    T_tab = linspace(300,1000,50)
    a_a_tab = zeros(50)
    for i in range(len(T_tab)):
        a_a_tab[i] = func_a_a(T_tab[i])
        
    fig =figure()
    ax = fig.add_subplot(111)
    ax.set_title("Air Thermal Diffusivity", fontsize=22,fontname='serif')

    ax.set_xlabel("$T\quad[K]$", fontsize=18)
    ax.set_ylabel("$a_a\quad[m^2/s]$", fontsize=18)
    ax.plot(T_tab, a_a_tab)


# Thermal Conductivity [W/m/K] - Variable
def func_cond_a(T):
    return func_mu_a(T)*func_cp_a(T)/Pr_a
    