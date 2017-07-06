# -*- coding: utf-8 -*-
"""
Created on Sat May 21 13:57:44 2016

@author: Ronan Vicquelin
"""

#-------------------------------
# Quartz Material
#-------------------------------

# Density [kg/m3]   `
rho_q = 2200.0

# Molar Weight [kg/mole] 
W_q   = 60.1e-3  

# Thermal Conductivity [W/m/K] - Constant
cond_q_ct = 0.45

# Heat Capacity [J/kg/K] - Constant
c_q_ct   = 1330.0

    
# Thermal Diffusivity [m2/s] - Constant
a_q_ct = 1.55e-7

#-------------------------------
# Air Material
#-------------------------------

# Prandtl Number [-]
Pr_a = 0.72 

# Heat Capacity [J/kg/K] - Constant
cp_a_ct = 1005.0    
    
# Density [kg/m3] - Constant  
rho_a_ct = 1.225
 

# Dynamic Viscosity [kg/m/s] - Constant
mu_a_ct = 1.716e-5
    
# Thermal Diffusivity [m2/s] - Constant
a_a_ct = mu_a_ct/rho_a_ct/Pr_a 
    
# Thermal Conductivity [W/m/K] - Constant
cond_a_ct = mu_a_ct*cp_a_ct/Pr_a