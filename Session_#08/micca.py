# -*- coding: utf-8 -*-
"""
Created on Fri May 27 15:39:38 2016

@author: Administrator
"""
from pylab import *

#def cool_uniform():

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



h=7.9
Ta=300
e=5*10**(-3)
dt=10*2
tFin=10**4
Nt=(int)(tFin/dt)

# 0D, h homogene
def cooling_0D(h,Nt):
    sol=zeros(Nt)
    sol[0]=890
    for i in range (Nt-1):
        sol[i+1]=2*h*dt/(rho_q*e*c_q_ct)*(Ta-sol[i])+sol[i]
    return sol

sol=cooling_0D(h,Nt)
x=linspace(0,tFin,N)
plot(x,sol)

# 1D
def hx(T,x):
    Grx=9.8*(2/(300+890))*(T-Ta)*x**3/(mu_a_ct)**2
    Rcx=Grx*Pr_a
    Nux=3/4*0.515*Rcx**(1/4)
    hx=cond_q_ct*Nux/200

Nx=200
hx=zeros(Nx)


   

