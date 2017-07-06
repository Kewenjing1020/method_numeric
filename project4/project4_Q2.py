# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 18:30:26 2016
@author: kewenjing

main_function of Question #2
"""
from pylab import *
from project_function import *

u=zeros((N+1,N+1))
v=zeros((N+1,N+1))
P=zeros((N+1,N+1))
T=zeros((N+1,N+1))

T_initial(T)

v=set_BC(v)
set_BC(u)
u=set_BC_U0(u,Utarget)

E=0 #energy

for step in range (300):    # time step
    convection(u,v,1)
    diffusion(u,v,1)  
    for i in range (5):     # time step of P
        P=divergence_free(P,u,v)
    rescale_P(P)
    update(u,v,P)
    Ut=cal_Ut(Utarget, step)
    v=set_BC(v)
    u=set_BC_U0(u,Ut)
    tFin=step*dt
    T,ave,var=scalar_eq(T,u,v,1)
    E=energy(E,u,Ut)
    while (var<0.1*315**2):
        break
    
#imshow(u,origin='lower')
#colorbar()
#show()
#
#imshow(v,origin='lower')
#colorbar()
#show()
#
#imshow(P,origin='lower')
#colorbar()
#show()

x=linspace(0,L,N+1)
y=linspace(0,L,N+1)
quiver(x,y,u,v)

imshow(T,origin='lower')
colorbar()
show()

print(E)
print(tFin)
print(var)

