# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 22:45:08 2016
@author: kewenjing

main_function of Question #3
change the moving lid: we have the moving lid both the upper surface and bottom: Ut,Ub 
"""
from pylab import *
from project_function import *

# change the boundary condition
def set_BC_Q3(u,Ut,Ub):
    for i in range (N+1):
        u[0,i]=Ub
        u[N,i]=Ut
    return u

Ut=U0
Ub=-U0

u=zeros((N+1,N+1))
v=zeros((N+1,N+1))
P=zeros((N+1,N+1))
T=zeros((N+1,N+1))

T_initial(T)

v=set_BC(v)
u=set_BC(u)
u=set_BC_Q3(u,Ut,Ub)

E=0 #energy

for step in range (300):    # time step
    convection(u,v,1)
    diffusion(u,v,1)  
    for i in range (5):     # time step of P, do the iteration to get stable P
        P=divergence_free(P,u,v)
    rescale_P(P)
    update(u,v,P)
    v=set_BC(v)
    set_BC(u)
    set_BC_Q3(u,Ut,Ub)
    tFin=step*dt
    T,ave,var=scalar_eq(T,u,v,1)
    E=energy(E,u,Ut)
    while (var<0.1*315**2):
        break
    
x=linspace(0,L,N+1)
y=linspace(0,L,N+1)
quiver(x,y,u,v)
#  
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

imshow(T,origin='lower')
colorbar()
show()

print(E)
print(tFin)
print(var)