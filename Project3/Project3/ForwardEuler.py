# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:26:09 2016

@author: Christopher
"""

from pylab import *
from initial_BC import *
close('all')

#convection diffusion equation + source term
#forward euler + upwind (advection) + centered 2nd order (diffusion)

def forwardEulerStep(c,dt,dx,u,D,A):
    dim = size(c)
#    tmp = zeros(dim)
    for i in range (1,dim-1):
        c[i] = c[i] + dt*(-u/dx*(c[i]-c[i-1])+D/(dx**2)*(c[i+1]-2*c[i]+c[i-1])+A*(c[i]**2)*(1-c[i]))
        
    BC_zeroGrad(c)
        
    return c
    
def forwardEuler(c,dt,dx,u,D,A,k): #k is number of steps
    for i in range(k):
        c = forwardEulerStep(c,dt,dx,u,D,A)
    
    return c
    
def forwardEuler_1stupwind(c0,dt,dx,u,D,A,dimTime):
    dimSpace=size(c0)
    sol=zeros((dimTime,dimSpace))
    sol[0,:]=c0
    for n in range(dimTime-1):
        BC_zeroGrad(sol[n,:])
        for i in range (1,dimSpace-1):
            sol[n+1,i]=sol[n,i]+ dt*(-u/dx*(sol[n,i]-sol[n,i-1])+D/(dx**2)*(sol[n,i+1]-2*sol[n,i]+sol[n,i-1])+A*(sol[n,i]**2)*(1-sol[n,i]))
     
    return sol      
            
def forwardEuler_centered2nd(c0,dt,dx,u,D,A,dimTime):
    dimSpace=size(c0)
#    dimTime=(int)(tFin/dt)
    sol=zeros((dimTime,dimSpace))
    sol[0,:]=c0
    for n in range(dimTime-1):
        BC_zeroGrad(sol[n,:])
        for i in range (1,dimSpace-1):
            sol[n+1,i]=sol[n,i]+ dt*(-u/dx*(sol[n,i+1]-sol[n,i-1])/2+D/(dx**2)*(sol[n,i+1]-2*sol[n,i]+sol[n,i-1])+A*(sol[n,i]**2)*(1-sol[n,i]))
     
    return sol      
    
def forwardEuler_1stdownwind(c0,dt,dx,u,D,A,dimTime):
    dimSpace=size(c0)
#    dimTime=(int)(tFin/dt)
    sol=zeros((dimTime,dimSpace))
    sol[0,:]=c0
    for n in range(dimTime-1):
        BC_zeroGrad(sol[n,:])
        for i in range (1,dimSpace-1):
            sol[n+1,i]=sol[n,i]+ dt*(-u/dx*(sol[n,i+1]-sol[n,i])+D/(dx**2)*(sol[n,i+1]-2*sol[n,i]+sol[n,i-1])+A*(sol[n,i]**2)*(1-sol[n,i]))
     
    return sol   
