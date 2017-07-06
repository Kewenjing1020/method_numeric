# -*- coding: utf-8 -*-
"""
Created on Mon May 23 21:17:35 2016

@author: Christopher
"""


#Runge Kutta methods (2 and 4 implemented)

from pylab import *
from initial_BC import *

#definition of dc/dt
def f(c,dx,u,D,A,):
    dim = size(c)
    tmp = zeros(dim)
    for i in range(1,dim-1):
        tmp[i] = (-u/dx*(c[i]-c[i-1])+D/(dx**2)*(c[i+1]-2*c[i]+c[i-1])+A*(c[i]**2)*(1-c[i]))
    tmp = BC_zeroGrad(tmp) #zero grad applies to the derivative BC aswell
    return tmp

#implement Runge Kutta 2 (mid-point method)
def halfPointStep(c,dt,dx,u,D,A):
    dim = size(c)
    k1 = f(c,dx,u,D,A)
    k2 = f(c+dt*k1/2,dx,u,D,A)
    return c+dt*k2

#implement Runge Kutta 4 (standard RK)
def RK4Step(c,dt,dx,u,D,A):
    dim = size(c)
    k1 = f(c,dx,u,D,A)
    k2 = f(c+dt*k1/2,dx,u,D,A)
    k3 = f(c+dt*k2/2,dx,u,D,A)
    k4 = f(c+dt*k3,dx,u,D,A)
    
    return c + dt/6*(k1+2*k2+2*k3+k4)

def halfPoint(c,dt,dx,u,D,A,nSteps):
    for k in range(nSteps):
        c = halfPointStep(c,dt,dx,u,D,A)
        
    return c
    
#def RK4(c,dt,dx,u,D,A,nSteps):
#    for k in range(nSteps):
#        c = halfPointStep(c,dt,dx,u,D,A)
#
#    return c

def RK4(c0,dt,dx,u,D,A,dimTime):     
    dimSpace=size(c0)
    sol=zeros((dimTime,dimSpace))
    sol[0,:]=c0
    for n in range(1,dimTime):
        sol[n,:]=halfPointStep(sol[n-1,:],dt,dx,u,D,A)
     
    return sol      
    
    