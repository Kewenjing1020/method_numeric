# -*- coding: utf-8 -*-
"""
Created on Wed May 25 16:03:55 2016

@author: Administrator
"""

#workshop 07: implicit methods for parabolic PDEs
#1D unsteady diffusion

from pylab import *
from TDMA import *
# sol is 1D array

def BC_1Ddiffusion(sol):
    dim=size(sol)
    sol[0]=0
    sol[dim-1]=0
    
    return sol
    
# T(x,o)=sin(pi*x/L)
def initial_1Ddiffusion(dx,L,dimSpace):
    sol=zeros(dimSpace)
    for i in range(dimSpace):
        sol[i]=sin(pi*i*dx/L)
    return sol

# implicit Crank-Nicolson method and the centered-2nd-order in space
# use TDMA.py
def implicitCN(sol0,F0,dimTime,dimSpace):
    sol=zeros((dimTime,dimSpace))
    sol[0,:]=sol0
    a=createVector(-F0/2,dimSpace-1)
    b=createVector(1+F0,dimSpace)
    c=createVector(-F0/2,dimSpace-1)
    c[0]=0
    b[0]=1
    a[-1]=0
    b[dimSpace-1]=1
    for n in range (dimTime-1):
#        sol[n,:]=BC_1Ddiffusion(sol[n,:])
        d=d_TDMA(sol[n,:],F0)
        sol[n+1,:]=TDMASolve(a, b, c, d)
        b=createVector(1+F0,dimSpace)
        b[0]=1
        b[dimSpace-1]=1
    
    return d,sol
    
def d_TDMA(sol,F0):
    dim=size(sol)
    d=zeros(dim)
    d[0]=0
    d[dim-1]=0
    for i in range (1,dim-2):
        d[i]=F0/2*sol[i-1]+(1-F0)*sol[i]+F0/2*sol[i+1]
    
    return d

def createVector(value,dim):
    x = zeros (dim)
    for i in range (dim-1):
        x[i]=value
    return x

      
a=10**(-4)
L=1
dx=0.05

F0=0.1
dt=F0*(dx)**2/a

dimSpace=(int)(L/dx)+1
dimTime=100

sol0=initial_1Ddiffusion(dx,L,dimSpace)
d,sol=implicitCN(sol0,F0,dimTime,dimSpace)

x=linspace(0,dx*dimSpace,dimSpace)
for i in range (dimTime):
    plot(x,sol[i])

