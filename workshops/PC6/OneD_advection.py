# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:05:36 2016

@author: wenjing ke
"""

# workshop 05/13
# solve linear advection equation

#import numpy as np
from pylab import *
#import matplotlib as plt 

x0=1
sigma=0.2
L=12
c=10

def advection_fE_1stOrder(x0,sigma,L,c,dx,tFin):
    dt=0.8*dx/c
    n1=(int)(L/dx)  #space
    n2=(int)(tFin/dt)   #time
    l=dt*c/dx
    sol=np.zeros((n1,n2))
    
    for i in range(n1-1):
        sol[i,0]=exp(-(i*dx)**2/(sigma)**2)   
#        
#    for j in range (n2-1):
#        sol[0,j]=0
#        for i in range (int(c*j*dt/dx),n1-1):
#            sol[i,j]=sol[(int)(i-c*j*dt/dx),0]
    
    for n in range(n2-1):
        for i in range(1,n1-1):
            sol[i,n+1]=sol[i,n]-l*(sol[i,n]-sol[i-1,n])
    
    return n1,n2,l,sol

dx=0.2
fig=figure()
n1,n2,l,sol= advection_fE_1stOrder(x0,sigma,L,c,dx,10)
x=linspace(x0,n1*dx,n1)
for i in range(0,1):  
    plot(x,sol[:,i])
legend()
show()

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    