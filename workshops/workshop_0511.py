# -*- coding: utf-8 -*-
"""
Created on Wed May 11 15:59:26 2016

@author: kewenjing
"""
import numpy as np
from pylab import *
import time
from mpl_toolkits.mplot3d.axes3d import Axes3D

#1d unsteady difussion

L=1
a=0.0001

def forwardEuler_2nd(a,L,dx,dt,tEnd):
    N1=(int)(L/dx)
    N2=(int)(tEnd/dt)
    F0=a*dt/(dx)**2
    print (N1)
    print (N2)
    
    sol = np.zeros((N1,N2))

    for i in range (0,N2-1):
        sol[0,i]=0
        sol[N1-1,i]=0
    
    for i in range(0,N1-2):
        sol[i,0]=sin(pi*i*dx/L)
        for j in range (0,N2-2):
            sol[i,j+1]=sol[i,j]+F0*(sol[i+1,j]-2*sol[i,j]+sol[i-1,j])
            
    return F0,sol
    
t0=time.clock()
F0,forwardEuler_2nd_sol=forwardEuler_2nd(a,L,0.01,0.001,1)

fig=plt.figure()
x=linspace(0,1,100)
for i in range (0,1000):
    plot(x,forwardEuler_2nd_sol[:,i],label='forwardEuler_2nd')

            
    
    
    