# -*- coding: utf-8 -*-
"""
Created on Wed May 18 13:52:03 2016

@author: Roxane
"""

from pylab import *
import matplotlib as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from math import *

#Paramètres du Problème
D=2*10**-5
A=8000
u=0.3
ErrorMax=0.1
e=D/u #epaisseur de la flamme calculée avec D/u

def f(c,i,n,deltax):
    return D*(c[i+1,n]+c[i-1,n]-2*c[i,n])+A*c[i,n]**2*(1-c[i,n])-u*(c[i,n]-c[i-1,n])/deltax 
    
def progress(t_end):
    #Grille espace
    Nx=100
    nx=Nx+1
    deltax=e 
    x=np.linspace(0,deltax*Nx,nx)
    
    #Grille temps
    C=0.5
    deltat=deltax*C/u
    Nt=int(t_end/deltat)
    nt=Nt+1
    time=np.linspace(0,deltat*Nt,nt)
    c=np.ones((nx,2)) 
    
    #conditions initiales
    
    c[0,0] = 0
    c[-1,0] = 1    
    
    x = linspace(-3*e,3*e,nx)
    for i in range(1,nx-1):
        c[i,0] = np.arctan(150000*x[i])/math.pi+0.5
        
    #forward euler
    for t in range(nt-1):
        print(t)
        for i in range(1,nx-1):
#FE :             
#            c[i,1]=c[i,0]+deltat*(D*(c[i+1,0]+c[i-1,0]-2*c[i,0])+A*c[i,0]**2*(1-c[i,0])-u*(c[i,0]-c[i-1,0])/deltax)       
#RK4 :
            k1=f(c,i,0,deltax)
            k2=f(c+deltat/2*k1,i,0,deltax)
            k3=f(c+deltat/2*k2,i,0,deltax)
            k4=f(c+deltat*k3,i,0,deltax)
            c[i,1]=c[i,0]+deltat/6*(k1+2*k2+2*k3+k4)
        c[:,0]=c[:,1]
    return c[:,0]

Nx=100
nx=Nx+1
deltax=e
x=np.linspace(0,deltax*Nx,nx)
progress=progress(0.011)

plot(x,progress)
xlabel('X')
ylabel('Progress variable C')
title('Flame Propagation')