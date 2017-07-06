# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 23:48:09 2016

@author: kewenjing
"""

import numpy as np
from pylab import *
from scipy.optimize import *
from numpy import *

q= 8*10**(-4)
 
def EulerForward(y0,t0,tEnd,dt,f,e):
      N=int ((tEnd-t0)/dt)+1
      time=np.linspace(t0,tEnd,N)
     
      sol=np.zeros((N,np.size(y0)))
      sol[0]=y0
     
      for i in range (N-1):
          sol[i+1,0]=sol[i,0]+ dt/e*(sol[i,0]*(1-sol[i,0])+f*(q-sol[i,0])/(q+sol[i,0])*sol[i,1] )
          sol[i+1,1]=sol[i,1]+ dt*(sol[i,0]-sol[i,1])
     
      return time, sol
      
def EulerBackward(f,e,y0,t0,tEnd,dt):
    N=int ((tEnd-t0)/dt)+1
    time=np.linspace(t0,tEnd,N)
    Y=zeros((N,2))
    Y[0]=y0
    for i in range(1,N):
        def G(x):
            x0=x[0]
            x1=x[1]
#            return [Y[i-1][0]-x0+dt*((x0*(1-x0)+f*(q-x0)/(q+x0)*x1)), Y[i-1][1]-x1+dt*(x0-x1)]
            return [Y[i-1][0]-x0+dt/e*((Y[i-1][0]*(1-Y[i-1][0])+f*(q-Y[i-1][0])/(q+Y[i-1][0])*Y[i-1][1])), Y[i-1][1]-x1+dt*(Y[i-1][0]-Y[i-1][1])]
        (Y[i][0],Y[i][1])=fsolve(G,(0,0))
    return time,Y


def Trapezoidal(f,e,y0,t0,tEnd,dt):
    N=int ((tEnd-t0)/dt)+1
    time=np.linspace(t0,tEnd,N)
    Y=zeros((N,2))
    Y[0]=y0
    for i in range(1,N):
        def G(x):
            x0=x[0]
            x1=x[1]
            temp1=( (Y[i-1][0]*(1-Y[i-1][0])+f*(q-Y[i-1][0])/(q+Y[i-1][0])*Y[i-1][1])+(x0*(1-x0)+f*(q-x0)/(q+x0)*x1) )/2
            temp2=( (x0-x1)+ Y[i-1][0]-Y[i-1][1] )/2
            return [ Y[i-1][0]-x0+dt/e*temp1, Y[i-1][1]-x1+dt*temp2 ]
        (Y[i][0],Y[i][1])=fsolve(G,(0,0))
    return time,Y
            
            
            
            