# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 23:49:17 2016

@author: kewenjing
"""


from functions import *
from pylab import *
from scipy.optimize import *
from numpy import *

y0=zeros(2)
y0[1]=0.25
y0[0]=0.5

f=2/3
e=0.01

t0=0.0
tEnd=15.0

dt=0.0001


time,solForward=EulerForward(y0,t0,tEnd,dt,f,e)
fig=plt.figure()
plot(time,solForward[:,0],label='x_Forward Euler')
plot(time,solForward[:,1],label='z_Forward Euler')
xlabel(r'$t [s]$',fontsize=20)
ylabel(r'$x [m]$',fontsize=20)
legend()
show()


time,solBackward=EulerBackward(f,e,y0,t0,tEnd,dt)
fig=plt.figure()
plot(time,solBackward[:,0],label='x_Backward Euler')
plot(time,solBackward[:,1],label='z_Backward Euler')
xlabel(r'$t [s]$',fontsize=20)
ylabel(r'$x [m]$',fontsize=20)
legend()
show()


time,solTrapezoidal=Trapezoidal(f,e,y0,t0,tEnd,dt)
fig=plt.figure()
plot(time,solTrapezoidal[:,0],label='x_Trapezoidal Euler')
plot(time,solTrapezoidal[:,1],label='z_Trapezoidal Euler')
xlabel(r'$t [s]$',fontsize=20)
ylabel(r'$x [m]$',fontsize=20)
legend()
show()

#fig=plt.figure()
#plot(time,solForward[:,1],label='z_Forward Euler')
#plot(time,solBackward[:,1],label='z_Backward Euler')
#plot(time,solTrapezoidal[:,1],label='z_Trapezoidal Euler')
#xlabel(r'$t [s]$',fontsize=20)
#ylabel(r'$z [m]$',fontsize=20)
#legend()
#show()



