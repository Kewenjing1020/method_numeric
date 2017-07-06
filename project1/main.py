# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 23:49:17 2016

@author: kewenjing
"""


from BZ_Implicit import *
from pylab import *
from math import *

y0=zeros(2)
y0[1]=0.25
y0[0]=0.5

f=2/3
e=0.4

t0=0.0
tEnd=50.0

time,solBackward=EulerBackward_BZ(y0,t0,tEnd,0.00001,f,e)

fig=plt.figure()
plot(time,solBackward[:,0],label='Backforward Euler')
xlabel(r'$t [s]$',fontsize=20)
ylabel(r'$x [m]$',fontsize=20)
legend()
show()


fig=plt.figure()
plot(time,solBackward[:,1],label='Backforward Euler')
xlabel(r'$t [s]$',fontsize=20)
ylabel(r'$z [m]$',fontsize=20)
legend()
show()


fig=plt.figure()
plot(solBackward[:,0],solBackward[:,1],label='Backforward Euler')
xlabel(r'$t [s]$',fontsize=20)
ylabel(r'$z [m]$',fontsize=20)
legend()
show()