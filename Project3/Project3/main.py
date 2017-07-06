# -*- coding: utf-8 -*-
"""
Created on Mon May 23 15:54:15 2016

@author: Christopher
"""

from pylab import *
from ForwardEuler import *
from initial_BC import *
from RungeKutta import *

close('all')

#main

#variables
C = 0.05 #CFL number
k = 100 #number of timesteps

#set fields
N = 100
c = zeros(N)
x = zeros(N)


#set variables
A = 8000 #1/s
D = 2*10**(-5) #m^2/s
u = 0.1 #m/s  
e = D/u #flame thickness (just checking)
dimTime=100

#set initial and bc
x,c = initialCond(c,u,D)
dx = x[1] - x[0]
dt = dx*C/u #cfl condition
c = BC_diri(c)

fig1=figure()
sol1= forwardEuler_1stupwind (c,dt,dx,u,D,A,dimTime)
for i in range (dimTime):
    plot(x,sol1[i,:])
title('forwardEuler_1st upwind, u=0.3')
show()

#fig2=figure()
#sol2=forwardEuler_centered2nd(c,dt,dx,u,D,A,dimTime)
#for i in range (dimTime):
#    plot(x,sol2[i,:])
#title('forwardEuler_centered 2nd, u=0.3,dimTime(500->1000)')
#show()
#
#
#fig3=figure()
#sol3=forwardEuler_1stdownwind(c,dt,dx,u,D,A,dimTime)
#for i in range (dimTime):
#    plot(x,sol3[i,:])
#title('forwardEuler_1st downwind, u=0.3')
#show()
#
#
#fig4=figure()
#sol4=RK4(c,dt,dx,u,D,A,dimTime)
#for i in range (dimTime):
#    plot(x,sol4[i,:])
#title('RK4, u=0.3')
#show()