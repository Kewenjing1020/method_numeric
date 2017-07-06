# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 23:48:09 2016

@author: kewenjing
"""

import numpy as np
from pylab import *


def EulerBackward_BZ(y0,t0,tEnd,dt,f,e):
	N=int ((tEnd-t0)/dt)+1
	time=np.linspace(t0,tEnd,N)

	sol=np.zeros((N,np.size(y0)))
	sol[0]=y0


	for i in range(1,N):
         b=1-f*dt/(1+dt)-e/dt
         c=-e/dt*sol[i-1,0]+f*sol[i-1,1]/(1+dt)
         sol[i,0]=(b+sqrt(b**2-4*c))*(1/2)
         sol[i,1]=(sol[i-1,1]+dt*sol[i,0])/(1+dt)

	return time,sol
