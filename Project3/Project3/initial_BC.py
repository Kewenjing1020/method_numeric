# -*- coding: utf-8 -*-
"""
Created on Mon May 23 15:02:46 2016

@author: Christopher
"""
from pylab import *
import numpy
close('all')
#boundary and initial conditions for advection diffusion problem

#set left bc to 0 set right bc to 1
def BC_diri(c):
    dim = size(c)
    c[0] = 0
    c[dim-1] = 1
    
    return c

#zero gradient bc
def BC_zeroGrad(c):
    dim = size(c)
    c[0] = c[1]
    c[dim-1] = c[dim-2]
    
    return c

def initialCond(c,u,D): #returns x array and fitting c array
    dim = size(c)
    t_f = D/u #flame thickness
    x = linspace(-10*t_f,10*t_f,dim)
    for i in range(dim):
        c[i] = numpy.arctan(10*x[i]/t_f)/math.pi+0.5
    
    return x,c
