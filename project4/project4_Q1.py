# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 17:49:48 2016
@author: kewenjing

main_function of Question #1
"""


from pylab import *
from project_function import *

u=zeros((N+1,N+1))
v=zeros((N+1,N+1))
T=zeros((N+1,N+1))

T_initial(T)

step=300    # timestep of T
T,ave,var=scalar_eq(T,u,v,step)
tFin=step*dt
imshow(T)
colorbar()