# -*- coding: utf-8 -*-
"""
Created on Tue May 17 20:41:24 2016

@author: kewenjing
"""

from pylab import *
import project2_function as func
close('all')


M=50
E=2
e=1
dx=E/M
dy=E/M
h=dx

lamda_s=50
lamda_w=0.6

hl=300
ht=30
hr=100

Tg=700
Tin=300

aw=10**(-6)

k1=dx*hl/lamda_s
k2=dx*hr/lamda_s
k3=dy*ht/lamda_s

w=2*(1-pi/M)

Uin=1 #to change
u=zeros((M/2+1,M+1))
u=func.water_velocity(u,Uin,dy)

Tsteel=zeros((M+1,M+1))
Twater=zeros((M/2+1,M+1)) 


Tsteel=func.initial(Tsteel,600)

Tsteel=func.heatEq_Steel(Tsteel,Tg,k1,k2,k3,20,w,300,M)

phi=zeros(M)
phi=func.phi_interface(phi,Tsteel,M,lamda_s,dx)


xArray = linspace(0,2,M+1)
yArray = linspace(0,2,M+1)
yArray2 = linspace(0,2,int(M/2)+1)
X,Y = meshgrid(xArray,yArray,indexing='ij')
X2,Y2= meshgrid(yArray2,xArray,indexing='ij')

fig = figure(figsize=(10,10))
axes = fig.add_subplot(111,projection='3d')
pl = axes.plot_surface(X,Y,Tsteel,rstride=1,cstride=1,cmap=cm.jet)
cb = fig.colorbar(pl)

Twater=func.initial(Twater,300)
Twater=func.setBC_water2(Twater,phi,lamda_w,dx,M )
#Twater=func.heatEq_water(Twater,h,u,20,phi,lamda_w,dx,aw)

fig2 = figure(figsize=(10,10))
axes = fig2.add_subplot(111,projection='3d')
#plots
pl2 = axes.plot_surface(X2,Y2,Twater,rstride=1,cstride=1,cmap=cm.jet)
cb2 = fig.colorbar(pl2)


    