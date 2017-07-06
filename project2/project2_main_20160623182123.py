 # -*- coding: utf-8 -*-
"""
Created on Tue May 17 20:41:24 2016

@author: kewenjing
"""

from pylab import *
import project2_function as func
close('all')


M=60
E=2*10**(-2)
e=1*10**(-3)
dx=E/M
dy=2*e/M
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

w=2*(1-pi/M+(pi/M)**2)

Uin=1.04 #to change

# calculate the velocity
u=zeros((int(M/2)+1,M+1))
u=func.water_velocity(u,Uin,dy,e)


Tsteel=zeros((M+1,M+1))
Twater=zeros((int(M/2)+1,M+1))

 # calculate the steel part
Tsteel=func.initial(Tsteel,350)
Tsteel=func.heatEq_Steel(Tsteel,Tg,k1,k2,k3,300,w,300,M)

# calculate the heat flux of interface
phi=zeros(M+1)
phi=func.phi_interface(phi,Tsteel,M,lamda_s,dx)

# calculate the water part
Twater=func.initial(Twater,300)
Twater=func.setBC_water2(Twater,phi,lamda_w,dy)
Twater=func.heatEq_water(Twater,h,u,300,phi,lamda_w,dx,dy,aw)

# calculate the reversal case
Twater_reversal=zeros((int(M/2)+1,M+1))
Twater_reversal=func.initial(Twater_reversal,300)
Twater_reversal=func.setBC_water_reversal(Twater_reversal,phi,lamda_w,dy)
Twater_reversal=func.heatEq_water_reversal(Twater_reversal,h,u,300,phi,lamda_w,dx,dy,aw)

## power heat
power = func.power_heat(phi,M,dx)
print(power)

print(Twater[30,60])

# print the image

#imshow(Tsteel,origin='lower')
#colorbar()
#show()
#
imshow(Twater,origin='lower')
colorbar()
show()

#imshow(Twater_reversal,origin='lower')
#colorbar()
#show()