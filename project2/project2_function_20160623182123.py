# -*- coding: utf-8 -*-
"""
Created on Tue May 17 22:41:50 2016

@author: kewenjing
"""

from pylab import *
import time
from mpl_toolkits.mplot3d.axes3d import Axes3D




def setBC_steel4(T,Tg,k1,k2,k3,M):    
    for i in range(M+1):
        T[i,0]=(4*T[i,1]-T[i,2]+Tg*k1*2)/(3+2*k1) #left side of steel
        T[i,M]=(4*T[i,M-1]-T[i,M-2]+Tg*k2*2)/(3+2*k2) #rignt side
        T[M,i]=(4*T[M-1,i]-T[M-2,i]+Tg*k3*2)/(3+2*k3) #top side
    return T

def setBC_steel2(T,value,M):
    for i in range (M+1):
        T[0,i]=300  #down side of steel
    return T

    
def phi_interface(phi,T,M,lamda_s,dx):
    for i in range (M+1):
        phi[i]=lamda_s*(4*T[1,i]-T[2,i]-3*T[0,i])/(2*dx)
    return phi
    
#heat equation of steel
def heatEq_Steel(T,Tg,k1,k2,k3,step,w,value,M): 
 
    T=setBC_steel2(T,value,M)
    for k in range (step):
        T=setBC_steel2(T,value,M)  
        T=setBC_steel4(T,Tg,k1,k2,k3,M)
        for i in range (1,M):
                for j in range(1,M):
                    temp = (T[i+1,j] + T[i-1,j] + T[i,j+1] +T[i,j-1])/4 
                    T[i,j]=(1-w)*T[i,j]+w*temp
        
    return T
                              


#calculate the velocity of water u(x,y)
def water_velocity(u,Uin,dy,e):
    n,m=u.shape
    for j in range(n):
        for i in range(m):
            y=dy*j
            u[j,i]=8*Uin*y/e*(1-y/e)
    return u

def initial(T,value):
    m,n=T.shape
    for i in range(m):
        for j in range (n):
           T[i,j]=value
    return T      

def heatEq_water(T,h,u,step,phi,lamda_w,dx,dy,aw):
    m,n=T.shape
       
    for k in range (step):
        #boundary condition
        T=setBC_water1(T)
        T=setBC_water2(T,phi,lamda_w,dy)
        for i in range(1,m-1):
            for j in range(1,n-1):
                temp=u[i,j]/(aw*dx)
                T[i,j]=((T[i+1,j]+T[i-1,j])/dy**2+(T[i,j+1]+T[i,j-1])/dx**2 +T[i,j-1]*temp)/(2/dx**2+2/dy**2+temp) 
    return T

def heatEq_water_reversal(T,h,u,step,phi,lamda_w,dx,dy,aw):
    m,n=T.shape
       
    for k in range (step):
        #boundary condition
        T=setBC_water_reversal(T,phi,lamda_w,dy)
        for i in range(m-2,0,-1):
            for j in range(n-2,0,-1):
                temp=u[i,j]/(aw*dx)
                T[i,j]=((T[i+1,j]+T[i-1,j])/dy**2+(T[i,j+1]+T[i,j-1])/dx**2 +T[i,j+1]*temp)/(2/dx**2+2/dy**2+temp) 
    return T

#boundary condition of water
#inlet temperature = 300K
def setBC_water1(T):
    m,n=T.shape
    for i in range (m):
        T[i,0]=300  # left side, Tin=300
    return T

def setBC_water2(T,phi,lamda_w,dy):
    m,n=T.shape
    k=dy/lamda_w
    for i in range (1,n):
        T[0,i]=(4*T[1,i]-T[2,i])/3  #bottom side
        T[m-1,i]=(phi[i]*2*k+4*T[m-2,i]-T[m-3,i])/3    #interface side       
    for i in range (m):
        T[i,n-1]=(4*T[i,n-2]-T[i,n-3])/3  #right side

    return T
    
def setBC_water_reversal(T,phi,lamda_w,dy):
    m,n=T.shape
    k=dy/lamda_w
    for i in range (m-1,-1,-1):
        T[i,n-1]=300  # in side, Tin=300
        T[i,0]=(4*T[i,1]-T[i,2])/3  #out side
    for i in range (n-1,-1,-1):
        T[0,i]=(4*T[1,i]-T[2,i])/3  #bottom side
        T[m-1,i]=(phi[i]*2*k+4*T[m-2,i]-T[m-3,i])/3    #interface side 
    return T

def power_heat(phi,M,dx):
    m=phi.shape
    sum=0
    for i in range (M+1):
        sum += phi[i]*dx
    return sum
 