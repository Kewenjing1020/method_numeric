# -*- coding: utf-8 -*-
"""
Created on Tue May 17 22:41:50 2016

@author: kewenjing
"""

from pylab import *
import time
from mpl_toolkits.mplot3d.axes3d import Axes3D



#boundarty condition of steel,1D
def setBC_steel1(T,Tg,k1,k2,k3,M):
    
    for i in range(M-1):
        T[i,0]=(T[i,1]+Tg*k1)/(1+k1) #left side of steel
        T[i,M-1]=(T[i,M-2]+Tg*k2)/(1+k2) #rignt side
        T[M-1,i]=(T[M-2,i]+Tg*k3)/(1+k3) #top side
    return T

def setBC_steel2(T,value,M):
    for i in range (M):
        T[0,i]=400  #down side of steel
    return T
    
 #boundarty condition of steel,2D  
def setBC_steel3(T,Tg,k1,k2,k3,M):
    for i in range (1,M-1):
        T[i,0]=(T[i,1]+T[i-1,0]+T[i+1,0]+Tg*k1)/(3+k1) #left side of steel
        T[i,M]=(T[i,M-1]+T[i-1,M]+T[i+1,M]+Tg*k2)/(3+k2) #rignt side
        T[M,i]=(T[M-1,i]+T[M,i-1]+T[M,i+1]+Tg*k3)/(3+k3) #top side
#    T[M,0]=(T[M,1]+T[M-1,0]+Tg*k1)/(2+k1)
#    T[M,M]=(T[M,M-1]+T[M-1,M]+Tg*k2)/(2+k2)
#        
#    T[M,0]=((k1+k2)*Tg+T[M-1,0]+T[M,1])/(k1+k2+2)
#    T[0,0]=(T[0,1]+T[1,0]+Tg*k1)/(3+k1) #left side of steel
#    T[0,M]=(T[0,M-1]+T[1,M]+Tg*k2)/(3+k2) #rignt side
#        
    T[0,0]=(k1*Tg+T[0,1])/(1+k1)
    return T
    
def phi_interface(phi,T,M,lamda_s,dx):
    for i in range (M-2):
        phi[i]=lamda_s*(T[1,i]+T[0,i-1]+T[0,i+1]-3*T[0,i])/dx
    return phi
    
#heat equation of steel
def heatEq_Steel(T,Tg,k1,k2,k3,step,w,value,M): 
 
    T=setBC_steel2(T,value,M)
    for k in range (step):
        T=setBC_steel2(T,value,M)
        T=setBC_steel3(T,Tg,k1,k2,k3,M)
        
        for i in range (1,M):
                for j in range(1,M):
                    temp = (T[i+1,j] + T[i-1,j] + T[i,j+1] +T[i,j-1])/4 
                    T[i,j]=(1-w)*T[i,j]+w*temp
    return T
                              


#calculate the velocity of water u(x,y)
def water_velocity(u,Uin,dy):
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

def heatEq_water(T,h,u,step,phi,lamda_w,dx,aw):
    m,n=T.shape
       
    for k in range (step):
        #boundary condition
        T=setBC_water1(T)
        T=setBC_water2(T,phi,lamda_w,dx,n,aw)
        for i in range(m-1):
            for j in range(n-1):
                temp=u[i,j]*aw/h
                T[i,j]=(temp*(T[i+1,j]+T[i-1,j]+T[i,j+1]+T[i,j-1]) +T[i,j-1])/(1+4*temp) 
    return T

#boundary condition of water
#inlet temperature = 300K
def setBC_water1(T):
    m,n=T.shape
    for i in range (1,m):
        T[i,0]=300
    return T

def setBC_water2(T,phi,lamda_w,dx,M):
    k=dx/lamda_w
    for i in range (1,M-2):
        T[0,i]=(phi[i]*k+T[0,i-1]+T[0,i]+T[0,i+1])/3    #interface side
#        T[M/2-1,i]=(T[M/2-1,i-1]+T[M/2-2,i]+T[M/2-1,i+1])/3 #down side
#    
#    for i in range (1,int(M/2)-2):
#        T[i-1,M-1]=(T[i-1,M-1]+T[i+1,M-1]+T[i,M-2])/3 #right side
    return T

#def setbBC3(T):
#    m,n =T.shape
#    for i in range (m):
#        
  
def interface_steel(Tsteel,Twater,h,lamda_s,lamda_w):
    m1,n1=Tsteel.shape
    m2,n2=Twater.shape  
    if (n1!=n2):
        print ("error: dimension of Tsteel and Twater not coherent")
        return NULL
    for i in range (n1-1):
        Tsteel[0,i]=Twater[0,i]
        Tsteel[1,i]=-lamda_w/lamda_s*(Twater[1,i]-Twater[0,i])+Tsteel[0,i]  
    return Tsteel

def interface_water(Tsteel,Twater,h,lamda_s,lamda_w):
    m1,n1=Tsteel.shape
    m2,n2=Twater.shape
    if (n1!=n2):
        print ("error: dimension of Tsteel and Twater not coherent")
        return NULL
    for i in range (n1-1):
        Twater[0,i]=Tsteel[0,i]
        Twater[1,i]=-lamda_s/lamda_w*(Tsteel[1,i]-Tsteel[0,i])+Twater[0,i] 
    return Twater
