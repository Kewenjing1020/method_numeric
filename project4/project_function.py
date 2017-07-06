# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 14:39:09 2016

@author: Administrator
"""

# project3: Temperatture mixing in a lid driven cavity

from pylab import *
#from TDMA import *


nu_a=6*10**(-6)
Re=50
L=0.06   # m
U0=Re*nu_a/L
Pr=7.0
Dth=nu_a/Pr
tau_lid=0.1
Utarget=Re*nu_a/L
rho=1000

N=30
u=zeros((N+1,N+1))
v=zeros((N+1,N+1))
P=zeros((N+1,N+1))

dx=L/N
FO=0.2
dt1=FO*(dx**2)/nu_a

F=0.2
dt2=2*F/nu_a

F_T=0.5
dt3=F_T*dx**2/Dth

dt=min(dt1,dt2,dt3)

F=nu_a*dt/2
FO=nu_a*dt/(dx**2)
F_T=Dth*dt/(dx**2)



def set_BC(u): 
    for i in range (N+1):
        u[0,i]=0
        u[N,i]=0
        u[i,0]=0
        u[i,N]=0       
    return u


def set_BC_U0(u,U0):
    for i in range (N+1):
        u[N,i]=U0
    return u
    

def diffusion(u,v,step):
    for n in range (step):
        for i in range (1,N):
            for j in range (1,N):
                u[j,i]=u[j,i]+FO*(u[j,i+1]-2*u[j,i]+u[j,i-1]) +FO*(u[j+1,i]-2*u[j,i]+u[j-1,i])
                v[j,i]=v[j,i]+FO*(v[j,i+1]-2*v[j,i]+v[j,i-1]) +FO*(v[j+1,i]-2*v[j,i]+v[j-1,i])
    
    return u,v

# to create the a,b in method TDMA
def createVect(value,dim):
    x = zeros (dim)
    for i in range (dim):
        x[i]=value
    return x
    
def cal_d_TDMA(array,d0,dN):
    d=zeros(N+1)
    d[0]=d0
    d[N]=dN
    for i in range (1,N):
        d[i]=F/2*array[i-1]+(1-F)*array[i]+F/2*array[i+1]
    return d
    
def cal_coeff_TDMA(array,F,d0,dN):
    a=createVect(-F/2,N)
    a[N-1]=0
    b=createVect(1-F,N+1)
    b[0]=1
    b[N]=1
    c=createVect(-F/2,N)
    c[N-1]=0
    d=cal_d_TDMA(array,d0,dN)
    return a,b,c,d
    

def convection(u,v,step):
    for n in range (step):
        for i in range (1,N):
            for j in range (1,N):
                u[j,i]=u[j,i]-dt/(2*dx)*(u[j,i+1]**2-u[j,i-1]**2)-dt/(2*dx)*(u[j+1,i]*v[j+1,i]-u[j-1,i]*v[j-1,i])
                v[j,i]=v[j,i]-dt/(2*dx)*(v[j+1,i]**2-v[j-1,i]**2)-dt/(2*dx)*(u[j,i+1]*v[j,i+1]-u[j,i-1]*v[j,i-1])
    return u,v

       
def set_BC_P(mat):
    for i in range (N+1):
        mat[i,0]=mat[i,1]
        mat[i,N]=mat[i,N-1]
        mat[0,i]=mat[1,i]
        mat[N,i]=mat[N-1,i]
    return mat

def rescale_P(P):
    ave=sum(P)/(N+1)**2
    for j in range(N+1):
        for i in range (N+1):
            P[j,i] -= ave
    return P
    

def divergence_free(P,u,v):
    set_BC_P(P)
    for j in range (1,N):
        for i in range (1,N):
            R_side = ( (u[j,i+1]-u[j,i-1])/(2*dx)+(v[j+1,i]-v[j-1,i])/(2*dx) )/dt
            P[j,i] = (P[j,i+1]+P[j,i-1]+P[j+1,i]+P[j-1,i]-R_side*dx**2)/4
    return P
    
def variance(T):
    sum1=0
    sum2=0
    for i in range (N+1):
        for j in range (N+1):
            sum1 += T[i,j]**2
            sum2 += T[i,j]
    ave= sum2/(N+1)**2
    var= abs(sum1-(sum2**2))/(N+1)**2
    return ave,var    


def update(u,v,P):
    for j in  range (1,N):
        for i in range (1,N):
            u[j,i] = u[j,i]-dt*(P[j,i+1]-P[j,i-1])/(2*dx)
            v[j,i] = v[j,i]-dt*(P[j+1,i]-P[j-1,i])/(2*dx)
    return u,v




def T_initial(T):
    for i in range (N+1):
        for j in range (int(N/2)+1):
            T[j,i]=300
        for j in range (int(N/2)+1,N+1):
            T[j,i]=330
    
def set_BC_T(T):
    for i in range(N+1):
        T[N,i]=T[N-1,i]
        T[0,i]=T[1,i]
        T[i,N]=T[i,N-1]
        T[i,0]=T[i,1]
    return T
    
def scalar_eq(T,u,v,step):
    ave,var=variance(T)
    for i in range (step):
        set_BC_T(T)
        for i in range (1,N):
            for j in range (1,N):
#                T[i,j]=T[i,j]+dt*( Dth/dx**2*(T[i+1,j]+T[i-1,j]+T[i,j+1]+T[i,j-1]-4*T[i,j]) - u[i,j]/dx*(T[i,j+1]-T[i,j]) - v[i,j]/dx*(T[i+1,j]-T[i,j]) )   
                T[i,j]=T[i,j]+dt*( Dth/dx**2*(T[i+1,j]+T[i-1,j]+T[i,j+1]+T[i,j-1]-4*T[i,j]) - cal_upwind(u[i,j],T[i,j-1],T[i,j],T[i,j+1],dx) - cal_upwind(v[i,j],T[i-1,j],T[i,j],T[i+1,j],dx) ) 
        ave,var=variance(T)
        if (var<1):
            break
        
    
    return T,ave,var

def cal_upwind(u,t0,t1,t2,dx):
    if(u<0):
        sol=u*(t2-t1)/dx
    else:
        sol=u*(t1-t0)/dx
    return sol
    


def cal_Ut(Utarget, step ):
    Ut= Utarget*(1-exp(-dt*step/tau_lid))
    return Ut


def power(u,Ut):
    power=0
    for i in range (N+1):
        power += rho* nu_a* abs(u[i,0]-u[i,1])*Ut
    return power

def energy(E,u,Ut):
    pow=power(u,Ut)
    E += pow * dt
    return E
        

#x=linspace(0,L,N+1)
#y=linspace(0,L,N+1)
#quiver(x,y,u,v)

