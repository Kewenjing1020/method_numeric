# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 16:29:38 2016

@author: Administrator
"""

# workshop 09, incompressible Navier Stokes equations- Explicit method
from pylab import *
from TDMA import *

mu_a=1.716e-5
rho_a=1.225
nu_a=mu_a/rho_a
Re=100
L=0.1   # m
U0=Re*nu_a/L
gamma=nu_a


N=50
u=zeros((N+1,N+1))
v=zeros((N+1,N+1))
P=zeros((N+1,N+1))

dx=L/N
FO=0.2
dt1=FO*(dx**2)/gamma

F=0.2
dt2=2*F/gamma
dt=min(dt1,dt2)
F=gamma*dt/2
FO=gamma*dt/(dx**2)


def set_BC(u): 
    for i in range (N+1):
        u[0,i]=0
        u[i,0]=0
        u[i,N]=0       
    return u


def set_BC_U0(u,U0):
    set_BC(u)
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
    
def diffusion_free_ADI(u,v,step):
    v_temp=zeros((N+1,N+1))
    u_temp=zeros((N+1,N+1))
    u_temp=set_BC_U0(u_temp,U0)
    for n in range (step):
        # x direction
        for i in range (1,N):
             # calculate each row in u
            a,b,c,d=cal_coeff_TDMA(u[i,:],F,0,0)
            u_temp[i,:]=TDMASolve(a, b, c, d)

             # calculate each row in v
            a,b,c,d=cal_coeff_TDMA(v[i,:],F,0,0)
            v_temp[i,:]=TDMASolve(a, b, c, d)
        # y direction
        for i in range (1,N):
            # calculate each column in u
            a,b,c,d=cal_coeff_TDMA(u_temp[:,i],F,0,U0)
            u[:,i]=TDMASolve(a, b, c, d)
            # calculate each column in v
            a,b,c,d=cal_coeff_TDMA(v_temp[:,i],F,0,0)
            v[:,i]=TDMASolve(a, b, c, d)
            
    return u,v

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
    

def divergence_free(u,v):
    set_BC_P(P)
    for j in range (1,N):
        for i in range (1,N):
            R_side = ( (u[j,i+1]-u[j,i-1])/(2*dx)+(v[j+1,i]-v[j-1,i])/(2*dx) )/dt
            P[j,i] = (P[j,i+1]+P[j,i-1]+P[j+1,i]+P[j-1,i]-R_side*dx**2)/4
    return P

def update(u,v,P):
    for j in  range (1,N):
        for i in range (1,N):
            u[j,i] = u[j,i]-dt*(P[j,i+1]-P[j,i-1])/(2*dx)
            v[j,i] = v[j,i]-dt*(P[j+1,i]-P[j-1,i])/(2*dx)
    return u,v

step = 50
v=set_BC(v)
u=set_BC_U0(u,U0)

for j in range (30):
    print(j)
    convection(u,v,1)
    diffusion(u,v,1)  
#    diffusion_free_ADI(u,v,1)
    for i in range (50):
        P=divergence_free(u,v)
        rescale_P(P)
    update(u,v,P)
    v=set_BC(v)
    u=set_BC_U0(u,U0)
    


imshow(u,origin='lower')
colorbar()
show()

imshow(v,origin='lower')
colorbar()
show()

imshow(P,origin='lower')
colorbar()
show()

x=linspace(0,L,N+1)
y=linspace(0,L,N+1)
quiver(x,y,u,v)

