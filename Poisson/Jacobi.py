# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 14:06:37 2016

@author: Christopher
"""

#implement function f and jacobi discretisation

from pylab import *
import time
from mpl_toolkits.mplot3d.axes3d import Axes3D


close('all')

N = 50 # hx = hy = h
M = N #randomly chose a number of points for the test
xArray = linspace(-1,1,N)
yArray = linspace(-1,1,M)



def f(x,y):
    return 2*(x**2+y**2-2)

def fArray(x,y):
    fSolArray = zeros((size(x),size(y)))
    for i in range(size(x)):
        for j in range(size(y)):    
            fSolArray[i,j] = f(x[i],y[j])
    return fSolArray

fSol = fArray(xArray,yArray)
X,Y = meshgrid(xArray,yArray,indexing='ij')

##3D figure
#fig = figure(figsize=(4,3))
#axes = fig.add_subplot(111,projection='3d')
##plots
#pl = axes.plot_surface(X,Y,fSol,rstride=1,cstride=1,cmap=cm.jet)
#cb = fig.colorbar(pl)

#implement jacobi algorithm now
k = 10 #number of time iterations to be executed
BcValue = 0 #Dirichlet at every boundary

#set boundary condidtion
def setBoundary(gridSol):
    n,m = shape(gridSol)
    for i in range(n-1):
        gridSol[i,0] = BcValue
        gridSol[i,m-1] = BcValue
    
    for j in range(m-1):
        gridSol[0,j] = BcValue
        gridSol[n-1,j] = BcValue
    
    return gridSol

solBound = setBoundary(fSol)

def GaussSeidel(phi,fSol,steps):
    h = 1/N
    n,m = phi.shape 
    
    for k in range(steps):
        setBoundary(phi)
        for i in range (1,n-1):
            for j in range(1,m-1):
                phi[i,j] = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] +phi[i,j-1])/4 - fSol[i,j]/4
            
    return phi
    

    
def Jacobi(phi,fSol,steps): 
    h = 1/N
    n,m = phi.shape 
    temp = zeros(phi.shape)
    
    for k in range(steps):
        phi = temp
        setBoundary(phi)
        for i in range (1,n-1):
            for j in range(1,m-1):
                temp[i,j] = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] +phi[i,j-1])/4 - fSol[i,j]/4
            
    setBoundary(temp)
    return temp
    
    
def compError(phi):
    n,m = phi.shape
    err = zeros(phi.shape)
    absErrSqr = 0
    for i in range (1,n-1):
        for j in range (1,m-1):
            err[i,j] = abs(phi[i,j] - ((phi[i+1,j] + phi[i-1,j] + phi[i,j+1] +phi[i,j-1])/4 - fSol[i,j]/4))
            absErrSqr = absErrSqr + err[i,j]**2
    
    return sqrt(absErrSqr/(n*m))

def JacobiWhile(phi,fSol,maxErr,w): #error has to be chosen smaller than 10
    n,m = phi.shape 
    temp = zeros(phi.shape)
    absErr = 10
    
    while absErr >= maxErr:
        absErr = compError(temp)
        phi = temp
        setBoundary(phi)
        for i in range (1,n-1):
            for j in range(1,m-1):
                a = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] +phi[i,j-1])/4 - fSol[i,j]/4
                temp[i,j]=(1-w)*temp[i,j]+w*a
            
    setBoundary(temp)
    return temp

def GaussSeidelWhile(phi,fSol,maxErr,w):
    n,m = phi.shape 
    absErr = 10
    
    while absErr >= maxErr:
        absErr = compError(phi)
        setBoundary(phi)
        for i in range (1,n-1):
            for j in range(1,m-1):
                temp = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] +phi[i,j-1])/4 - fSol[i,j]/4
                phi[i,j]=(1-w)*phi[i,j]+w*temp
            
    return phi

    
maxErr1 = 0.001

test = zeros(fSol.shape)
GS_tested = GaussSeidel(test,fSol,5)
Jac_tested = Jacobi(test,fSol,5)


t0 = time.clock() # test time of the algorithm
Jac_testedWhile = JacobiWhile(test,fSol,maxErr1,1.9)
GS_testedWhile = GaussSeidelWhile(test,fSol,maxErr1,1.9)
deltaT = time.clock()-t0
print(deltaT)
#err = compError(Jac_tested)

#3D figure
fig = figure(figsize=(4,3))
axes = fig.add_subplot(111,projection='3d')
#plots
pl = axes.plot_surface(X,Y,GS_tested,rstride=1,cstride=1,cmap=cm.jet)
cb = fig.colorbar(pl)

#3D figure
fig2 = figure(figsize=(4,3))
axes = fig2.add_subplot(111,projection='3d')
#plots
pl2 = axes.plot_surface(X,Y,Jac_tested,rstride=1,cstride=1,cmap=cm.jet)
cb2 = fig.colorbar(pl2)

#3D figure
fig2 = figure(figsize=(4,3))
axes = fig2.add_subplot(111,projection='3d')
#plots
pl2 = axes.plot_surface(X,Y,Jac_testedWhile,rstride=1,cstride=1,cmap=cm.jet)
cb2 = fig.colorbar(pl2)
