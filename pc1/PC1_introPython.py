# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 15:00:34 2016

@author: kewenjing
PC: introduction to  python
"""

from pylab import *
#ex1.1
def func1(t):
    return x**3+x**2+x+1
 
N=10   
x=linspace(-1,1,N)
y=func1(x)
xlabel('X')
ylabel('Y')
title(r"$f(x)=x^3+x^2+x+1$")
pl=plot (x,y,'r--')
setp(pl,'linewidth',3)
setp(pl,'color','r')

#ex1.2
def func2(m,n):
    return m**2+n**2-2*m*n
    
x,y=mgrid[-1:1:0.1,-1:1:0.1]
g=func2(x,y)
fig,axes=subplots(1,2,figsize=(8,4)) #plots
axes[0].pcolor(x,y,g) 
axes[0].set_title('pcolor') 
axes[1].contour(x,y,g) 
axes[1].set_title('contour')


#figure
fig2=figure(figsize=(4,3))
#axes 
axes=fig2.add_subplot(1,1,1,projection='3d')
 #plots 
pl2=axes.plot_surface(x,y,g,cmap=cm.coolwarm) 
cb=fig2.colorbar(pl)