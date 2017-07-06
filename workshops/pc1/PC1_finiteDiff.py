# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 15:57:26 2016

@author: kewenjing
PC1: Finite Differences & ODE
"""

from pylab import *

#ex1.1
def f(t):
    return sin(t)/(t**3)
    
def df_num(t,h):
    return (f(t+h)-f(t))/h


def df_exacte(t):
    return (t*cos(t)-3*sin(t))/(t**4)
    
def d2f_exacte(x):
    return (-(x**2)*sin(x)-6*x*cos(x)+12*sin(x))/(x**5)
    
def error_1st(t,h):
    return abs(df_num(t,h) - df_exacte(t))
    
def df_num_2nd(t,h):
    return (f(t+h)-f(t-h))/(2*h)
  
def error_2nd(x,h):
    return abs(ddf_num_2nd(x,h)-df_exacte(x))
    
def df_num_4th(x,h):
    return (f(x-2*h)-8*f(x-h))
    
    
h=linspace(10**-4,1)
err1=error_1st(4,h)
err2=error_2nd(4,h)
plot(h,err1)
plot(h,err2)
xscale('log')
yscale('log')


#ex1.2
def dev_2nd_2h(x,h):
    return (f(x+h)-2*f(x)+f(x-h))/(h**2)

def dev_2nd_4h(x,h):
    return (f(x+h)-2*f(x)+f(x-2*h))/(4*h**4)
    
def dev_4th(x,h):
    return (f(x-2*h)+16*f(x-h)-30*f(x)+16*f(x+h)-f(x+2*h))
    
fx=df_exacte(4)
plot(h, abs(dev_2nd_2h(4,h)-))






#ex2

































