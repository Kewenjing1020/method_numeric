# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 15:49:14 2016

@author: kewenjing
"""
from pylab import*

def forward_euler(y0,t_ini,t_end,t_step,w):
    n=(int)((t_end-t_ini)/t_step);
    time=linspace(t_ini,t_end,n);
    y=np.zeros((n+2,2));    
    #y=array[n+2];
    y[0][0]=y0[0];
    y[0][1]=y0[1];
    for i in range(n):
        y[i+1][0]=y[i][0]+t_step*y[i][1];
        y[i+1][1]=y[i][1]-t_step*w**2*y[i][0];
    return time,y;
        
    
y0=np.zeros(2);
y0[0]=1;
y0[1]=0;
(time,y)=forward_euler(y0,0,10,0.001,4);
plot(y[][0],time)
#print y
#plot(time,y[0]);