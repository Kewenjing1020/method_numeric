
from math import *
import numpy as np

######################################################################
# Forward Euler method for the specific problem of the spring mass
######################################################################
def EulerForward(y0,t0,tEnd,deltat,w):
	# Value of w is an argument of the function
	wSquare = w**2

	# Number of points in the time mesh
	N = int((tEnd-t0)/deltat)+1

	# Time array
	time = np.linspace(t0,tEnd,N)

	# Solution arrays
	sol = np.zeros((N,2))
	sol[0] = y0

	# Resolution of the system
	for i in range(1,N):
		sol[i,0] = sol[i-1,0]+deltat*sol[i-1,1]
		sol[i,1] = sol[i-1,1]-deltat*wSquare*sol[i-1,0]
	
	return time,sol

