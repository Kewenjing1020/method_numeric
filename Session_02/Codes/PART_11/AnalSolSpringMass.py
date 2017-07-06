
#from math import sqrt, atan, pi, cos
from numpy import *
from math import *


######################################################################
# Analytical solution for the spring mass system 
######################################################################
def analSol(y0,t0,tEnd,w):
	# Initial position and velocity
	x0 = y0[0]
	u0 = y0[1]

	# Parameters A and phi
	A = sqrt(x0**2+(u0/w)**2)
	if (x0 > 0):
		phi = arctan(-u0/(x0*w))
	else:
		# Initial position equal to 0 : value of phi is +/-pi/2 
             #depending on the sign of u0
		if (u0 > 0):
			phi = -pi/2.0
		else:
			phi = pi/2.0

	# Time array where the solution will be defined 
    (must have enough point to capture the oscillation)
	T = (2.0*pi)/w             # Period of the solution
	N = 20*int((tEnd-t0)/T)+1    # 20 point per period
	time = linspace(t0,tEnd,N)

	# Solution arrays
	sol = zeros(N)

	# Fill the solution array
	for i in range(N):
		sol[i] = A*cos(w*time[i]+phi)		

	return time,sol
