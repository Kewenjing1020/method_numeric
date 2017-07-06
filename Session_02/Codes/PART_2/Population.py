
import numpy as np
from math import exp

################################################################################
# Right hand side function for the population evolution system
################################################################################
def FrhsPop(t,p):
	return p*(p-1.0)
	


################################################################################
# Jacobian of the right hand side function for the population evolution system
################################################################################
def JacFrhsPop(t,p):
	return 2.0*p-1.0
	


################################################################################
# Analytical solution for the population evolution system 
################################################################################
def analSolPop(p0,t0,tEnd):
	# Parameter A
	A = 1.0-1.0/p0

	# Number of points in the solution (step of 0.01)
	N = 1+int(100*(tEnd-t0))

	# Time array
	time = np.linspace(t0,tEnd,N)

	# Solution array
	sol = np.zeros(N)

	# Fill the solution array
	for i in xrange(N):
		sol[i] = 1.0/(1.0-A*exp(time[i]))

	return time,sol
