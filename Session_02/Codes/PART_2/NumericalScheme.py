
import numpy as np
import numpy.linalg as alg
from scipy.optimize import fsolve


#####################################################################################
# Forward Euler method for the system with right hand side function Frhs 
#####################################################################################
def EulerForward(y0,t0,tEnd,dt):
	# Number of points in the time mesh
	N = int((tEnd-t0)/dt)+1

	# Time array
	time = np.linspace(t0,tEnd,N)

	sol = np.zeros(N)
	sol[0] = y0

	# Resolution 
	for i in xrange(1,N):
		sol[i] = sol[i-1]-dt*sol[i-1]*(1.-sol[i-1])
	
	return time,sol



#####################################################################################
# Backward Euler method for the system with right hand side function Frhs
#####################################################################################
def EulerBackward(y0,t0,tEnd,dt):
	# Number of points in the time mesh
	N = int((tEnd-t0)/dt)+1

	# Time array
	time = np.linspace(t0,tEnd,N)

	sol = np.zeros(N)
	sol[0] = y0

	# Resolution
	for i in xrange(1,N):
           delta= (1.+dt)**2-4.*dt*sol[i-1]
           sol[i]=1./2./dt*(1.+dt-delta**0.5)
	return time,sol



#####################################################################################
# Trapezoidal method for the system with right hand side function Frhs
#####################################################################################
def Trapezoidal(y0,t0,tEnd,dt):
	# Number of points in the time mesh
	N = int((tEnd-t0)/dt)+1

	# Time array
	time = np.linspace(t0,tEnd,N)

	sol = np.zeros(N)
	sol[0] = y0
 
	# Resolution 
	for i in xrange(1,N):
		#sol[i] = sol[i-1]-dt/2.*sol[i-1]*(1.-sol[i-1])-dt/2.*sol[i]*(1.-sol[i])
		#sol[i] = sol[i-1]-dt/2.*sol[i-1]*(1.-sol[i-1])-dt/2.*sol[i]*(1.-sol[i])
           delta= (1.+dt/2.)**2-4.*dt*(sol[i-1]-dt/2.*sol[i-1]*(1.-sol[i-1]))
           sol[i]=1./dt*(1.+dt/2.-delta**0.5)
	return time,sol



#####################################################################################
# Backward linearized method for the system with right hand side function Frhs
#####################################################################################
def EulerBackwardLin(y0,t0,tEnd,dt):
	# Number of points in the time mesh
	N = int((tEnd-t0)/dt)+1

	# Time array
	time = np.linspace(t0,tEnd,N)

	sol = np.zeros(N)
	sol[0] = y0

	# Resolution 
	for i in xrange(1,N):
		sol[i] = sol[i-1]-dt*sol[i-1]*(1.-sol[i-1])-dt*(1.-2*sol[i-1])*(-sol[i-1])
		sol[i] = sol[i]/(1.+dt*(1.-2*sol[i-1]))
	return time,sol



#####################################################################################
# Trapezoidal linearized method for the system with right hand side function Frhs
#####################################################################################
def EulerTrapezoidalLin(y0,t0,tEnd,dt):
	# Number of points in the time mesh
	N = int((tEnd-t0)/dt)+1

	# Time array
	time = np.linspace(t0,tEnd,N)

	sol = np.zeros(N)
	sol[0] = y0

	# Resolution 
	for i in xrange(1,N):
		sol[i] = sol[i-1]-dt*sol[i-1]*(1.-sol[i-1])-dt/2.*(1.-2*sol[i-1])*(-sol[i-1])
		sol[i] = sol[i]/(1.+dt/2.*(1.-2*sol[i-1]))

	return time,sol
