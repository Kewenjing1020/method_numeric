import numpy as np
from pylab import *

###########################################################################
# Forward Euler method for the system with right hand side function Frhs 
###########################################################################
def EulerForward_OsciHarm(y0,t0,tEnd,dt,w):
	# Number of points in the time mesh
	N = int((tEnd-t0)/dt)+1

	# Time array
	time = np.linspace(t0,tEnd,N)

	if (type(y0) is float):
		# Scalar ODE
		sol = np.zeros(N)
		sol[0] = y0
	else:
		# Vectorial ODE
		sol = np.zeros((N,np.size(y0)))
		sol[0] = y0
      
	# Resolution of the system
	for i in range(1,N):
		sol[i,0] = sol[i-1,0]+dt*sol[i-1,1]
		sol[i,1] = sol[i-1,1]-dt*w**2*sol[i-1,0]
	
	return time,sol



###########################################################################
# Backward Euler method for the system with right hand side function Frhs
###########################################################################
def EulerBackward_OsciHarm(y0,t0,tEnd,dt,w):
	# Number of points in the time mesh
	N = int((tEnd-t0)/dt)+1

	# Time array
	time = np.linspace(t0,tEnd,N)

	if (type(y0) is float):
		# Scalar ODE
		sol = np.zeros(N)
		sol[0] = y0
	else:
		# Vectorial ODE
		sol = np.zeros((N,np.size(y0)))
		sol[0] = y0

	# Resolution of the system
	for i in range(1,N):
          sol[i,0]=(sol[i-1,0]+dt*sol[i-1,1])/(1.+dt**2*w**2)
          sol[i,1]=sol[i-1,1]-dt*w**2*sol[i,0]

	return time,sol



###########################################################################
# Trapezoidal method for the system with right hand side function Frhs
###########################################################################
def Trapezoidal_OsciHarm(y0,t0,tEnd,dt,w):
	# Number of points in the time mesh
	N = int((tEnd-t0)/dt)+1

	# Time array
	time = np.linspace(t0,tEnd,N)

	if (type(y0) is float):
		# Scalar ODE
		sol = np.zeros(N)
		sol[0] = y0
	else:
		# Vectorial ODE
		sol = np.zeros((N,np.size(y0)))
		sol[0] = y0

	# Resolution of the system
	for i in range(1,N):
          sol[i,0]=(sol[i-1,0]+dt*sol[i-1,1]-dt**2*w**2/4.*sol[i-1,0])/(1.+dt**2*w**2/4.)
          sol[i,1]=sol[i-1,1]-dt*w**2*sol[i,0]/2.-dt*w**2*sol[i-1,0]/2.

	return time,sol
