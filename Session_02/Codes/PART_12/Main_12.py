#########################################################
# Main program for the part 1.2 : Implicit method
#########################################################

from AnalSolSpringMass import *
from NumericalScheme import *
from pylab import *
from math import *

#close all figures
close('all')

# Initial solution
y0 = zeros(2)
y0[0] = 1.0
y0[1] = 0.0

# Pulsation
w = 4.0


# Start and final time
t0 = 0.0
tEnd = 10.0

# Analytical solution
timeAnal,solAnal = analSol(y0,t0,tEnd,w)



# Solution for a time step of 0.1 with the differents method
timeForward,solForward = EulerForward_OsciHarm(y0,t0,tEnd,0.1,w)
timeBackward,solBackward = EulerBackward_OsciHarm(y0,t0,tEnd,0.1,w)
timeTrapez,solTrapez = Trapezoidal_OsciHarm(y0,t0,tEnd,0.1,w)

# Plot of the different solution
fig=plt.figure()
plot(timeAnal,solAnal,'k',label='Analytic')
plot(timeForward,solForward[:,0],label='Forward Euler')
plot(timeBackward,solBackward[:,0],label='Backward Euler')
plot(timeTrapez,solTrapez[:,0],label='Trapezoidal')
ylim([-2.0,2.0])
xlabel(r'$t [s]$',fontsize=20)
ylabel(r'$x [m]$',fontsize=20)
legend()
show()
if not os.path.exists('./Figures'):
  os.makedirs('./Figures')
fig.savefig("./Figures/HarmOsci_SolVStime_all.pdf", format="pdf")




# Time steps to test
dtTab = array([0.1,0.01,0.001,0.0001])

# Numerical solution for different time step
errorForward = zeros(size(dtTab))
errorBackward = zeros(size(dtTab))
errorTrapez = zeros(size(dtTab))
for i in xrange(size(dtTab)):
	# Computation of the solution
	timeForward, solForward = EulerForward_OsciHarm(y0,t0,tEnd,dtTab[i],w)
	timeBackward, solBackward = EulerBackward_OsciHarm(y0,t0,tEnd,dtTab[i],w)
	timeTrapez, solTrapez = Trapezoidal_OsciHarm(y0,t0,tEnd,dtTab[i],w)

	# Computation of the error
	errorForward[i]  = abs(solAnal[-1]-solForward[-1,0])
	errorBackward[i] = abs(solAnal[-1]-solBackward[-1,0])
	errorTrapez[i]   = abs(solAnal[-1]-solTrapez[-1,0])

# Plot of the error
fig=plt.figure()
plot(dtTab,errorForward,label='Forward Euler')
plot(dtTab,errorBackward,label='Backward Euler')
plot(dtTab,errorTrapez,label='Trapezoidal')
xlabel(r'$\Delta t [s]$',fontsize=20)
ylabel('Error [m]',fontsize=20)
xscale('log')
yscale('log')
legend(bbox_to_anchor=(0.4, 1.0))
show()

fig.savefig("./Figures/HarmOsci_Error_All.pdf", format="pdf")