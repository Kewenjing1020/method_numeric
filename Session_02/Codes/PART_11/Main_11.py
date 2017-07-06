#########################################################
# Main program for the part 1.1 : Explicit method
#########################################################

from AnalSolSpringMass import *
from EulerForwardSpringMass import *
from pylab import *
import os

#close all figures
close('all')

# Initial condition, and value of w
y0 = zeros(2)
y0[0] = 1.0  # Initial position
y0[1] = 0.0  # Initial velocity
w = 4.0

# Start and final time
t0 = 0.0
tEnd = 10.0

# Analytical solution
timeAnal,solAnal = analSol(y0,t0,tEnd,w)



# Plot of analytical solution
fig=plt.figure()
plot(timeAnal,solAnal,'k',label=r'$Analytic$')

# Time steps to test
dtTab = array([0.01,0.001,0.0001,0.00001])

# Numerical solution for different time step
error = zeros(size(dtTab))



for i in range(size(dtTab)):
	# Computation of the solution
	timeNum, solNum = EulerForward(y0,t0,tEnd,dtTab[i],w)

	# Computation of the error
	error[i] = abs(solAnal[-1]-solNum[-1,0])

	# Plot of numerical solution (only position)
	plot(timeNum,solNum[:,0],label=r'$\Delta t='+str(dtTab[i])+'$')

# Label and show the plot
xlabel(r'$t [s]$',fontsize=20)
ylabel(r'$x [m]$',fontsize=20)
ylim((-3,3))
legend(bbox_to_anchor=(0.35, 1.2))
show()

#create the Figure file
if not os.path.exists('./Figures'):
    os.makedirs('./Figures')
fig.savefig("./Figures/HarmOsci_SolVStime_ForwardEuler.pdf", format="pdf")


# Plot of the error
fig=plt.figure()
plot(dtTab,error)
xlabel(r'$\Delta t [s]$',fontsize=20)
ylabel('Error [m]',fontsize=20)
xscale('log')
yscale('log')
show()
fig.savefig("./Figures/HarmOsci_Error_ForwardEuler.pdf", format="pdf")
