#########################################################
# Main program for the part 2 : Evolution of a population
#########################################################

from Population import *
from NumericalScheme import *
from pylab import *

# close all figures
close('all')

# Initial population (between 0 and 1)
p0 = 0.5

# Start and final time
t0 = 0.0
tEnd = 20.0

# Analytical solution
timeAnal,solAnal = analSolPop(p0,t0,tEnd)



# Time steps to test
dtTab = array([3.0,1.0,0.1,0.01,0.001])

# Numerical solution for different time step
error = zeros(size(dtTab))
for i in xrange(size(dtTab)):
	# Computation of the solution
	timeNum, solNum = EulerForward(p0,t0,tEnd,dtTab[i])

	# Computation of the error
	error[i] = abs(solAnal[-1]-solNum[-1])

# Plot of the error
fig=plt.figure()
plot(dtTab,error,label='Forward')
xlabel(r'$\Delta t$',fontsize=20)
ylabel('Error',fontsize=20)
xscale('log')
yscale('log')
legend()
show()
if not os.path.exists('./Figures'):
  os.makedirs('./Figures')
fig.savefig("./Figures/Population_SolVStime_ForwardEuler.pdf", format="pdf")



# Time evolution with different time step
for i in xrange(size(dtTab)):
	# Solution with the different scheme
	timeForward,solForward = EulerForward(p0,t0,tEnd,dtTab[i])
	timeBackward,solBackward = EulerBackward(p0,t0,tEnd,dtTab[i])
	timeTrapez,solTrapez = Trapezoidal(p0,t0,tEnd,dtTab[i])
	timeBackwardLin,solBackwardLin = EulerBackwardLin(p0,t0,tEnd,dtTab[i])
	timeTrapezLin,solTrapezLin = EulerTrapezoidalLin(p0,t0,tEnd,dtTab[i])

	# Plot of the different solution
	fig=plt.figure()
	plot(timeAnal,solAnal,'k',label='Analytic')
	plot(timeForward,solForward,label='Euler Forward')
	plot(timeBackward,solBackward,label='Euler Backward')
	plot(timeTrapez,solTrapez,label='Trapezoidal')
	plot(timeBackwardLin,solBackwardLin,label='Euler Backward Linearized')
	plot(timeTrapezLin,solTrapezLin,label='Trapezoidal Linearized')
	xlabel(r'$t$',fontsize=20)
	ylabel(r'$P$',fontsize=20)
	ylim((-0.5,0.5))
	title(r'$\Delta t = '+str(dtTab[i])+'$',fontsize=30)
	legend()
	show()
	fig.savefig("./Figures/Population_SolVStime_"+str(dtTab[i])+"_all.pdf", format="pdf")
	
