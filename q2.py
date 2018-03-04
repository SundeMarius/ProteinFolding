import numpy as np
import random
import Protein as prot
from matplotlib import pyplot as plt

# Define some things for plotting
font = {'family': 'normal', 'weight': 'bold', 'size': 16}
plt.rc('font', **font)

def exercise2(protein,T):
	"""
	:param protein: The polymer to twist
	:param T: temperature-interval
	:return: non (creates a plot that solves exercise 2.1)
	"""
	kb = 1.38e-23  # J/K (Boltzmann's constant)
	Z, Ems, Boltz, Eprob = (0,0,0,0) #Partition-sum, Energy in current microstate, boltzmann-factor, weighted energy (eq. 4)
	def d(dmax,s,T):
		"""
		:param dmax: d at T = 0
		:param s: decreasing-rate (choose a staisfactory value, s>0)
		:param T: temperature in K
		:return: the number of necessary iterations d
		"""
		return int(dmax*np.exp(-s*T))

	E = np.zeros(len(T))
	t = T[0]
	it = 1
	while it < len(T):

		#Twist protein d times (d+1 microstates), calculate Z
		B = 1/(kb*t)

		for i in range(d(11000,0.5,t)):
			#Choose rotation-way and monom-pivot randomly
			rotate = random.randint(0, 1)
			pivot = random.randint(2, 14)

			if protein.isLegalTwist(pivot,rotate):
				protein.twist(pivot,rotate)
				protein.draw()
				Ems = protein.calculateEnergy() #Energy in current microstate
				Boltz = np.exp(-Ems*B) #Boltzmann-factor in current microstate
				Eprob += Ems*Boltz
				Z += Boltz
			else:
				continue
		#Calculate mean potential-energy, add to plot-array
		#E[it] = (1/Z)*Eprob
		#Update iterator and temperature
		t = T[it]
		it += 1

	plt.xlabel(r"Temperature $T[K]$")
	plt.ylabel(r"Mean energy $E$ of polymer")
	plt.grid(), plt.legend(loc="best")
	plt.plot(T,E,color="r",label="$E(T)$")
	plt.show()


def main():
	#Create a protein with length 15, and misc. variables/intervals
	protein = prot.Protein(15)
	T = np.linspace(1e-12,1500,1500)

	exercise2(protein,T)

main()