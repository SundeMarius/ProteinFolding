import numpy as np
import random
import Protein as prot
from matplotlib import pyplot as plt

kb = 1.38e-23  # J/K (Boltzmann's constant)

# Define some things for plotting
font = {'family': 'normal', 'weight': 'bold', 'size': 16}
plt.rc('font', **font)

def meanEnergy(protein,T):
	"""
	:param protein: The polymer to twist
	:param T: temperature-interval
	:return: non (creates a plot that solves exercise 2.1)
	"""
	def d(dmax,s,T):
		"""
		:param dmax: d at T = 0
		:param s: decreasing-rate (choose a staisfactory value, s>0)
		:param T: temperature in K
		:return: the number of necessary iterations d
		"""
		return int(dmax*np.exp(-s*T))

	E = np.zeros(len(T))

	for k in range(len(T)):
		# Update temperature and boltzmann's paramter beta
		t = T[k]
		B = 1/(kb*t)

		Eprob, Z = (0,0) # Zero out mean energy and Partitioning-sum (eq. 4)

		#Create a copy-protein
		prot2 = prot.Protein(protein.n)
		Ems1 = 0

		# Twist protein d times (d+1 microstates), calculate mean energy
		for i in range(d(11000,0.006,t)):

			#Choose rotation-way and monom-pivot randomly
			rotate = prot.randomBool()
			pivot = prot2.randomMonomer()

			#twist the copy-protein
			twisted = prot2.twist(pivot,rotate)

			if twisted:

				#Calculate energy in copy-protein
				Ems2 = prot2.calculateEnergy() #Energy in current microstate
				#Check if new energy-state is less than previous, update protein if it is (energy-minimization)
				if Ems2 < Ems1:
					protein = prot2 #Update protein
					Ems1 = Ems2 #Update energy of protein
				#Add thermal fluctuations (see worksheet)
				elif random.randint(0,1) < np.exp(-B*(Ems2-Ems1)):
					protein = prot2 #Update protein
					Ems1=Ems2 #Update energy of protein
			else:
				#Try again until twisting is possible
				continue

			#Calculate contribution to mean energy
			Ems = Ems2
			Boltz = np.exp(Ems * B)  # Boltzmann-factor in current microstate
			Eprob += Ems * Boltz
			Z += Boltz

		#Add result to energy-array
		if Z > 0:
			E[k] = (1/Z)*Eprob
			print(E[k])

	#Plot results
	plt.xlabel(r"Temperature $T[K]$")
	plt.ylabel(r"Mean energy $E$ of polymer")
	plt.grid(), plt.legend(loc="best")
	plt.plot(T,E,color="r",label="$E(T)$")
	plt.show()

def BindingEnergy():

def main():
	#Create a protein with length 15, and misc. variables/intervals
	protein = prot.Protein(15)
	Nt = 700 #number of steps in temp-interval
	To = 1e-2
	Tf = 1500
	T = np.linspace(To,Tf,Nt+1)

	#Calculate mean energy of protein
	meanEnergy(protein,T)

main()