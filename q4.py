import Protein as prot
import numpy as np
import matplotlib.pyplot as plt

# Define some things for plotting
font = {'family': 'normal', 'weight': 'bold', 'size': 16}
plt.rc('font', **font)

numTwists = 600  # The number of twists at each temperature
Tincrement = 30  # The distance between each temperature
protLength = 15 # Length of protein
acc = 30  # Must not exceed 30000

def quest_4_1(acc,protLength):
	# Define a protein
	P = prot.Protein(protLength)
	# x-axis of plot
	twistValues = np.linspace(0, 30000, acc)
	# y-axis of plot
	Evalues = np.zeros(acc)

	dist = int(twistValues[1] - twistValues[0])

	T = 1500 + 10e-12 # Initializing temperature
	counter = 0  # Initializing the increment between each temperature level

	# Plot while temperature is positive
	while T > 0:
		# For every plotting point
		for i in range(acc):
			# Inbetween plotting points
			for j in range(dist):
				# Check if temperature needs to be updated
				if counter == numTwists:
					T -= Tincrement
					counter = 0
				P = prot.randomTwist(P, T)
				counter += 1
			Evalues[i] = P.E

			print("Iter: ", i + 1, "/", acc)

	# Plot
	plt.plot(twistValues, Evalues)
	plt.xlabel(r"Twists")
	plt.ylabel(r"E")
	plt.show()

def quest_4_2(acc,protLength):
	#Define protein
	protein = prot.Protein(protLength)

	#x and y-values to plot (temperature and avg. energy)
	temp = np.linspace(1500,1e-12,acc)
	Evalues = np.zeros(acc)

	for i in range(acc):
		#perform 600 twists at every temperature, record avg-energy at every temp
		T = temp[i]
		energy = np.zeros(numTwists)
		for j in range(numTwists):
			#Twisting
			protein = prot.randomTwist(protein,T)
			#record energy in current micro-state
			energy[j] = protein.E
		#Add average energy in this temperature to plotting-vector
		Evalues[i] = np.average(energy)
		#Print current state:
		print("Iter: ",i+1,"/",acc, "\tTemp: ",T)

	plt.plot(temp, Evalues, lw='2', c='g')
	plt.xlabel(r"Temperature $T$ [K]")
	plt.ylabel(r"Mean energy $\langle E \rangle$  [J]")
	plt.legend(loc="best")
	plt.grid()
	plt.show()

def quest_4_3(acc,protLength):
	# Define protein
	protein = prot.Protein(protLength)

	# x and y-values to plot (temperature and avg. proteinlength)
	temp = np.linspace(1500, 1e-12, acc)
	Lvalues = np.zeros(acc)

	for i in range(acc):
		# perform 600 twists at every temperature, record avg-energy at every temp
		T = temp[i]
		lengths = np.zeros(numTwists)
		for j in range(numTwists):
			# Twisting
			protein = prot.randomTwist(protein, T)
			# record max-diameter in current micro-state
			lengths[j] = protein.D
		# Add average max-diameter in this temperature to plotting-vector
		Lvalues[i] = np.average(lengths)
		# Print current state:
		print("Iter: ", i + 1, "/", acc, "\tTemp: ", T)

	plt.plot(temp, Lvalues, lw='2', c='g')
	plt.xlabel(r"Temperature $T$ [K]")
	plt.ylabel(r"Mean protein-length $\langle L \rangle$  [1]")
	plt.legend(loc="best")
	plt.grid()
	plt.show()

def quest_4_4(acc,protLength,numbOfCoolDowns):
	#Do several cool-downs on a protein, and draw the protein

	# Define protein
	protein = prot.Protein(protLength)

	# Temperature-range
	temp = np.linspace(1500, 1e-12, acc)

	for k in range(numbOfCoolDowns):
		for i in range(acc):
			# perform 600 twists at every temperature, record avg-energy at every temp
			T = temp[i]
			lengths = np.zeros(numTwists)
			for j in range(numTwists):
				# Twisting
				protein = prot.randomTwist(protein, T)
			# Print current state:
			print("Iter: ", i + 1, "/", acc, "\tTemp: ", T)
	protein.draw()

#quest_4_1(acc,protLength)
#quest_4_2(acc,protLength)
#quest_4_3(acc,protLength)
#quest_4_4(acc,protLength,3)