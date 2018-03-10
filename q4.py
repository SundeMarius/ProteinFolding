import Protein as prot
import numpy as np
import matplotlib.pyplot as plt

# Define some things for plotting
font = {'family': 'normal', 'weight': 'bold', 'size': 16}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\boldmath']

numTwists = 600  # The number of twists at each temperature
Tincrement = 30  # The distance between each temperature
protLength = 15 # Length of protein
acc = 1500  # Must not exceed 30000

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

			print("Iter: ",i + 1,"/",acc)

	# Plot
	plt.plot(twistValues, Evalues, lw=1, label=r"15 monomers", color="crimson")

	changeTempPoints = np.linspace(0, 30000, int(1500/Tincrement))
	for point in changeTempPoints:
		plt.axvline(x=point)

	plt.xlabel(r"\textbf{Twists}", size=20)
	plt.ylabel(r"$E$  [J]", size=20)
	plt.legend(loc="best")
	plt.grid()
	plt.show()

def quest_4_2(acc,protLength):
	#Define protein
	protein = prot.Protein(protLength)

	#x and y-values to plot (temperature and avg. energy)
	temp = np.linspace(1500,1e-12,acc)
	Evalues = np.zeros(len(temp))

	for i in range(len(temp)):
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

	plt.plot(temp, Evalues, lw=0.5, label=r"15 monomers", color="crimson")
	plt.xlabel(r"$T$  [K]", size=20)
	plt.xlim(1500,0)
	plt.ylabel(r"$\langle E \rangle$  [J]", size=20)
	plt.legend(loc="best")
	plt.grid()
	plt.show()

def quest_4_3(acc,protLength):
	# Define protein
	protein = prot.Protein(protLength)

	# x and y-values to plot (temperature and avg. proteinlength)
	temp = np.linspace(1500, 1e-12, acc)
	Lvalues = np.zeros(len(temp))

	for i in range(len(temp)):
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
		print("Iter: ",i + 1,"/",acc,"\tTemp: ",T)

	plt.plot(temp, Lvalues, lw=0.6, label=r"%s monomers" %protein.n, color="crimson")
	plt.xlabel(r"$T$  [K]", size=20)
	plt.xlim(1500,0)
	plt.ylabel(r"$\langle L \rangle$  [1]", size=20)
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
		for i in range(len(temp)):
			# perform 600 twists at every temperature, record avg-energy at every temp
			T = temp[i]
			lengths = np.zeros(numTwists)
			for j in range(numTwists):
				# Twisting
				protein = prot.randomTwist(protein, T)
			# Print current state:
			print("Iter: ",i + 1,"/",len(temp),"\tTemp: ",T)
	protein.draw()

#quest_4_1(acc,protLength)
#quest_4_2(acc,protLength)
#quest_4_3(acc,protLength)
quest_4_4(acc,protLength,3)