import Protein as prot
import numpy as np
import matplotlib.pyplot as plt

numTwists = 600  # The number of twists at each temperature
Tincrement = 30  # The distance between each temperature
acc = 40  # Must not exceed 30000

def quest_4_1(acc):
	# Define a protein
	P = prot.Protein(15)
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


quest_4_1(acc)
