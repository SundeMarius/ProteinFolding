import numpy as np
import Protein as prot
from matplotlib import pyplot as plt

# Define some things for plotting
font = {'family': 'normal', 'weight': 'bold', 'size': 16}
plt.rc('font', **font)

Trange = np.array([10e-12, 1500])  # Range of possible temperatures
acc = 20  # Number of plotting points


def quest_2_1(acc):
	# Parameters for protein of length 15
	dmax = 10000
	s = 0.002

	# x-axis of plot
	Tvalues = np.linspace(Trange[0], Trange[1], acc)
	# y-axis of plot
	Evalues = np.zeros(acc)

	for i in range(acc):
		# Find number of twists necessary
		D = int(dmax * np.exp(-s * Tvalues[i]))
		# Define a new protein for each temperature
		P = prot.Protein(15)

		energies = np.zeros(D)
		print("\n\nTemperature:", Tvalues[i], "|| # Twists:", D, "|| Iter:", i + 1, "of", acc, "\n\n")

		for j in range(D):
			# Perform a single twist with energy and thermal fluctuation conserns
			P = prot.randomTwist(P, Tvalues[i])
			# Add energy to list
			energies[j] = P.E
		# Find the average given the temperature
		Evalues[i] = np.average(energies)

	# Plot
	plt.plot(Tvalues, Evalues)
	plt.xlabel(r"T")
	plt.ylabel(r"[E]")
	plt.show()


def quest_2_2(acc):
	# x-axis of plot
	twistValues = np.linspace(0, 5000, acc)
	# y-axis of plot
	Evalues = np.zeros(acc)

	dist = int(twistValues[1] - twistValues[0])


	### T = 0 ###


	# Define a protein
	P = prot.Protein(15)

	for i in range(acc):
		for j in range(dist):
			P = prot.randomTwist(P, 10e-12)  # Cannot make temperature == 0 because of 0 division error
		Evalues[i] = P.E

		print("Iter: ", i+1, "/", acc)

	# Plot
	plt.plot(twistValues, Evalues, label=r"T = 0K")
	plt.legend()
	plt.xlabel(r"Twists")
	plt.ylabel(r"E")
	plt.show()


	### T = 500 ###


	# Define another protein
	P = prot.Protein(15)

	for i in range(acc):
		for j in range(dist):
			P = prot.randomTwist(P, 500)
		Evalues[i] = P.E

		print("Iter: ", i+1, "/", acc)

	# Plot
	plt.plot(twistValues, Evalues, label=r"T = 500K")
	plt.legend()
	plt.xlabel(r"Twists")
	plt.ylabel(r"E")
	plt.show()


def quest_2_5(acc):
	# Parameters for protein of length 30
	dmax = 10000	### DON'T KNOW THESE YET
	s = 0.002

	# x-axis of plot
	Tvalues = np.linspace(Trange[0], Trange[1], acc)
	# y-axis of plot
	Evalues = np.zeros(acc)

	for i in range(acc):
		# Find number of twists necessary
		D = int(dmax * np.exp(-s * Tvalues[i]))
		# Define a new protein for each temperature
		P = prot.Protein(30)

		energies = np.zeros(D)
		print("\n\nTemperature:", Tvalues[i], "|| # Twists:", D, "|| Iter:", i + 1, "of", acc, "\n\n")

		for j in range(D):
			# Perform a single twist with energy and thermal fluctuation conserns
			P = prot.randomTwist(P, Tvalues[i])
			# Add energy to list
			energies[j] = P.E
		# Find the average given the temperature
		Evalues[i] = np.average(energies)

	# Plot
	plt.plot(Tvalues, Evalues)
	plt.xlabel(r"T")
	plt.ylabel(r"[E]")
	plt.show()


# Comment out the ones you don't want to run

#quest_2_1(acc)
#quest_2_2(acc)
quest_2_5(acc)
