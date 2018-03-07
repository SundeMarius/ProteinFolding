import numpy as np
import Protein as prot
import matplotlib.pyplot as plt

# Define some things for plotting
font = {'family': 'normal', 'weight': 'bold', 'size': 16}
plt.rc('font', **font)

Trange = np.array([10e-12, 1500])  # Range of possible temperatures
acc = 20


def quest_3_short(acc):
	# Parameters for protein of length 15
	dmax = 10000
	s = 0.002

	# x-axis of plot
	Tvalues = np.linspace(Trange[0], Trange[1], acc)
	# y-axis of plot
	Dvalues = np.zeros(acc)

	# For each temperature
	for i in range(acc):
		# Find number of twists necessary
		D = int(dmax * np.exp(-s * Tvalues[i]))
		diameters = np.zeros(D)

		# Define a new protein for each temperature
		P = prot.Protein(15)
		print("\n\nTemperature:", Tvalues[i], "|| # Twists:", D, "|| Iter:", i + 1, "of", acc, "\n\n")

		# Perform twists and find the diameter after each twist
		for j in range(D):
			# Perform a single twist with energy and thermal fluctuation conserns
			P = prot.randomTwist(P, Tvalues[i])
			# Add to diameter list
			diameters[j] = P.D
		# Find the average diameter for the given temperature
		Dvalues[i] = np.average(diameters)

	# Plot
	plt.plot(Tvalues, Dvalues)
	plt.xlabel(r"T")
	plt.ylabel(r"L")
	plt.show()


def quest_3_long(acc):
	# Parameters for protein of length 30
	dmax = 10000  ### DON'T ACTUALLY KNOW THESE YET
	s = 0.002

	# x-axis of plot
	Tvalues = np.linspace(Trange[0], Trange[1], acc)
	# y-axis of plot
	Dvalues = np.zeros(acc)

	# For each temperature
	for i in range(acc):
		# Find number of twists necessary
		D = int(dmax * np.exp(-s * Tvalues[i]))
		diameters = np.zeros(D)

		# Define a new protein for each temperature
		P = prot.Protein(30)
		print("\n\nTemperature:", Tvalues[i], "|| # Twists:", D, "|| Iter:", i + 1, "of", acc, "\n\n")

		# Perform twists and find the diameter after each twist
		for j in range(D):
			# Perform a single twist with energy and thermal fluctuation conserns
			P = prot.randomTwist(P, Tvalues[i])
			# Add to diameter list
			diameters[j] = P.D
		# Find the average diameter for the given temperature
		Dvalues[i] = np.average(diameters)

	# Plot
	plt.plot(Tvalues, Dvalues)
	plt.xlabel(r"T")
	plt.ylabel(r"L")
	plt.show()


# Comment out the ones you don't want to run

quest_3_short(acc)
#quest_3_long(acc)
