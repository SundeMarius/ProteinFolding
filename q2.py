import numpy as np
import Protein as prot
import matplotlib.pyplot as plt

# Define some things for plotting
font = {'family': 'normal', 'weight': 'bold', 'size': 16}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\boldmath']

Trange = np.array([10e-12, 1500])  # Range of possible temperatures
acc = 2  # Number of plotting points


def quest_2_1(acc):
	# Parameters for protein of length 15
	dmax = 13000
	s = 5e-4

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
	plt.plot(Tvalues, Evalues, lw=3, label=r"15 monomers", color="crimson")
	plt.xlabel(r"$T$  [K]", size=20)
	plt.ylabel(r"$\langle E \rangle$  [J]", size=20)
	plt.legend(loc="best")
	plt.grid()
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
	plt.plot(twistValues, Evalues, lw=3, label=r"$T$ = 0K", color="crimson")
	plt.xlabel(r"\textbf{Twists}", size=20)
	plt.ylabel(r"$E$  [J]", size=20)
	plt.legend(loc="best")
	plt.grid()
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
	plt.plot(twistValues, Evalues, lw=3, label=r"$T$ = 500K", color="crimson")
	plt.xlabel(r"\textbf{Twists}", size=20)
	plt.ylabel(r"$E$  [J]", size=20)
	plt.legend(loc="best")
	plt.grid()
	plt.show()


def quest_2_5(acc):
	# Parameters for protein of length 30
	dmax = 26000
	s = 0.001

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
	plt.plot(Tvalues, Evalues, lw=3, label=r"30 monomers", color="crimson")
	plt.xlabel(r"Temperature $T$  [K]", size=20)
	plt.ylabel(r"Mean energy $\langle E \rangle$  [J]", size=20)
	plt.grid()
	plt.show()


# Comment out the ones you don't want to run

quest_2_1(acc)
#quest_2_2(acc)
#quest_2_5(acc)
