import numpy as np
import scipy.spatial.distance as ssd
import matplotlib.pyplot as plt
import Protein as prot

# Define constants:
# s <= 0.01 for sufficiently convergent results for high temperatures
dmax = 11000
s = 0.001
acc = 10

def diameter(protein):
	"""

	:param protein: Protein class member; an instance of protein
	:return: Float; the diameter of the protein, defined as the longest distance between two monomers
	"""

	# Add all non-zero elements as points to an array
	points = np.zeros([protein.n, 2])
	currentCoords = np.array(protein.G.findElement(1))
	for i in range(protein.n):
		points[i] = currentCoords
		currentCoords = protein.G.searchAdjacent(currentCoords, i+2)

	# Do some stackoverflow magic
	D = ssd.pdist(points)
	D = ssd.squareform(D)
	N, [I_row, I_col] = np.nanmax(D), np.unravel_index(np.argmax(D), D.shape)

	return N

def d(dmax, s, T):
	"""

	:param dmax: Float; d at T = 0 (initial value)
	:param s: Float; Decreasing-rate (choose a staisfactory value, s>0)
	:param T: Float; Temperature in K
	:return: Int; The number of necessary iterations d
	"""
	return int(dmax * np.exp(-s * T))

T = np.array([1, 1500])
Tvalues = np.linspace(T[0], T[1], acc)
diametervalues = np.zeros(acc)

# For each temperature
for i in range(acc):
	# Find number of twists necessary
	D = d(dmax, s, Tvalues[i])
	diameters = np.zeros(D)
	# Define a new protein for each temperature
	P = prot.Protein(15)
	print("\n\nTemperature:" , Tvalues[i], "|| # Twists:" , D, "|| Iter:", i+1, "of" , acc, "\n\n")
	# Perform twists and find the diameter after each twist
	for j in range(D):
		P.randomTwist()
		diameters[j] = diameter(P)
		#print(j+1)
	# Find the average diameter for the given temperature
	diametervalues[i] = np.average(diameters)

# Plot
plt.plot(Tvalues, diametervalues)
plt.show()
