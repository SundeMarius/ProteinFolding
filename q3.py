import numpy as np
import scipy.spatial.distance as ssd

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

