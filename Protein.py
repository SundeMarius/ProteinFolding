import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.colors as col

#Additional functions
def weakBinding():
	return random.uniform(-10.4e-21,-3.47e-21)

def kronDelta(i,j):
	"""
	:param i: Monom number i
	:param j: Monom number j
	:return: 0 if i and j is subsequential, 1 otherwise
	"""
	if (abs(i-j) <= 1):
		return 0
	else:
		return 1

def randomBool():
	return random.choice([True, False])


# Classes
class Grid:

	def __init__(self, gridSize):
		self.gSize = gridSize
		self.grid = np.zeros((gridSize, gridSize)).astype(np.int16)

	def draw(self):
		"""

		:return: None; prints a visualisation of the grid to the screen
		"""
		print(self.grid)

	def findElement(self, x):
		"""

		:param x: Int; Number of the monomer we want
		:return: np.array; The position of the given x ([i, j])
		"""

		a = np.argwhere(self.grid == x)
		return a[0]

	def searchAdjacent(self, pivotCoords, targetNr):
		"""

		:param pivotCoords: np.array; The coordinates of the point we are searching next to
		:param targetNr: Int; The number assigned to the point we are searching for
		:return: np.array; The coordinates of the target number
		"""

		# Check column:
		if (pivotCoords[0] != 0) and (self.grid[pivotCoords[0] - 1][pivotCoords[1]] == targetNr):
			return np.array([pivotCoords[0] - 1, pivotCoords[1]])

		if (pivotCoords[0] != self.gSize) and (self.grid[pivotCoords[0] + 1][pivotCoords[1]] == targetNr):
			return np.array([pivotCoords[0] + 1, pivotCoords[1]])

		# Check row:
		if (pivotCoords[1] != 0) and (self.grid[pivotCoords[0]][pivotCoords[1] - 1] == targetNr):
			return np.array([pivotCoords[0], pivotCoords[1] - 1])
		if (pivotCoords[1] != self.gSize) and (self.grid[pivotCoords[0]][pivotCoords[1] + 1] == targetNr):
			return np.array([pivotCoords[0], pivotCoords[1] + 1])

		# If target number has not been found:
#		print("\n### ERROR ###\nfunction searchAdjacent cannot find target number", targetNr, "next to initial number\n")
		# The previous line was removed because it falsely displayed an error in the diameter function
		return np.array([0, 0])

	def revealAdjacent(self, pivotCoords, side):
		"""

        :param pivotCoords: coordinates of the monom of interest
        :param side: 'Over' is y+1, 'Below' is y-1, 'Left' is x-1, 'Right' is x+1
        :return: the potential monom on the side of interest
        """
		x = pivotCoords[0]
		y = pivotCoords[1]

		if side == 'Over':
			return self.grid[x][y + 1]
		elif side == 'Below':
			return self.grid[x][y - 1]
		elif side == 'Left':
			return self.grid[x - 1][y]
		elif side == 'Right':
			return self.grid[x + 1][y]
		else:

			print("**ERROR** enter a valid side.")

	def present(self):
		"""

		:return: None; plots a pretty picture of the grid, suited to the protein
		"""
		colors = ['white', 'crimson']
		bounds = [0, 1, np.Inf]
		cmap = col.ListedColormap(colors)
		norm = col.BoundaryNorm(bounds, cmap.N)

		fig, ax = plt.subplots()
		ax.matshow(self.grid, cmap=cmap, norm=norm)

		for i in range(self.gSize):
			for j in range(self.gSize):
				c = self.grid[j, i]
				if c != 0:
					ax.text(i, j, str(c), va='center', ha='center', fontsize=20)

		plt.show()


class Protein:

	def __init__(self, length):
		self.n = length  # Length of protein
		self.N = length + 2  # Length of grid
		self.G = Grid(self.N)
		self.G.grid[int(np.round(self.N / 2)), int(np.round(self.N / 2 - self.n / 2)): int(np.round(self.N / 2 - \
			self.n / 2)) + self.n] = np.linspace(1, self.n, self.n)  # .astype(np.int16)
		self.midValue = int((length + 1) / 2)

	def isAboveMiddle(self, x):
		"""

		:param x: Int, the number assigned to the pivot point
		:return: Bool, True if x > middle value (average). Randomly True or False if x == middle. Else False
		"""
		if (x == self.midValue):
			return randomBool()
		else:
			return (x > self.midValue)

	def draw(self):
		"""

		:return: None; prints a visualisation of the grid to the screen
		"""
		self.G.draw()

	def present(self):
		"""

		:return: None; plots a pretty picture of the protein
		"""
		self.G.present()

	def calculateEnergy(self):
		"""

		:return: Float; Total energy E in the grid for a given microstate ms
		"""
		E = 0
		Sides = ['Over', 'Below', 'Left', 'Right']
		# Iterate over the whole protein
		for i in range(1, self.n + 1):
			# find monom i
			pos = self.G.findElement(i)
			# Search for extra neighbours, add energy-contribution if there is
			for k in Sides:
				j = self.G.revealAdjacent(pos, k)
				if j != 0:
					E += weakBinding() * kronDelta(i, j)
				else:
					continue
		# Return E when done.
		return E

	def twist(self, x, clockwise):
		"""

		:param x: int, the number of the pivot monomer
		:param clockwise: bool, True if the rotation is clockwise
		:return: bool, True if the rotation is legal. If true, the member grid is changed according to the twist
		"""

		# Find the coordinates of the pivot monomer
		pivotCoords = self.G.findElement(x)

		# Find out if you should iterate up or down
		itUp = self.isAboveMiddle(x)

		# Make a copy to avoid doing damage to the actual grid
		G = self.G.grid.__deepcopy__(self)

		# Make a smaller matrix (g) to be rotated.
		# How we define that matrix' size must depend on the distance from the end of the protein
		if itUp:
			# Length of the rotation matrix
			n = self.n - x
		else:
			# Length of the rotation matrix
			n = x - 1
		g = np.zeros((2*n + 1, 2*n + 1))

		if itUp:
			# Place the monomer numbers higher than x in the rotation matrix and remove them from the stationary copy
			Gcoords = pivotCoords
			for i in range(x+1, self.n+1):
				# Find the next monomer to place in g
				Gcoords = self.G.searchAdjacent(Gcoords, i)
				# Find the relative coordinates between this point and the pivot point
				relCoordRow = pivotCoords[0] - Gcoords[0]
				relCoordCol = Gcoords[1] - pivotCoords[1]
				gcoords = np.array([n - relCoordRow, n + relCoordCol])
				# Place correctly in g
				g[gcoords[0]][gcoords[1]] = i
				# Remove the monomer from G
				G[Gcoords[0]][Gcoords[1]] = 0
		else:
			# Place the monomer number lower than x in the rotation matrix and remove the from the stationary copy
			Gcoords = pivotCoords
			for i in range(x-1, 0, -1):
				# Find the next monomer to place in g
				Gcoords = self.G.searchAdjacent(Gcoords, i)
				# Find the relative coordinates between this point and the pivot point
				relCoordRow = pivotCoords[0] - Gcoords[0]
				relCoordCol = Gcoords[1] - pivotCoords[1]
				gcoords = np.array([n - relCoordRow, n + relCoordCol])
				# Place correctly in g
				g[gcoords[0]][gcoords[1]] = i
				# Remove the monomer from G
				G[Gcoords[0]][Gcoords[1]] = 0

		# Rotate g
		if clockwise:
			g = np.rot90(g, -1)
		else:
			g = np.rot90(g, 1)

		if itUp:
			gPivotCoords = np.array([n, n])
			for i in range(x+1, self.n+1):
				# Find the next monomer to check in g
#				gcoords = self.G.searchAdjacent(gcoords, i)
				gcoords = np.argwhere(g == i)[0]
				# Find the relative coordinates between this point and the pivot point
				relCoordRow = gPivotCoords[0] - gcoords[0]
				relCoordCol = gcoords[1] - gPivotCoords[1]
				Gcoords = pivotCoords + [-relCoordRow, relCoordCol]
				# Check if the corresponding element in g is sent to 0. If yes: execute the transfer
				if G[Gcoords[0]][Gcoords[1]] == 0:
					G[Gcoords[0]][Gcoords[1]] = g[gcoords[0]][gcoords[1]]
				else:
					return False
		else:
			gPivotCoords = np.array([n, n])
			for i in range(x-1, 0, -1):
				# Find the next monomer to check in g
#				gcoords = self.G.searchAdjacent(gcoords, i)
				gcoords = np.argwhere(g == i)[0]
				# Find the relative coordinates between this point and the pivot point
				relCoordRow = gPivotCoords[0] - gcoords[0]
				relCoordCol = gcoords[1] - gPivotCoords[1]
				Gcoords = pivotCoords + [-relCoordRow, relCoordCol]
				# Check if the corresponding element in g is sent to 0. If yes: execute the transfer
				if G[Gcoords[0]][Gcoords[1]] == 0:
					G[Gcoords[0]][Gcoords[1]] = g[gcoords[0]][gcoords[1]]
				else:
					return False

		# If it was never denied, we now change the grid. Return True to mark the twist was successful
		self.G.grid = G
		return True

	def randomMonomer(self):
		"""

		:return: Int; A random monomer in the polymer
		"""
		return random.randint(2, self.n - 1)

	def randomTwist(self):
		"""

		:return: None; Performs a successful twist
		"""
		isLegal = False
		while not isLegal:
			randMono = self.randomMonomer()
			randBool = randomBool()
			isLegal = self.twist(randMono, randBool)

