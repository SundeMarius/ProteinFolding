import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.colors as col
import gc



class Grid:

	def __init__(self, gridSize):
		self.gSize = gridSize
		self.grid = np.zeros((gridSize, gridSize)).astype(np.int16)

	def draw(self):
		print(self.grid)

	def findElement(self, x):
		'''

		:param x: Int; Number of the monomer we want
		:return: np.array; The position of the given x ([i, j])
		'''

		a = np.argwhere(self.grid == x)
		return np.array(a[0])

	def searchAdjacent(self, pivotCoords, targetNr):
		'''

		:param pivotCoords: np.array; The coordinates of the point we are searching next to
		:param targetNr: Int; The number assigned to the point we are searching for
		:return: np.array; The coordinates of the target number
		'''

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
		print("\n### ERROR ###\nfunction searchAdjacent cannot find target number", targetNr, "next to initial number\n")
		return np.array([0, 0])

	def present(self):
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


		'''
		colors = ['white', 'red']
		bounds = [0, 0.5, np.Inf]
		cmap = col.ListedColormap(colors)
		norm = col.BoundaryNorm(bounds, cmap.N)

		plt.imshow(self.grid, interpolation='nearest', cmap=cmap, norm=norm)
		plt.show()
		'''




class Protein:

	def __init__(self, length):
		self.n = length  # Length of protein
		self.N = length + 2  # Length of grid
		self.G = Grid(self.N)
		self.G.grid[int(np.round(self.N / 2)), int(np.round(self.N / 2 - self.n / 2)): int(np.round(self.N / 2 - \
			self.n / 2)) + self.n] = np.linspace(1, self.n, self.n)  # .astype(np.int16)
		self.midValue = (length + 1) / 2

	def clockwiseTransformation(self, reducedCoords):
		'''

		:param reducedCoords: np.array; The relative coordinates [m, n]
		:return: np.array; the transformed relative coordinates in clockwise transformation [-m, n]
		'''

		return np.array([-reducedCoords[1], reducedCoords[0]])

	def counterClockwiseTransformation(self, reducedCoords):
		'''

		:param reducedCoords: np.array; The relative coordinates [m, n]
		:return: np.array; the transformed relative coordinates in clockwise transformation [m, -n]
		'''

		return np.array([reducedCoords[1], -reducedCoords[0]])

	def isAboveMiddle(self, x):
		'''

		:param x: Int, the number assigned to the pivot point
		:return: Bool, True if x > middle value (average). Randomly True or False if x == middle. Else False
		'''
		if (x == self.midValue):
			return random.randint(0, 1)
		else:
			return (x > self.midValue)

	def isLegalTwist(self, x, clockwise):
		'''

		:param x: Int; The number assigned to the pivot point
		:param clockwise: Bool; True if clockwise rotation
		:return: Bool; True if twist is legal, else False
		'''

		pivotCoords = self.G.findElement(x)

		# Make a copy of the swapping area and the non-swapping area

		# Stationary copy
		GG = self.G.grid.__deepcopy__(self)
		GG2 = self.G.grid.__deepcopy__(self)
		currentCoords = pivotCoords

		if self.isAboveMiddle(x):
			# Rotation copy
			nn = self.n - x
			gg = Grid(2 * nn + 1)
			gg.grid = GG[pivotCoords[0] - nn: pivotCoords[0] + nn + 1, pivotCoords[1] - nn: pivotCoords[1] + nn + 1]
			gCurrentCoords = np.array([nn, nn])
			gg.grid[gCurrentCoords[0], gCurrentCoords[1]] = 0
			for i in range(x - 1, x - nn - 1, -1):
				gCurrentCoords = gg.searchAdjacent(gCurrentCoords, i)
				gg.grid[gCurrentCoords[0]][gCurrentCoords[1]] = 0
			# Do the rotation
			if clockwise:
				gg.grid = np.rot90(gg.grid, -1)
			else:
				gg.grid = np.rot90(gg.grid, 1)
			# Remove excess elements from original grid
			for i in range(x + 1, self.n + 1):
				currentCoords = self.G.searchAdjacent(currentCoords, i)
				GG2[currentCoords[0]][currentCoords[1]] = 0
			# Check if the area the protein now wants to go to is empty
			gCurrentCoords = np.array([nn, nn])
			for i in range(x + 1, self.n + 1):
				gCurrentCoords = gg.searchAdjacent(gCurrentCoords, i)
				if (GG2[pivotCoords[0] + gCurrentCoords[0] - nn][pivotCoords[1] + gCurrentCoords[1] - nn] != 0):
					return False

		else:
			# Rotation copy
			nn = x - 1
			gg = Grid(2 * nn + 1)
			gg.grid = GG[pivotCoords[0] - nn: pivotCoords[0] + nn + 1, pivotCoords[1] - nn: pivotCoords[1] + nn + 1]
			gCurrentCoords = np.array([nn, nn])
			gg.grid[gCurrentCoords[0], gCurrentCoords[1]] = 0
			for i in range(x + 1, x + nn + 1):
				gCurrentCoords = gg.searchAdjacent(gCurrentCoords, i)
				gg.grid[gCurrentCoords[0], gCurrentCoords[1]] = 0
			# Do the rotation
			if clockwise:
				gg.grid = np.rot90(gg.grid, -1)
			else:
				gg.grid = np.rot90(gg.grid, 1)
			# Remove excess elements from original grid
			for i in range(x - 1, 0, -1):
				currentCoords = self.G.searchAdjacent(currentCoords, i)
				GG2[currentCoords[0], currentCoords[1]] = 0
			# Check if the area the protein now wants to go to is empty
			gCurrentCoords = np.array([nn, nn])
			for i in range(x - 1, 0, -1):
				gCurrentCoords = gg.searchAdjacent(gCurrentCoords, i)
				if (GG2[pivotCoords[0] + gCurrentCoords[0] - nn][pivotCoords[1] + gCurrentCoords[1] - nn] != 0):
					return False
		return True

	def twist(self, x, clockwise):
		'''

		:param x: Int; The number corresponding to the pivot monome
		:param clockwise: Bool; True if rotation is clockwise, else False
		:return: Grid; Grid class object with Grid.grid updated with the new twist
		'''

		assert self.isLegalTwist(x, clockwise)

		pivotCoords = self.G.findElement(x)

		# Make a copy of the swapping area and the non-swapping area

		# Stationary copy
		GG = self.G.grid.__deepcopy__(self)
		currentCoords = pivotCoords

		if self.isAboveMiddle(x):
			# Rotation copy
			nn = self.n - x
			gg = Grid(2*nn + 1)
			gg.grid = GG[pivotCoords[0]-nn: pivotCoords[0]+nn+1, pivotCoords[1]-nn: pivotCoords[1]+nn+1]
			gCurrentCoords = np.array([nn, nn])
			gg.grid[gCurrentCoords[0], gCurrentCoords[1]] = 0
			for i in range(x-1, x-nn-1, -1):
				gCurrentCoords = gg.searchAdjacent(gCurrentCoords, i)
				gg.grid[gCurrentCoords[0]][gCurrentCoords[1]] = 0
			# Do the rotation
			if clockwise:
				gg.grid = np.rot90(gg.grid, -1)
			else:
				gg.grid = np.rot90(gg.grid, 1)
			#Remove excess elements from original grid
			for i in range(x+1, self.n + 1):
				currentCoords = self.G.searchAdjacent(currentCoords, i)
				self.G.grid[currentCoords[0]][currentCoords[1]] = 0
			# Insert rotated matrix into updated grid
			gCurrentCoords = np.array([nn, nn])
			for i in range(x+1, self.n + 1):
				gCurrentCoords = gg.searchAdjacent(gCurrentCoords, i)
				self.G.grid[pivotCoords[0]+gCurrentCoords[0]-nn][pivotCoords[1]+gCurrentCoords[1]-nn] = gg.grid[gCurrentCoords[0]][gCurrentCoords[1]]

		else:
			# Rotation copy
			nn = x - 1
			gg = Grid(2*nn + 1)
			gg.grid = GG[pivotCoords[0]-nn: pivotCoords[0]+nn+1, pivotCoords[1]-nn: pivotCoords[1]+nn+1]
			gCurrentCoords = np.array([nn, nn])
			gg.grid[gCurrentCoords[0], gCurrentCoords[1]] = 0
			for i in range(x+1, x + nn + 1):
				gCurrentCoords = gg.searchAdjacent(gCurrentCoords, i)
				gg.grid[gCurrentCoords[0], gCurrentCoords[1]] = 0
			# Do the rotation
			if clockwise:
				gg.grid = np.rot90(gg.grid, -1)
			else:
				gg.grid = np.rot90(gg.grid, 1)
			# Remove excess elements from original grid
			for i in range(x-1, 0, -1):
				currentCoords = self.G.searchAdjacent(currentCoords, i)
				self.G.grid[currentCoords[0], currentCoords[1]] = 0
			# Insert rotated matrix into updated grid
			gCurrentCoords = np.array([nn, nn])
			for i in range(x-1, 0, -1):
				gCurrentCoords = gg.searchAdjacent(gCurrentCoords, i)
				self.G.grid[pivotCoords[0]+gCurrentCoords[0]-nn][pivotCoords[1]+gCurrentCoords[1]-nn] = gg.grid[gCurrentCoords[0]][gCurrentCoords[1]]
			# Release unreferenced memory
			gc.collect()

	def draw(self):
		self.G.draw()

	def present(self):
		self.G.present()

