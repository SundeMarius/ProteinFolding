import numpy as np
import random


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

		a, b = np.where(self.grid == x)
		return np.array([a[0], b[0]])


class Protein:

	def __init__(self, length):
		self.n = length  # Length of protein
		self.N = length + 2  # Length of grid
		self.G = Grid(self.N)
		self.G.grid[int(np.round(self.N / 2)), int(np.round(self.N / 2 - self.n / 2)): int(np.round(self.N / 2 - \
			self.n / 2)) + self.n] = np.linspace(1, self.n, self.n)  # .astype(np.int16)
		self.midValue = (length + 1) / 2

	def searchAdjacent(self, pivotCoords, targetNr):
		'''

		:param pivotCoords: np.array; The coordinates of the point we are searching next to
		:param targetNr: Int; The number assigned to the point we are searching for
		:return: np.array; The coordinates of the target number
		'''

		# Check column:
		if (pivotCoords[0] != 0) and (self.G.grid[pivotCoords[0] - 1, pivotCoords[1]] == targetNr):
			return np.array([pivotCoords[0] - 1, pivotCoords[1]])

		if (pivotCoords[0] != self.G.gSize) and (self.G.grid[pivotCoords[0] + 1, pivotCoords[1]] == targetNr):
			return np.array([pivotCoords[0] + 1, pivotCoords[1]])

		# Check row:
		if (pivotCoords[1] != 0) and (self.G.grid[pivotCoords[0], pivotCoords[1] - 1] == targetNr):
			return np.array([pivotCoords[0], pivotCoords[1] - 1])
		if (pivotCoords[1] != self.G.gSize) and (self.G.grid[pivotCoords[0], pivotCoords[1] + 1] == targetNr):
			return np.array([pivotCoords[0], pivotCoords[1] + 1])

		# If target number has not been found:
		print("\n\n### ERROR ###\nfunction searchAdjacent cannot find target number next to initial number\n\n")
		return np.array([0, 0])

	def nearestNeighbours(self,x,y):

	def counterClockwiseTransformation(self, reducedCoords):
		'''

		:param reducedCoords: np.array; The relative coordinates [m, n]
		:return: np.array; the transformed relative coordinates in clockwise transformation [-m, n]
		'''

		return np.array([-reducedCoords[1], reducedCoords[0]])

	def clockwiseTransformation(self, reducedCoords):
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
		:param grid: np.array; Square matrix representation of the current molecule
		:param clockwise: Bool; True if clockwise rotation
		:return: Bool; True if twist is legal, else False
		'''

		increasingIncrement = self.isAboveMiddle(x)
		pivotCoords = self.G.findElement(x)
		currentCoords = pivotCoords
		if (increasingIncrement) and (x < self.G.gSize):
			for i in range(x + 1, self.G.gSize-1):
				# Go to next monom:
				currentCoords = self.searchAdjacent(currentCoords, i)
				# Get reduced coords (relative to pivot coords):
				reducedCoords = np.array([pivotCoords[0] - currentCoords[0], currentCoords[1] - pivotCoords[1]])
				# Transform the reduced coords with clockwise or counterclockwise rotation:
				if clockwise:
					transformedCoords = self.clockwiseTransformation(reducedCoords)
				else:
					transformedCoords = self.counterClockwiseTransformation(reducedCoords)
				# Get the endpoint for the given x, and check if it is avaliable:
				endCoords = pivotCoords + transformedCoords
				if (self.G.grid[endCoords[0]][endCoords[1]] != 0):
					return False

		if (not increasingIncrement) and (x > 0):
			for i in range(x - 1, 0, -1):
				# Go to next monom:
				currentCoords = self.searchAdjacent(currentCoords, i)
				# Get reduced coords (relative to pivot coords):
				reducedCoords = np.array([pivotCoords[0] - currentCoords[0], currentCoords[1] - pivotCoords[1]])
				# Transform the reduced coords with clockwise or counterclockwise rotation:
				if clockwise:
					transformedCoords = self.clockwiseTransformation(reducedCoords)
				else:
					transformedCoords = self.counterClockwiseTransformation(reducedCoords)
				# Get the endpoint for the given x, and check if it is avaliable:
				endCoords = pivotCoords + transformedCoords
				if (self.G.grid[endCoords[0]][endCoords[1]] != 0):
					return False

		return True

	def twist(self, x, clockwise):
		'''

		:param x: Int; The number corresponding to the pivot monome
		:param clockwise: Bool; True if rotation is clockwise, else False
		:return: Grid; Grid class object with Grid.grid updated with the new twist
		'''

		assert self.isLegalTwist(x, clockwise)

		increasingIncrement = self.isAboveMiddle(x)
		pivotCoords = self.G.findElement(x)
		currentCoords = pivotCoords

		if (increasingIncrement) and (x < self.G.gSize):
			for i in range(x + 1, self.G.gSize-1):
				# Go to next monom:
				currentCoords = self.searchAdjacent(currentCoords, i)
				# Get reduced coords (relative to pivot coords):
				reducedCoords = np.array([pivotCoords[0] - currentCoords[0], currentCoords[1] - pivotCoords[1]])
				# Transform the reduced coords with clockwise or counterclockwise rotation:
				if clockwise:
					transformedCoords = self.clockwiseTransformation(reducedCoords)
				else:
					transformedCoords = self.counterClockwiseTransformation(reducedCoords)
				# Get the endpoint for the given x, and swap values:
				endCoords = pivotCoords + transformedCoords
				self.G.grid[endCoords[0]][endCoords[1]], self.G.grid[currentCoords[0], currentCoords[1]] =\
					self.G.grid[currentCoords[0], currentCoords[1]], self.G.grid[endCoords[0]][endCoords[1]]

		if (not increasingIncrement) and (x > 0):
			for i in range(x - 1, 0, -1):
				# Go to next monom:
				currentCoords = self.searchAdjacent(currentCoords, i)
				# Get reduced coords (relative to pivot coords):
				reducedCoords = np.array([pivotCoords[0] - currentCoords[0], currentCoords[1] - pivotCoords[1]])
				# Transform the reduced coords with clockwise or counterclockwise rotation:
				if clockwise:
					transformedCoords = self.clockwiseTransformation(reducedCoords)
				else:
					transformedCoords = self.counterClockwiseTransformation(reducedCoords)
				# Get the endpoint for the given x, and swap values:
				endCoords = pivotCoords + transformedCoords
				self.G.grid[endCoords[0]][endCoords[1]], self.G.grid[currentCoords[0], currentCoords[1]] = \
					self.G.grid[currentCoords[0], currentCoords[1]], self.G.grid[endCoords[0]][endCoords[1]]

		return self.G



