import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.colors as col
import scipy.spatial.distance as ssd
import copy

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

def randomTwist(protein, T):
    """
    Will generate a random, legal twist, and attempt it. If the resulting energy  is less, the twist is executed.
    Else, the twist is executed anyway if and only if the thermal fluctuations allow it. Else, nothing happens.

    :param protein: Protein; the protein object we want to twist
    :param T: Float; the temperature the twist attempts happens at
    :return: Protein; the (possibly) twisted object
    """
    newProtein = copy.deepcopy(protein)
    kb = 1.38e-23  # J/K (Boltzmann's constant)
    B = 1 / (kb * T)  # Boltzmann-parameter

    isLegal = False
    while not isLegal:
        randMono = newProtein.randomMonomer()
        randBool = randomBool()
        isLegal = newProtein.twist(randMono, randBool)
    # Legal twist has now occured on newProtein
    E1 = protein.E
    # Update the energy of the copy
    newProtein.E = newProtein.calculateEnergy()
    # Update the diameter of the copy
    newProtein.D = newProtein.calculateDiameter()
    E2 = newProtein.E
    r = random.uniform(0, 1)
    if E2 <= E1:
        return newProtein
    elif r < np.exp(-B * (E2 - E1)):
        return newProtein
    else:
        return protein


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
        self.E = 0
        self.D = length
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
        #return energy E
        return E

    def twist(self, x, clockwise):
        """
        This updates the Protein object with a legal twist specified by inputs.

        :param x: int, the number of the pivot monomer
        :param clockwise: bool, True if the rotation is clockwise
        :return: bool, True if the rotation occurred (only if legal)
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
        #self.D = self.calculateDiameter()
        #self.E = self.calculateEnergy()
        return True

    def randomMonomer(self):
        """

        :return: Int; A random monomer in the polymer
        """
        return random.randint(2, self.n - 1)

    def calculateDiameter(self):
        """
        :return: Float; the diameter of the protein, defined as the longest distance between two monomers
        """

        # Add all non-zero elements as points to an array
        points = np.zeros([self.n, 2])
        currentCoords = np.array(self.G.findElement(1))
        for i in range(self.n):
            points[i] = currentCoords
            currentCoords = self.G.searchAdjacent(currentCoords, i + 2)

        # Do some stackoverflow magic
        D = ssd.pdist(points)
        D = ssd.squareform(D)

        return np.nanmax(D)



def isNearestNeighbours(protein, i, j):
    '''
    :param i: Int; monomer nr. i
    :param j: Int; monomer nr. j
    :return: Bool; True if the two monomers are nearest neighbours
    '''
    coordMonomerI = protein.G.findElement(i)
    if protein.G.searchAdjacent(coordMonomerI, j)[0] == protein.G.searchAdjacent(coordMonomerI, j)[1] == 0:
        return False
    else:
        return True

def doubleSumProteinEnergy(l, u, summand, protein):
    # summand: To argumenter (i, j)
    # l = lower limit
    # u = upper limit
    summ = 0
    i = l
    j = l
    while i <= u:
        while j <= u:
            summ += summand(protein, i, j)
            j += 1
        i += 1
        j = l
    return summ

def n(protein, i, j):
    '''
    :param i: Int; monomer nr. i
    :param j: Int; monomer nr. j
    :param protein: Protein;
    :return: Int; 1 if nearest neighbors and |i-j|>1, 0 otherwise
    '''
    if abs(i-j)>1 and isNearestNeighbours(protein, i, j):
        return 1
    else:
        return 0

def energySummand(protein, i, j):
    binding = weakBinding()
    E_ij = binding * n(protein, i, j)    
    return E_ij

def calculateBindingEnergyE(T, protein):
    '''
    :param T: Float; Temperature
    :param polymer: Protein; The polymer
    :return: Float; Binding energy E of the polymer
    '''
    N = protein.n
    E = doubleSumProteinEnergy(1, N, energySummand, protein)
    return E

##Quest 4_1
protein15 = Protein(15)

f=1     ## Sett lik 1 for å plotte med 30000 twists (plotter for 30000*f twists)

def makePlottingVector4_1(protein):
    T = 1500
    twistList = np.arange(0, int(30000*f), 1)
    energyList = np.empty(int(30000*f))
    for i in range(int(30000*f)):
        if i>0 and i%int(600*f)==0:
            T -= 30
        energyList[i] = calculateBindingEnergyE(T, protein)
        protein = randomTwist(protein, T)
        print(i)
        #protein.draw()
    return energyList


#Plot:
twistList = np.arange(0, 30000*f, 1)
tempChangeX = np.arange(0, 30000*f, 600*f)
tempChangeY = np.zeros(50)
energyList = makePlottingVector4_1(protein15)

plt.grid()
plt.xlabel(r'$Twists$')
plt.ylabel(r'$E$ (J)')
plt.suptitle('Exercise 4.1')

plot1, = plt.plot(twistList, energyList, lw=0.5)
plot2, = plt.plot(tempChangeX, tempChangeY, '.')
plt.legend([plot1, plot2],['Energy', 'Temperature changes'])
#plt.xlim(0,max(xliste)*1.05)  # Setter grensene for x-aksen.
#plt.ylim(min(energyList)*(1.05),-min(energyList)*(1.05)) # Setter grensene for y-aksen.
plt.show()
""""""