import Protein as prot

'''
We chose some specific rotations to better show off that the protein folding code worked
To chose randomly, we would write something like:

import random
num = random.randint(1, SIZE)
clockwise = random.randint(0, 1)

Inserting num and clockwise as parameters for the protein twist, then generating some new num and
clockwise, and inserting them again. 
'''

SIZE = 10

# Initialize
P = prot.Protein(SIZE)

# Show the initial state. present() plots the protein with matplotlib, and draw() prints the matrix to output
P.present()
P.draw()
print()

# Show state after 1 twist
P.twist(8, False)
P.present()
P.draw()
print()

# Show state after 2 twists
P.twist(6, True)
P.present()
P.draw()
print()
