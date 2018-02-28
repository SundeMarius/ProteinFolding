import Protein as prot

P = prot.Protein(10)

P.draw()
P.present()
print()

P.twist(8, False)
P.draw()
P.present()
print()

P.twist(6, True)
P.draw()
P.present()
print()

