# https://openmm.github.io/openmm-cookbook/latest/notebooks/cookbook/Merging%20Molecules.html

from openmm.app import *

pdb1 = PDBFile('villin.pdb')
pdb2 = PDBFile('ala_ala_ala.pdb')

modeller = Modeller(pdb1.topology, pdb1.positions)
modeller.add(pdb2.topology, pdb2.positions)
mergedTopology = modeller.topology
mergedPositions = modeller.positions

#print(mergedTopology)
#print(mergedPositions)