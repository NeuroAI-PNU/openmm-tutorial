import torch
from openmmtorch import TorchForce
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout


class ForceModule(torch.nn.Module):
    def forward(self, positions):
        return torch.sum(positions**2)

module = torch.jit.script(ForceModule())
force = TorchForce(module)
pdb = PDBFile('input.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
system.addForce(force)

# importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.