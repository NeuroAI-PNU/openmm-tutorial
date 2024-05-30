from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from datetime import datetime

# http://docs.openmm.org/latest/userguide/application/02_running_sims.html#a-first-example

tik = datetime.now()
pdb = PDBFile('input.pdb')

forcefield = ForceField('amber14-all.xml','amber14/tip3pfb.xml')
   # use https://github.com/openmm/openmmforcefields
   # see also: https://ambermd.org/AmberTools.php and https://simtk.org/projects/sander_openmm
   
                        

system = forcefield.createSystem(pdb.topology, 
                                 nonbondedMethod=PME, # particle mesh ewald fro long range electrostatics interaction
                                 nonbondedCutoff=1*nanometer, #1 nm (10 angstrom) cutoff for nonbonded interactions
                                 constraints=HBonds) # all bonds involving hydrogen are constrained

# for Langevin dynamics
# Temperature, friction coefficient, time step
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

# for CUDA
platform = Platform.getPlatformByName('CUDA')


simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy() # minimize energy at a starting point
simulation.reporters.append(PDBReporter('output.pdb', 
                                        1000)) # every 1000 steps, write the coordinates to a file
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)
tok = datetime.now()
print("Time for gpu: ", tok - tik)


tik = datetime.now()
pdb = PDBFile('input.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)
tok = datetime.now()
print("Time for cpu: ", tok - tik)
