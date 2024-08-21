# https://openmm.github.io/openmm-cookbook/latest/notebooks/cookbook/Changing%20Temperature%20and%20Pressure.html
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

pdb = PDBFile('./data/input.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
#system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
#                                 nonbondedCutoff=1*nanometer, constraints=HBonds)
''' save the the forcefield system to a file
with open('./results/system.xml', 'w') as f:
    f.write(XmlSerializer.serialize(system))
'''
# load the forcefield system from a file
with open('./results/system.xml', 'r') as f:
    system = XmlSerializer.deserialize(f.read())

platform = Platform.getPlatformByName('CUDA')

integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

# Change the temperature
for i in range(10):
    integrator.setTemperature((300-30*i)*kelvin)
    integrator.step(10)

# When we use barostat, MC acceptance criteria depends on temperature.
# So, we need to input the temperature to the barostat and integrator.
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
simulation.context.reinitialize(preserveState=True)

# We can change the temperature and pressure with setParameter method.
for i in range(10):
    temperature = (300-30*i)*kelvin
    pressure = (0.5+0.01*i)*bar
    simulation.context.setParameter(MonteCarloBarostat.Temperature(), temperature)
    simulation.context.setParameter(MonteCarloBarostat.Pressure(), pressure)
    integrator.setTemperature(temperature)
    integrator.step(10)

simulation.minimizeEnergy()
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)