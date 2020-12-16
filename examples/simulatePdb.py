from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time as t

t_start = t.perf_counter()
print("pdb = PDBFile('input.pdb')")

pdb = PDBFile('input.pdb')

print(t.perf_counter()-t_start)
t_start = t.perf_counter()
print("forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')")

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

print(t.perf_counter()-t_start)
t_start = t.perf_counter()
print("system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)")

system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

print(t.perf_counter()-t_start)
t_start = t.perf_counter()
print("integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)")

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

print(t.perf_counter()-t_start)
t_start = t.perf_counter()
print("simulation = Simulation(pdb.topology, system, integrator)")

simulation = Simulation(pdb.topology, system, integrator)

print(t.perf_counter()-t_start)
t_start = t.perf_counter()
print("simulation.context.setPositions(pdb.positions)")

simulation.context.setPositions(pdb.positions)

print(t.perf_counter()-t_start)
t_start = t.perf_counter()
print("simulation.minimizeEnergy()")

simulation.minimizeEnergy()

print(t.perf_counter()-t_start)
t_start = t.perf_counter()
print("simulation.reporters.append(PDBReporter('output.pdb', 1000))")

simulation.reporters.append(PDBReporter('output.pdb', 1000))

print(t.perf_counter()-t_start)
t_start = t.perf_counter()
print("simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))")

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

print(t.perf_counter()-t_start)
t_start = t.perf_counter()
print("simulation.step(10000)")

simulation.step(10000)

print(t.perf_counter()-t_start)
print("Done!")
