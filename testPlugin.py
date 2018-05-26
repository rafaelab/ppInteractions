from crpropa import *
import ProtonProtonInteraction as pp


dataDir =  '../../data/ProtonProtonInteractions/'
field = 'PowerLaw-2.50'

obs = Observer()
obs.add(ObserverPoint())
output = TextOutput('test.txt', Output.Event1D)
output.setEnergyScale(eV)
obs.onDetection(output)

sim = ModuleList()

sim.add(SimplePropagation())
sim.add(MaximumTrajectoryLength(100 * Mpc))
sim.add(MinimumEnergy(1e9 * eV))
sim.add(pp.ProtonProtonInteraction(field, dataDir, 1e-6, True, True, True))
sim.add(obs)

source = Source()
source.add(SourceEnergy(1e18 * eV))
source.add(SourcePosition(Vector3d(100 * Mpc, 0, 0)))
source.add(SourceDirection(Vector3d(-1, 0, 0)))
source.add(SourceParticleType(nucleusId(1, 1)))

sim.setShowProgress(True)
sim.run(source, 2000, True)
