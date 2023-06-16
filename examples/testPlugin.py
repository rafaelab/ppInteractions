from crpropa import *
import NucleusNucleusInteraction as ppint

photons = True
neutrinos = True
electrons = True
muons = True
thinning = .99
density = 1e13
d = 1. * kpc # source distance
e = 1e16 * eV # proton energy


obs = Observer()
obs.add(Observer1D())
output = TextOutput('test.txt', Output.Event1D)
output.setEnergyScale(eV)
output.enable(output.WeightColumn)
output.disable(output.CandidateTagColumn)
obs.onDetection(output)

pd = ppint.ParticleDecay(photons, neutrinos, electrons, muons, thinning) # must be called
nn = ppint.NucleusNucleusInteraction(density, thinning)

sim = ModuleList()
sim.add(SimplePropagation())
sim.add(MinimumEnergy(1e10 * eV)) # code is not tested below this energy
sim.add(nn)
sim.add(pd)
sim.add(obs)

source = Source()
source.add(SourceEnergy(e))
source.add(SourcePosition(Vector3d(d, 0, 0)))
source.add(SourceDirection(Vector3d(-1, 0, 0)))
source.add(SourceParticleType(nucleusId(1, 1)))

sim.setShowProgress(True)
sim.run(source, 2000, True)
