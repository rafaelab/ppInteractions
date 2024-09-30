from crpropa import *
# import NucleusNucleusInteraction as ppint
import ppInteractions as ppint

photons = True
neutrinos = True
electrons = True
muons = True
thinning = .99
d = 1. * kpc # source distance
e = 1e16 * eV # proton energy


## example constant density
# density = ConstantDensity(1e13, 2, .5)

# example density on grid
filename = '/Users/rab/Dropbox/analyses/interactions_implementation/data/H2_dens_mean_BEG03.txt'
origin = Vector3d(-15.96875, -15.96875, -0.46875) * kpc
nX, nY, nZ = 512, 512, 16
spacing = 0.0625 * kpc
grid = Grid1f(origin, nX, nY, nZ, spacing)
grid.setClipVolume(True) # return 0 outside of the volume
loadGridFromTxt(grid, filename)
density = DensityGrid(grid)


obs = Observer()
obs.add(Observer1D())
output = TextOutput('test.txt', Output.Event1D)
output.setEnergyScale(eV)
output.enable(output.WeightColumn)
output.disable(output.CandidateTagColumn)
obs.onDetection(output)

pd = ppint.ParticleDecay(photons, neutrinos, electrons, muons, thinning) # must be called
nn = ppint.NucleusNucleusInteraction(density, thinning)
nn.setNormDensityField(1e-5)

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
