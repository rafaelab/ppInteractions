#include "ppInteractions/DecayEtaMeson.h"

using namespace crpropa;

DecayEtaMeson::DecayEtaMeson(bool photons, double thinning, double limit) {
	setHavePhotons(photons);
	setLimit(limit);
	setThinning(thinning);
}

void DecayEtaMeson::setHavePhotons(bool b) {
	havePhotons = b;
}

void DecayEtaMeson::setLimit(double l) {
	limit = l;
}

void DecayEtaMeson::setThinning(double thinning) {
	thinning = thinning;
}

double DecayEtaMeson::lossLength(const double& lf) const {
	return c_light * tauEtaMeson * lf;
}

double DecayEtaMeson::energyFractionPhoton(Random& random) const {
	return 0.5;
}

void DecayEtaMeson::performInteraction(Candidate* candidate) const {

	double E = candidate->current.getEnergy();  

	Random& random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	// particle disappears after decay
	candidate->setActive(false);

	double f = energyFractionPhoton(random);
	if (havePhotons) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			candidate->addSecondary(22, E * (1 - f), pos, w);
		} 
		if (random.rand() < pow(f, thinning)) {
			double w = 1. / pow(f, thinning);
			candidate->addSecondary(22, E * f, pos, w);
		} 
	}
}