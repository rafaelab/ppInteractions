#include "ppInteractions/DecayNeutralPion.h"

using namespace crpropa;

DecayNeutralPion::DecayNeutralPion(bool photons, double thinning, double limit) {
	setHavePhotons(photons);
	setLimit(limit);
	setThinning(thinning);
}

void DecayNeutralPion::setHavePhotons(bool b) {
	havePhotons = b;
}

void DecayNeutralPion::setLimit(double l) {
	limit = l;
}

void DecayNeutralPion::setThinning(double thinning) {
	thinning = thinning;
}

double DecayNeutralPion::lossLength(double lf) const {
	// Returns the loss length in the lab frame.
	double lifetime;
	lifetime = tauNeutralPion;
	return c_light * lifetime * lf;
}

double DecayNeutralPion::energyFractionPhoton() const {
	return 0.5;
}

void DecayNeutralPion::performInteraction(Candidate *candidate) const {

	double E = candidate->current.getEnergy();  
	double w0 = candidate->getWeight();

	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	// particle disappears after decay
	candidate->setActive(false);

	double f = energyFractionPhoton();
	if (havePhotons) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = w0 / pow(1 - f, thinning);
			candidate->addSecondary(22, E * (1 - f), pos, w);
		} 
		if (random.rand() < pow(f, thinning)) {
			double w = w0 / pow(f, thinning);
			candidate->addSecondary(22, E * f, pos, w);
		} 
	}
}