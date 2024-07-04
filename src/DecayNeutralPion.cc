#include "ppInteractions/DecayNeutralPion.h"

using namespace crpropa;

DecayNeutralPion::DecayNeutralPion(bool photons, double thinning, double limit) {
	setHavePhotons(photons);
	setLimit(limit);
	setThinning(thinning);
	setInteractionTag("Dec111");
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

void DecayNeutralPion::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string DecayNeutralPion::getInteractionTag() const {
	return interactionTag;
}

double DecayNeutralPion::lossLength(const double& lf) const {
	return c_light * tauNeutralPion * lf;
}

double DecayNeutralPion::energyFractionPhoton(Random& random) const {
	return 0.5;
}

void DecayNeutralPion::performInteraction(Candidate* candidate) const {
 
	double E = candidate->current.getEnergy();  
	
	Random& random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	// particle disappears after decay
	candidate->setActive(false);

	double f = energyFractionPhoton(random);
	if (havePhotons) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			candidate->addSecondary(22, E * (1 - f), pos, w, interactionTag);
		} 
		if (random.rand() < pow(f, thinning)) {
			double w = 1. / pow(f, thinning);
			candidate->addSecondary(22, E * f, pos, w, interactionTag);
		} 
	}
}