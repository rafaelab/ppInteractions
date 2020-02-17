#include "ppInteractions/DecayChargedPion.h"

using namespace crpropa;

DecayChargedPion::DecayChargedPion(bool muons, bool neutrinos, double thinning, double limit) {
	setHaveMuons(muons);
	setHaveNeutrinos(neutrinos);
	setThinning(thinning);
	setLimit(limit);
}

void DecayChargedPion::setHaveMuons(bool b) {
	haveMuons = b;
}

void DecayChargedPion::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void DecayChargedPion::setLimit(double l) {
	limit = l;
}

void DecayChargedPion::setThinning(double thinning) {
	thinning = thinning;
}

double DecayChargedPion::lossLength(double lf) const {
	// Returns the loss length in the lab frame.
	double lifetime;
	lifetime = tauChargedPion;
	return c_light * lifetime * lf;
}

double DecayChargedPion::energyFractionMuon() const {
	// Random &random = Random::instance();
	// double r = 1 - pow(mMuon / mChargedPion, 2);
	// return random.randUniform(0., 1 - r);
	return pow(mMuon / mChargedPion, 2.);
}

void DecayChargedPion::performInteraction(Candidate *candidate) const {

	int id = candidate->current.getId();
	double sign = (id > 0) ? 1 : ((id < 0) ? -1 : 0);
	double E = candidate->current.getEnergy();  
	double w0 = candidate->getWeight();

	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	// particle disappears after decay
	candidate->setActive(false);

	double f = energyFractionMuon();
	if (haveMuons) {
		if (random.rand() < pow(1 - f, thinning)) {
			double w = w0 / pow(1 - f, thinning);
			candidate->addSecondary(sign * 13, E * (1 - f), pos, w);
		} 
	}
	if (haveNeutrinos) {
		if (random.rand() < pow(f, thinning))  {
			double w = w0 / pow(f, thinning);
			candidate->addSecondary(-sign * 14, E * f, pos, w);
		} 
	}
}
