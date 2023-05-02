#include "ppInteractions/ParticleDecay.h"

using namespace crpropa;


ParticleDecay::ParticleDecay(bool photons, bool electrons, bool neutrinos, bool muons, double thinning, const std::vector<int>& pList, double limit) : Module() {
	setHaveElectrons(electrons);
	setHaveNeutrinos(neutrinos);
	setHavePhotons(photons);
	setHaveMuons(muons);
	setDecayingParticles(pList);
	setLimit(limit);
	setDescription("ParticleDecay");
	muonDecay = new DecayMuon(neutrinos, electrons, thinning, limit);
	chargedPionDecay = new DecayChargedPion(muons, neutrinos, thinning, limit);
	neutralPionDecay = new DecayNeutralPion(photons, thinning, limit);
	etaMesonDecay = new DecayEtaMeson(photons, thinning, limit);
}

void ParticleDecay::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void ParticleDecay::setHaveMuons(bool b) {
	haveMuons = b;
}

void ParticleDecay::setHavePhotons(bool b) {
	havePhotons = b;
}

void ParticleDecay::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void ParticleDecay::setLimit(double l) {
	limit = l;
}

void ParticleDecay::setDecayingParticles(std::vector<int> v) {
	if (v.size() == 0) {
		decayParticles.push_back( -13);
		decayParticles.push_back(  13);
		decayParticles.push_back( 111);
		decayParticles.push_back(-211);
		decayParticles.push_back( 211);
		decayParticles.push_back( 221);
	} else {
		decayParticles = v;
	}
}

void ParticleDecay::setThinning(double thinning) {
	thinning = thinning;
}

double ParticleDecay::lossLength(const int& id, const double& lf) const {
	// Returns the loss length in the lab frame.
	double lifetime;

	switch (id) {
		case 13:
		case -13:
			return muonDecay->lossLength(lf);
		case 111: 
			return neutralPionDecay->lossLength(lf);
		case 211: 
		case -211:
			return chargedPionDecay->lossLength(lf);
		case 221:
			return etaMesonDecay->lossLength(lf);
		default:
			lifetime = std::numeric_limits<double>::max();
			return c_light * lifetime * lf;
	}
}

void ParticleDecay::process(Candidate *candidate) const {
	
	int id = candidate->current.getId();

	// available particle decay information	
	if (id != 111 && id != 221 && fabs(id) != 211 && fabs(id)!= 13)
		return;

	bool process = false;
	for (int i = 0; i < decayParticles.size(); i++) {
		if (id == decayParticles[i]) 
			process = true;
	}
	if (process == false) 
		return;

	double E = candidate->current.getEnergy();
	double z = candidate->getRedshift();

	double mass;
	switch (id) {
		case 221:
			mass = mEtaMeson;
			break;
		case -211:
		case 211:
			mass = mChargedPion;
			break;
		case 111:
			mass = mNeutralPion;
			break;
		case -13:
		case 13:
			mass = mMuon;
			break;
		default:
			std::cout << "Unknow particle Id for ParticleDecay " << id << std::endl;
			break;
	}

	// For some reason, lorentz factors are not being automatically computed (lack of mass info in CRPropa)
	double lf = E / (mass * c_squared);
	double rate = 1. / lossLength(id, lf);
	rate *= pow_integer<2>(1 + z);

	// check for interaction
	Random &random = Random::instance();
	double randDistance = -log(random.rand()) / rate;
	if (candidate->getCurrentStep() > randDistance)
		performInteraction(candidate);
	else
		candidate->limitNextStep(limit / rate);
}

void ParticleDecay::performInteraction(crpropa::Candidate *candidate) const {
	int id = candidate->current.getId();
	switch (id) {
		case 13:
		case -13:
			muonDecay->performInteraction(candidate);
			break;
		case 111: 
			neutralPionDecay->performInteraction(candidate);
			break;
		case 211: 
		case -211:
			chargedPionDecay->performInteraction(candidate);
			break;
		case 221:
			etaMesonDecay->performInteraction(candidate);
			break;
	}
}