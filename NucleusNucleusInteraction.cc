#include "NucleusNucleusInteraction.h"

using namespace crpropa;

NucleusNucleusInteraction::NucleusNucleusInteraction(double normMatterField, bool photons, bool electrons, bool neutrinos, double thinning, double limit) : Module() {
	setFieldNorm(normMatterField);
	setHaveElectrons(electrons);
	setHaveNeutrinos(neutrinos);
	setHavePhotons(photons);
	setLimit(limit);
	setIsDensityConstant(true);
	initMesonSpectra();
	setDescription("NucleusNucleusInteraction");
}

NucleusNucleusInteraction::NucleusNucleusInteraction(ref_ptr<ScalarGrid> grid, double normMatterField, bool photons, bool electrons, bool neutrinos, double limit) : Module() {
	// initSpectra();
	setFieldNorm(normMatterField);
	setHaveElectrons(electrons);
	setHaveNeutrinos(neutrinos);
	setHavePhotons(photons);
	setLimit(limit);
	setDescription("NucleusNucleusInteraction");
	setIsDensityConstant(false);
	if (normMatterField != 1) 
		scaleGrid(grid, normMatterField);
	initMesonSpectra();
	densityGrid = grid;
}

void NucleusNucleusInteraction::setIsDensityConstant(bool b) {
	isDensityConstant = b;
}

void NucleusNucleusInteraction::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void NucleusNucleusInteraction::setHavePhotons(bool b) {
	havePhotons = b;
}

void NucleusNucleusInteraction::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void NucleusNucleusInteraction::setFieldNorm(double x) {
	if (x == 0) 
		x = std::numeric_limits<double>::min();
	normMatterField = x;
}

void NucleusNucleusInteraction::setLimit(double l) {
	limit = l;
}

void NucleusNucleusInteraction::setThinning(double thinning) {
	thinning = thinning;
}

void NucleusNucleusInteraction::decayAll(bool b) {
	this->decayNeutralPion(b);
	this->decayChargedPion(b);
	this->decayMuon(b);
}

void NucleusNucleusInteraction::decayNeutralPion(bool b) {
	doDecayNeutralPion = b;
}

void NucleusNucleusInteraction::decayChargedPion(bool b) {
	doDecayChargedPion = b;
}

void NucleusNucleusInteraction::decayMuon(bool b) {
	doDecayMuon = b;
}

void NucleusNucleusInteraction::initMesonSpectra() {
	// Initialise pion spectra.
	// Follows Eqs. 12 and 13 from:
	//  Kelner et al. 74 (2006) 034018.
	// Note that this parametrisation is implemented for older versions of SIBYLL.
	// QGSJet is not implemented as SYBILL provides a better fit.

	// clear previously loaded tables
	logIncidentEnergy.clear();
	neutralPionFraction.clear();
	chargedPionFraction.clear();
	etaMesonFraction.clear();
	probabilities.clear();

	const int nx = 100000; // number of points to sample energy fraction
	const int np = 100000; // number of points to sample probability array


	// Create an array with incident particle energies for sampling
	const int ne = 120; // number of points to sample
	const double emin = 1e10 * eV; // minimum proton energy
	const double emax = 1e22 * eV; // maximum proton energy
	double dle = log10(emax / emin) / ne;
	for (int i = 0; i < ne; i++) {
		double en = log10(emin) + i * dle;
		logIncidentEnergy.push_back(en);
	}

	// Create an array with energy fractions
	const double xmin = 1e-8; // minimum energy fraction
	const double xmax = 1; // maximum energy fraction
	double dlx = log10(xmax / xmin) / nx;
	for (int i = 0; i < nx + 1; i++) {
		double x_ = log10(xmin) + i * dlx;
		logFraction.push_back(x_);
	}

	// Create an array with probabilities for sampling
	const double pmin = 1e-4; // minimum energy fraction
	const double pmax = 1; // maximum energy fraction
	double dlp = log10(pmax / pmin) / np;
	for (int i = 0; i < np; i++) {
		double p_ = log10(pmin) + i * dlp;
		probabilities.push_back(pow(10, p_));
	}

	for (int i = 0; i < ne; i++) {

		// For a fixed incident proton energy, get distribution of energy fractions.
		std::vector<double> fCharged;
		std::vector<double> fNeutral;
		std::vector<double> fEtaMeson;
		std::vector<double> logFCharged;
		std::vector<double> logFNeutral;
		std::vector<double> logFEtaMeson;
		double ei = pow(10, logIncidentEnergy[i]);

		for (int j = 0; j < nx + 1; j++) {

			// Pions
			double x = pow(10, logFraction[j]);
			double epi = x * ei;
			double l = log(ei / TeV);
			double a = 3.67 + 0.83 * l + 0.075 * l * l;
			double alpha = 0.98 / sqrt(a);
			double bpi = a + 0.25;
			double r = 2.6 / sqrt(a);
			double xalpha = pow(x, alpha);
			double pf_F = 4 * alpha * bpi * xalpha / x;
			double f1_F = pow((1. - xalpha) / (1. + r * xalpha * (1. - xalpha)), 4.);
			double f2_F = 1. / (1. - xalpha) + r * (1. - 2. * xalpha) / (1. + r * xalpha * (1. - xalpha));
			double f3_F_a = sqrt(std::max(0., 1. - (mNeutralPion * c_squared) / epi));
			double f3_F_b = sqrt(std::max(0., 1. - (mChargedPion * c_squared) / epi));
			double f_a = pf_F * f1_F * f2_F * f3_F_a;
			double f_b = pf_F * f1_F * f2_F * f3_F_b; 
			fCharged.push_back(f_b * x);
			fNeutral.push_back(f_a * x);

			// Eta mesons
			double f4_F = (0.55 + 0.028 * log(x)) * std::max(0., 1. - mEtaMeson * c_squared / epi);
			fEtaMeson.push_back(f4_F * fNeutral[j] / x);

			logFCharged.push_back(0);
			logFNeutral.push_back(0);
			logFEtaMeson.push_back(0);
		}
		fCharged.resize(nx);
		fNeutral.resize(nx);
		fEtaMeson.resize(nx);
		logFCharged.resize(nx);
		logFNeutral.resize(nx);
		logFEtaMeson.resize(nx);


		// Compute cumulative distribution
		std::partial_sum(fNeutral.begin(), fNeutral.end(), fNeutral.begin());
		std::partial_sum(fCharged.begin(), fCharged.end(), fCharged.begin());
		std::partial_sum(fEtaMeson.begin(), fEtaMeson.end(), fEtaMeson.begin());


		// Normalise distribution
		double normNeutral = 1. / *std::max_element(fNeutral.begin(), fNeutral.end());
		double normCharged = 1. / *std::max_element(fCharged.begin(), fCharged.end());
		double normEtaMeson = 1. / *std::max_element(fEtaMeson.begin(), fEtaMeson.end());

		for (int j = 0; j < fNeutral.size(); j++) {
			fNeutral[j] *= normNeutral;
			fCharged[j] *= normCharged;
			fEtaMeson[j] *= normEtaMeson;
			logFNeutral[j] = log10(fNeutral[j]);
			logFCharged[j] = log10(fCharged[j]);
			logFEtaMeson[j] = log10(fEtaMeson[j]);
		}

		// Invert distribution and sample from it
		for (int j = 0; j < probabilities.size(); j++) {
			double y1 = interpolate(log10(probabilities[j]), logFNeutral, logFraction);
			double y2 = interpolate(log10(probabilities[j]), logFCharged, logFraction);
			double y3 = interpolate(log10(probabilities[j]), logFEtaMeson, logFraction);
			neutralPionFraction.push_back(pow(10, y1));
			chargedPionFraction.push_back(pow(10, y2));
			etaMesonFraction.push_back(pow(10, y3));
		}
	}
}


double NucleusNucleusInteraction::crossSection(double en) const {
	// Parametrisation from:
	//   Kafexhiu et al. PRD 90 (2014) 123014.
	// Note: only works for en >> Ethr
	// Values are given in mb, hence the 1e-31 factor
	double ethr = 2.797e8 * eV;
	double x = (en - mass_proton * c_squared) / ethr;
	double res = (30.7 - 0.96 * log(x) + 0.18 * log(x) * log(x)) * pow(1. - pow(1. / x, 1.9), 3);
	return res * 1e-31; 
}

double NucleusNucleusInteraction::energyFractionChargedPion(double energy, double xmin, double xmax) const {
	Random &random = Random::instance();
	return interpolate2d(log10(energy), random.randUniform(xmin, xmax), logIncidentEnergy, probabilities, chargedPionFraction);
}

double NucleusNucleusInteraction::energyFractionNeutralPion(double energy, double xmin, double xmax) const {
	Random &random = Random::instance();
	return interpolate2d(log10(energy), random.randUniform(xmin, xmax), logIncidentEnergy, probabilities, neutralPionFraction);
}

double NucleusNucleusInteraction::energyFractionEtaMeson(double energy, double xmin, double xmax) const {
	Random &random = Random::instance();
	return interpolate2d(log10(energy), random.randUniform(xmin, xmax), logIncidentEnergy, probabilities, etaMesonFraction);
}

double NucleusNucleusInteraction::lossLength(int id, double energy) const {
	double rate = crossSection(energy) * normMatterField;
	return 1. / rate;
}

double NucleusNucleusInteraction::lossLength(int id, double energy, Vector3d position) const {
	double rate = crossSection(energy) * densityGrid->interpolate(position);
	return 1. / rate;
}

void NucleusNucleusInteraction::process(Candidate *candidate) const {

	// check if nucleus
	int id = candidate->current.getId();
	if (not (isNucleus(id)))
		return;

	// the loop should be processed at least once for limiting the next step
	double step = candidate->getCurrentStep();
	do {
		double energy = candidate->current.getEnergy();
		double z = candidate->getRedshift();
		double rate;
		if (isDensityConstant == false) {
			Vector3d pos = candidate->current.getPosition();
			rate = 1. / lossLength(id, energy, pos);
		} else
			rate = 1. / lossLength(id, energy);
		// std::cout << kpc / rate << " " << energy / eV << std::endl;

		// find interaction mode with minimum random decay distance
		Random &random = Random::instance();
		double randDistance = std::numeric_limits<double>::max();       
		double d = -log(random.rand()) / rate;
		if (d > randDistance)
			continue;
		randDistance = d;

		// check if interaction doesn't happen
		if (step < randDistance) {
			// limit next step to a fraction of the mean free path
			candidate->limitNextStep(limit / rate);
			return;
		}
		// interact and repeat with remaining step
		performInteraction(candidate);
		step -= randDistance;
	} while (step > 0);
}

void NucleusNucleusInteraction::performInteraction(Candidate *candidate) const {
	
	double energy = candidate->current.getEnergy();
	double w0 = candidate->getWeight();

	Random &random = Random::instance();

	// Select position of secondary within step.
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	// Check if in tabulated energy range.
	if (log10(energy) < logIncidentEnergy.front() or (log10(energy) > logIncidentEnergy.back()))
		return;

	// double x1 = energyFractionNeutralPion(energy, 0., 1.);
	// double x2 = energyFractionChargedPion(energy, 0., 1.);
	// double x3 = energyFractionChargedPion(energy, 0., 1.);
	// double x4 =    energyFractionEtaMeson(energy, 0., 1.);
	double x1 = energyFractionNeutralPion(energy, 1e-10, 1.);
	double x2 = energyFractionChargedPion(energy, 1e-10, 1.);
	double x3 = energyFractionChargedPion(energy, 1e-10, 1.);
	double x4 = energyFractionEtaMeson(energy, 1e-10, 1.);
	double y = (x1 + x2 + x3 + x4);
	double x = 0;
	double r = random.rand();

	// std::cout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << y << std::endl;

	int id;
	if (r < x1 / y) {
		id = 111;
		x = x1;
	}
	else if (r >= x1 / y && r < (x1 + x2) / y) {
		id = 211;
		x = x2;
	}
	else if (r >= (x1 + x2) / y && r < (x1 + x2 + x3) / y) {
		id = -211;
		x = x3;
	}
	else {
		id = 221;
		x = x4;
	}

	if (random.rand() < pow(x, thinning)) {
		double w = w0 / pow(x, thinning);
		candidate->addSecondary(id, x * energy, pos, w);
	}

	candidate->current.setEnergy((1 - x) * energy);
}

///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////
ParticleDecay::ParticleDecay(bool photons, bool electrons, bool neutrinos, bool muons, const std::vector<int>& pList, double thinning, double limit) : Module() {
	setHaveElectrons(electrons);
	setHaveNeutrinos(neutrinos);
	setHavePhotons(photons);
	setHaveMuons(muons);
	setDecayParticles(pList);
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
	//
	havePhotons = b;
}

void ParticleDecay::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void ParticleDecay::setLimit(double l) {
	limit = l;
}

void ParticleDecay::setDecayParticles(std::vector<int> v) {
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

double ParticleDecay::lossLength(int id, double lf) const {
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
	rate *= pow(1 + z, 2);

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


///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////
DecayChargedPion::DecayChargedPion(bool muons, bool neutrinos, double thinning, double limit) {
	setHaveMuons(muons);
	setHaveNeutrinos(neutrinos);
	setThinning(thinning);
	setLimit(limit);
	// setDescription("ChargedPionDecay");
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
	if (f < 1e-10) f = 0.;
	if (f > 1.) f = 1.;
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

///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////
DecayNeutralPion::DecayNeutralPion(bool photons, double thinning, double limit) {
	setHavePhotons(photons);
	setLimit(limit);
	setThinning(thinning);
	// setDescription("NeutralPionDecay");
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

///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////
DecayEtaMeson::DecayEtaMeson(bool photons, double thinning, double limit) {
	setHavePhotons(photons);
	setLimit(limit);
	setThinning(thinning);
	// setDescription("EtaMesonDecay");
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

double DecayEtaMeson::lossLength(double lf) const {
	// Returns the loss length in the lab frame.
	double lifetime;
	lifetime = tauEtaMeson;
	return c_light * lifetime * lf;
}

double DecayEtaMeson::energyFractionPhoton() const {
	return 0.5;
}

void DecayEtaMeson::performInteraction(Candidate *candidate) const {

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

///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////
DecayMuon::DecayMuon(bool neutrinos, bool electrons, double thinning, double limit) {
	setHaveNeutrinos(neutrinos);
	setHaveElectrons(electrons);
	setLimit(limit);
	setThinning(thinning);
	initSpectra();
	// setDescription("MuonDecay");
}

void DecayMuon::initSpectra() {

	// clear previously loaded tables
	logFraction.clear();
	probabilities.clear();
	electronNeutrinoFraction.clear();
	electronFraction.clear();
	muonNeutrinoFraction.clear();

	const int nx = 10000; // number of points to sample energy fraction
	const int np = nx; // number of points to sample probability array

	// Create an array with energy fractions
	const double xmin = 1e-10; // minimum energy fraction
	const double xmax = 1; // maximum energy fraction
	for (int i = 0; i < nx; i++) {
		double dlx = log10(xmax / xmin) / nx;
		double x_ = log10(xmin) + i * dlx;
		logFraction.push_back(x_);
	}

	// Create an array with probabilities for sampling
	const double pmin = xmin; // minimum energy fraction
	const double pmax = 1; // maximum energy fraction
	for (int i = 0; i < np; i++) {
		double dlp = log10(pmax / pmin) / np;
		double p_ = log10(pmin) + i * dlp;
		probabilities.push_back(pow(10, p_));
	}

	// For a fixed incident proton energy, get distribution of energy fractions.
	std::vector<double> fElectron; // from muon decay
	std::vector<double> fElectronNeutrino; // from muon decay
	std::vector<double> fMuonNeutrino; // from muon decay
	std::vector<double> logFElectron; // from muon decay
	std::vector<double> logFElectronNeutrino; // from muon decay
	std::vector<double> logFMuonNeutrino; // from muon decay

	for (int j = 0; j < nx; j++) {
		// Equations from Gaisser's book, 2016, tab. 6.2
		double x = pow(10, logFraction[j]); 

		// electrons/muon neutrinos from muon decay
		double g0 = 5. / 3. - 3. * pow(x, 2.) + 4. / 3. * pow(x, 3.);
		double g1 = 1. / 3. - 3. * pow(x, 2.) + 8. / 3. * pow(x, 3.);
		fElectron.push_back(g0 + g1);
		fMuonNeutrino.push_back(g0 + g1);

		// electron neutrinos from muon decay
		double h0 = 2. - 6. * pow(x, 2.) + 4. * pow(x, 3.);
		double h1 = -2. + 12. * x - 18. * pow(x, 2.) + 8. * pow(x, 3.);
		fElectronNeutrino.push_back(h0 + h1);

		// Simplified expressions
		// double g = 2. + 2. * x * x * (2 * x - 3);
		// double h = 2. * x * (1 - 2 * x + x * x);
		// fElectron.push_back(g);
		// fMuonNeutrino.push_back(g);
		// fElectronNeutrino.push_back(h);
	}
	fElectron.resize(nx);
	fElectronNeutrino.resize(nx);
	fMuonNeutrino.resize(nx);

	// Compute cumulative distributions
	std::partial_sum(fElectron.begin(), fElectron.end(), fElectron.begin());
	std::partial_sum(fElectronNeutrino.begin(), fElectronNeutrino.end(), fElectronNeutrino.begin());
	std::partial_sum(fMuonNeutrino.begin(), fMuonNeutrino.end(), fMuonNeutrino.begin());

	// Normalise distribution
	double normElectron = 1. / *std::max_element(fElectron.begin(), fElectron.end());
	double normElectronNeutrino = 1. / *std::max_element(fElectronNeutrino.begin(), fElectronNeutrino.end());
	double normMuonNeutrino = 1. / *std::max_element(fMuonNeutrino.begin(), fMuonNeutrino.end());

	for (int j = 0; j < nx; j++) {
		fElectron[j] *= normElectron;
		fElectronNeutrino[j] *= normElectronNeutrino;
		fMuonNeutrino[j] *= normMuonNeutrino;
		logFElectron.push_back(log10(fElectron[j]));
		logFElectronNeutrino.push_back(log10(fElectronNeutrino[j]));
		logFMuonNeutrino.push_back(log10(fMuonNeutrino[j]));
	}
	logFElectron.resize(nx);
	logFElectronNeutrino.resize(nx);
	logFMuonNeutrino.resize(nx);

	// Invert distribution and sample from it
	for (int j = 0; j < probabilities.size(); j++) {
		double y1 = interpolate(log10(probabilities[j]), logFMuonNeutrino, logFraction);
		double y2 = interpolate(log10(probabilities[j]), logFElectron, logFraction);
		double y3 = interpolate(log10(probabilities[j]), logFElectronNeutrino, logFraction);
		muonNeutrinoFraction.push_back(pow(10, y1));
		electronFraction.push_back(pow(10, y2));
		electronNeutrinoFraction.push_back(pow(10, y3));
	}
}

void DecayMuon::setHaveElectrons(bool b) {
	//
	haveElectrons = b;
}

void DecayMuon::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void DecayMuon::setLimit(double l) {
	limit = l;
}

void DecayMuon::setThinning(double thinning) {
	thinning = thinning;
}

double DecayMuon::lossLength(double lf) const {
	// Returns the loss length in the lab frame.
	double lifetime;
	lifetime = tauMuon;
	return c_light * lifetime * lf;
}

double DecayMuon::energyFractionElectron(double xmin, double xmax) const {
	Random &random = Random::instance();
	return interpolate(random.randUniform(xmin, xmax), probabilities, electronFraction);
}

double DecayMuon::energyFractionElectronNeutrino(double xmin, double xmax) const {
	Random &random = Random::instance();
	return interpolate(random.randUniform(xmin, xmax),probabilities, electronNeutrinoFraction);
}

double DecayMuon::energyFractionMuonNeutrino(double xmin, double xmax) const {
	Random &random = Random::instance();
	return interpolate(random.randUniform(xmin, xmax), probabilities, muonNeutrinoFraction);
}

void DecayMuon::performInteraction(Candidate *candidate) const {

	double E = candidate->current.getEnergy();    
	double w0 = candidate->getWeight();
	int id = candidate->current.getId();
	double sign = (id > 0) ? 1 : ((id < 0) ? -1 : 0);

	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	// particle disappears after decay
	candidate->setActive(false);

	// /***********************************************/
	// // This is the logical way to do it, using a low tolerance.
	// double tol = .8; // hard-coded for now.
	// double fe = 0;
	// double fnue = 0;
	// double fnumu = 0;
	// double ftot = 100;
	// while (fabs(1 - ftot) >= tol) {
	// 	fe = energyFractionElectron(0, 1);
	// 	fnue = energyFractionElectronNeutrino(0, 1);
	// 	fnumu = energyFractionMuonNeutrino(0, 1);
	// 	ftot = fe + fnue + fnumu;
	// }

	// /***********************************************/
	// // This is the fast way to do it.
	// // It doesn't conserve energy, but it hold statistically.
	double fe = energyFractionElectron(0, 1);
	double fnue = energyFractionElectronNeutrino(0, 1);
	double fnumu = energyFractionMuonNeutrino(0, 1);
	double ftot = fe + fnue + fnumu;

	// std::cout << fe << " " << fnue << " " << fnumu << " " << ftot <<  std::endl;

	if (haveElectrons) {
		if (random.rand() < pow(fe, thinning)) {
			double w = w0 / pow(fe, thinning);
			candidate->addSecondary(sign * 11, E * fe, pos, w);
		}
	}
	if (haveNeutrinos) {
		if (random.rand() < pow(fnue, thinning)) {
			double w = w0 / pow(fnue, thinning);
			candidate->addSecondary(- sign * 12, E * fnue, pos, w);
		}
		if (random.rand() < pow(fnumu, thinning)) {
			double w = w0 / pow(fnumu, thinning);
			candidate->addSecondary(sign * 14, E * fnumu, pos, w);	
		}
	}

}
