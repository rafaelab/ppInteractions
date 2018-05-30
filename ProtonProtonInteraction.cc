#include "ProtonProtonInteraction.h"

using namespace crpropa;

ProtonProtonInteraction::ProtonProtonInteraction(std::string fieldName, std::string dataDir, double normBaryonField, bool photons, bool electrons, bool neutrinos, double limit) : Module() {
	std::string filename1 = dataDir + "rate_" + fieldName + ".txt";
	std::string filename2 = dataDir + "enFracPi_" + fieldName + ".txt";
	gammaSpecFile = dataDir + "enFracGamma_" + fieldName + ".txt";

	initRate(filename1);
	initFraction(filename2);
	initSpectra();
	setFieldNorm(normBaryonField);
	setHaveElectrons(electrons);
	setHaveNeutrinos(neutrinos);
	setHaveNeutrinos(photons);
	setLimit(limit);
	setDescription("ProtonProtonInteraction");
}

void ProtonProtonInteraction::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void ProtonProtonInteraction::setHavePhotons(bool b) {
	havePhotons = b;
}

void ProtonProtonInteraction::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void ProtonProtonInteraction::setFieldNorm(double x) {
	normBaryonField = x;
}

void ProtonProtonInteraction::setLimit(double l) {
	limit = l;
}

void ProtonProtonInteraction::initRate(std::string filename) {

	// clear previously loaded tables
	tabEnergy.clear();
	tabRate.clear();

	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("ProtonProtonInteraction: could not open file " + filename);

	while (infile.good()) {
		if (infile.peek() == '#') {
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			continue;
		}
		double a, b;
		infile >> a >> b;
		if (!infile)
			break;
		tabEnergy.push_back(a * eV);
		tabRate.push_back(b / Mpc);
	}
	infile.close();
}

void ProtonProtonInteraction::initSpectra() {
	electronSpec = new ElectronSpectrum();
	electronNuSpec = new ElectronNeutrinoSpectrum();
	muonNuSpec1 = new MuonNeutrinoSpectrum1();
	muonNuSpec2 = new MuonNeutrinoSpectrum2();
	gammaSpec = new GammaSpectrum();
	gammaSpec->initSpectrum(gammaSpecFile);
}

void ProtonProtonInteraction::initFraction(std::string filename) {

	// clear previously loaded tables
	tabFrac.clear();
	tabFracProb.clear();
	tabFracEnergy.clear();

	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("ProtonProtonInteraction: could not open file " + filename);

	int nl = 2049;
	int nc = 401;
	double pmin = 1e-6;
	double pmax = 1;
	for (int k = 0; k < nc; k++) {
		double dp = log10(pmin) + k * (log10(pmax) - log10(pmin)) / nc;
		tabFracProb.push_back(pow(10, dp));
	}
	double entries[nl][nc+1];

	int i = 0;
	while (!infile.eof()) {
		for (int j = 0; j < nc + 1; j++) {
			double entry = 0;
			infile >> entry;
			entries[i][j] = entry;
		}
		i++;
	}
	infile.close();

	for (int i = 0; i < nl; i++)
		tabFracEnergy.push_back(entries[i][0] * eV);

	for (int j = 1; j < nc; j++) {
		for (int i = 0; i < nl; i++) {
			tabFrac.push_back(entries[i][j]);
		}
	}

}

void ProtonProtonInteraction::process(Candidate *candidate) const {
	// the loop should be processed at least once for limiting the next step
	double step = candidate->getCurrentStep();
	do {
		// check if nucleus
		int id = candidate->current.getId();
		if (not (isNucleus(id)))
			return;

		double E = candidate->current.getEnergy();
		double z = candidate->getRedshift();
		double rate = interpolate(E, tabEnergy, tabRate) * normBaryonField;
		// rate /= (1 + z);  // rate per light travel distance -> rate per comoving distance

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

void ProtonProtonInteraction::performInteraction(Candidate *candidate) const {
	
	double en = candidate->current.getEnergy();
	Random &random = Random::instance();

	// check if in tabulated energy range
	if (en < tabFracEnergy.front() or (en > tabFracEnergy.back()))
		return;

	/* multiplicity considerations go here */

	// decide which meson is produced
	int mesonId;
	double r = random.rand();
	if (r < 1 / 3) {
		mesonId = 1; // pi^+ charged pion
	} else {
		mesonId = 2; // neutral pion
	}

	switch(mesonId) {
		case 1:
			posChargedPionChannel(candidate);
			break;
		case 2:
			neutralPionChannel(candidate);
			break;
		default:
			neutralPionChannel(candidate); // throw
	}
}

double ProtonProtonInteraction::pionMultiplicity(Candidate *candidate) {
	// Parametrisation from Kafexhiu et al. 2014 for the pi0 multiplicity
	// Parameters are for GEANT4 (will be changed in future releases).

	double en = candidate->current.getEnergy() - mass_proton * c_squared;
	double a1 = 0.728;
	double a2 = 0.596;
	double a3 = 0.491;
	double a4 = 0.2503;
	double a5 = 0.117;
	double csi = (en - 3e9 * eV) / (mass_proton * c_squared);
	return a1 * pow(csi, a4) * (1 + exp(-a2 * pow(csi, a5))) * (1 - exp(-a3 * pow(csi, 0.25)));
}

void ProtonProtonInteraction::neutralPionChannel(Candidate *candidate) const {
	double en = candidate->current.getEnergy();
	Random &random = Random::instance();
	double x = interpolate2d(random.rand(), en, tabFracProb, tabFracEnergy, tabFrac);

	double enPion = x * en;
	double enProton = (1 - x) * en;

	if (x < 0 || x > 1)
		std::cerr << "neutralPionChannel: pion energy fraction = " << x << ". Should be between 0-1." << std::endl;

	candidate->current.setEnergy(enProton);


	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	if (havePhotons) {
		std::cout << "photons neutralPion" << std::endl;
		double f = gammaSpec->energyFraction(enPion);
		if (f < 0 || f > 1)
			std::cerr << "neutralPionChannel: gamma-ray energy fraction outside range (0,1) (f=" << f << ")." << std::endl;
		candidate->addSecondary(22, f * enPion, pos);
		candidate->addSecondary(22, (1 - f) * enPion, pos);
	}
}

void ProtonProtonInteraction::posChargedPionChannel(Candidate *candidate) const {
	double en = candidate->current.getEnergy();
	Random &random = Random::instance();
	double x = interpolate2d(random.rand(), en, tabFracProb, tabFracEnergy, tabFrac);

	if (x < 0 || x > 1)
		std::cerr << "chargedPionChannel: pion energy fraction = " << x << ". Should be between 0-1." << std::endl;

	double enPion = x * en;
	double enNeutron = (1 - x) * en;

	candidate->current.setId(nucleusId(1, 0)); // neutron
	candidate->current.setEnergy(enNeutron);

	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	if (haveNeutrinos) {
		std::cout << "neutrinos chargedPion" << std::endl;
		double f1 = electronNuSpec->energyFraction();
		double f2 = muonNuSpec2->energyFraction();
		double f3 = muonNuSpec1->energyFraction();
		if (f1 < 0 || f1 > 1)
			std::cerr << "neutralPionChannel: gamma-ray energy fraction outside range (0,1) (f1=" << f1 << ")." << std::endl;
		if (f2 < 0 || f2 > 1)
			std::cerr << "neutralPionChannel: gamma-ray energy fraction outside range (0,1) (f2=" << f2 << ")." << std::endl;
		if (f3 < 0 || f3 > 1)
			std::cerr << "neutralPionChannel: gamma-ray energy fraction outside range (0,1) (f3=" << f3 << ")." << std::endl;
		candidate->addSecondary( 12, f1 * enPion, pos);
		candidate->addSecondary(-14, f2 * enPion, pos);
		candidate->addSecondary( 14, f3 * enPion, pos);
	}
	if (haveElectrons) {
		std::cout << "electrons chargedPion" << std::endl;
		double f = electronSpec->energyFraction();
		candidate->addSecondary( 11, f * enPion, pos);
	}
}

void ProtonProtonInteraction::negChargedPionChannel(Candidate *candidate) const {
	double en = candidate->current.getEnergy();
	Random &random = Random::instance();
	double x = interpolate2d(random.rand(), en, tabFracProb, tabFracEnergy, tabFrac);

	if (x < 0 || x > 1)
		std::cerr << "chargedPionChannel: pion energy fraction = " << x << ". Should be between 0-1." << std::endl;

	double enPion = x * en;
	double enNeutron = (1 - x) * en;

	candidate->current.setId(nucleusId(1, 0)); 
	candidate->current.setEnergy(enNeutron);

	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	if (haveNeutrinos) {
		std::cout << "neutrinos chargedPion" << std::endl;
		double f1 = electronNuSpec->energyFraction();
		double f2 = muonNuSpec2->energyFraction();
		double f3 = muonNuSpec1->energyFraction();
		if (f1 < 0 || f1 > 1)
			std::cerr << "neutralPionChannel: gamma-ray energy fraction outside range (0,1) (f1=" << f1 << ")." << std::endl;
		if (f2 < 0 || f2 > 1)
			std::cerr << "neutralPionChannel: gamma-ray energy fraction outside range (0,1) (f2=" << f2 << ")." << std::endl;
		if (f3 < 0 || f3 > 1)
			std::cerr << "neutralPionChannel: gamma-ray energy fraction outside range (0,1) (f3=" << f3 << ")." << std::endl;
		candidate->addSecondary(-12, f1 * enPion, pos);
		candidate->addSecondary( 14, f2 * enPion, pos);
		candidate->addSecondary(-14, f3 * enPion, pos);
	}
	if (haveElectrons) {
		std::cout << "electrons chargedPion" << std::endl;
		double f = electronSpec->energyFraction();
		candidate->addSecondary(-11, f * enPion, pos);
	}
}
void ProtonProtonInteraction::etaMesonChannel(Candidate *candidate) const {
	// for now unused

	double en = candidate->current.getEnergy();
	Random &random = Random::instance();
	double x = interpolate2d(random.rand(), en, tabFracProb, tabFracEnergy, tabFrac);

	double enPion = x * en;
	double enProton = (1 - x) * en;

	candidate->current.setEnergy(enProton);

	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	if (havePhotons) {
		candidate->addSecondary(22, 0.5 * enPion, pos);
		candidate->addSecondary(22, 0.5 * enPion, pos);
	}
}

void ElectronSpectrum::initSpectrum() {
	// Electron spectrum; follows Eq. 36 from Kelner et al. 2006.

	size_t n = 1000;
	for (int i = 0; i < n; i++) {
		double x = 0.001 * (i + 1);
		frac.push_back(x);
	}

	double r = pow(mMuon / mChargedPion, 2);
	double res = 0;
	for (auto it = prob.begin(); it != prob.end(); it++) {
		double x = *it;
		double g_1 = (3 - 2 * r) / (9 * (1 - r) * (1 - r));
		double g_2 = 9 * x * x - 6 * log(x) - 4 * x * x * x - 5;
		double g = g_1 * g_2;
		double h1_1 = (3 - 2 * r) / (9 * (1 - r) * (1 - r));
		double h1_2 = 9 * r * r - 6 * log(r) - 4 * r * r * r - 5;
		double h1 = h1_1 * h1_2;
		double h2_1 = (1 + 2 * r) * (r - x) / (9 * r * r);
		double h2_2 = 9 * (r + x) - 4 * (r * r + r * x + x * x);
		double h2 = h2_1 * h2_2;
		if (x <= r) {
			res = res + h1 + h2;
			prob.push_back(res);
		} else {
			res += g;
			prob.push_back(res);
		}
	}
}

double ElectronSpectrum::energyFraction() {
	Random &random = Random::instance();
	return interpolate(random.rand(), prob, frac);
}

void ElectronNeutrinoSpectrum::initSpectrum() {
	// Electron neutrino spectrum; follows Eq. 40 from Kelner et al. 2006.
	// Note that there is a typo in this reference. The correct version is presented
	// in the erratum: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.79.039901

	size_t n = 1000;
	for (int i = 0; i < n; i++) {
		double x = 0.001 * (i + 1);
		frac.push_back(x);
	}

	double r = pow(mMuon / mChargedPion, 2);
	double res = 0;
	for (auto it = prob.begin(); it != prob.end(); it++) {
		double x = *it;
		double g_1 = 2 / (3 * (1 - r) * (1 - r)); 
		double g_2 = (6 * (1 - x) * (1 - x) + r * (5 + 5 * x - 4 * x * x)) * (1 - x) + 6 * r * log(x);
		double g = g_1 * g_2;
		double h1_1 = 2 / (3 * (1 - r) * (1 - r));
		double h1_2 = (1 - r) * (6 - 7 * r + 11 * r * r - 4 * r * r * r) + 6 * r * log(r);
		double h1 = h1_1 * h1_2;
		double h2_1 = 2 * (r - x) / (3 * r * r);
		double h2_2 = (7 * r * r - 4 * r * r * r + 7 * x * r - 4 * x * r * r - 2 * x * x - 4 * r * x * x);
		double h2 = h2_1 * h2_2;
		if (x >= r) {
			res = res + h1 + h2;
			prob.push_back(res);
		} else {
			res += g;
			prob.push_back(res);
		}
	}
}

double ElectronNeutrinoSpectrum::energyFraction() {
	Random &random = Random::instance();
	return interpolate(random.rand(), prob, frac);
}

void MuonNeutrinoSpectrum1::initSpectrum() {
	// Initiates the table computed following Eq. 36 from Kelner et al. 2006.

	ElectronSpectrum *es = new ElectronSpectrum();
	frac = es->frac;
	prob = es->prob;
}

double MuonNeutrinoSpectrum1::energyFraction() {
	Random &random = Random::instance();
	return interpolate(random.rand(), prob, frac);
}

void MuonNeutrinoSpectrum2::initSpectrum() {
	// Muon neutrino due to pion+ decay (pi+ -> mu+ + nu_mu)
	// See Eqs. 23-26 from Kelner et al. 2006.

	size_t n = 1000;
	for (int i = 0; i < n; i++) {
		double x = 0.001 * (i + 1);
		frac.push_back(x);
	}
	double r = pow(mMuon / mChargedPion, 2);
	double res = 0;
	for (auto it = prob.begin(); it != prob.end(); it++) {
		double x = *it;
		if (x <= r) {
			res = res + 2.367;
			prob.push_back(res);
		} else {
			prob.push_back(res);
		}
	}
}

double MuonNeutrinoSpectrum2::energyFraction() {
	Random &random = Random::instance();
	return interpolate(random.rand(), prob, frac);
}

void GammaSpectrum::initSpectrum(std::string filename) {
	// Initiates the table computed following Eq. 11 of Kafexhiu et al. 2014.

	// clear previously loaded tables
	frac.clear();
	prob.clear();

	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("ProtonProtonInteraction: could not open file " + filename);

	int nl = 2049;
	int nc = 401;
	double pmin = 1e-6;
	double pmax = 1;
	for (int k = 0; k < nc; k++) {
		double dp = log10(pmin) + k * (log10(pmax) - log10(pmin)) / nc;
		prob.push_back(pow(10, dp));
	}
	double entries[nl][nc+1];

	int i = 0;
	while (!infile.eof()) {
		for (int j = 0; j < nc + 1; j++) {
			double entry = 0;
			infile >> entry;
			entries[i][j] = entry;
		}
		i++;
	}
	infile.close();

	for (int i = 0; i < nl; i++)
		energy.push_back(entries[i][0] * eV);

	for (int j = 1; j < nc; j++) {
		for (int i = 0; i < nl; i++) {
			frac.push_back(entries[i][j]);
		}
	}
}

double GammaSpectrum::energyFraction(double en) {
	Random &random = Random::instance();
	return interpolate2d(random.rand(), en, prob, energy, frac);
}

