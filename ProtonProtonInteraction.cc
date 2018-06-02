#include "ProtonProtonInteraction.h"

using namespace crpropa;

ProtonProtonInteraction::ProtonProtonInteraction(std::string fieldName, std::string dataDir, double normBaryonField, bool photons, bool electrons, bool neutrinos, double limit) : Module() {
	std::string filename1 = dataDir + "rate_" + fieldName + ".txt";
	std::string filename2 = dataDir + "enFracPi_" + fieldName + ".txt";
	pionSpecFile = dataDir + "enFracGamma_" + fieldName + ".txt";

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
	LeptonSpectrum *positronSpec = new LeptonSpectrum(-11);
	LeptonSpectrum *muonNuSpec = new LeptonSpectrum(14);
	LeptonSpectrum *muonAntiNuSpec = new LeptonSpectrum(-14);
	LeptonSpectrum *electronNuSpec = new LeptonSpectrum(12);
	pionSpec = PionSpectrum(pionSpecFile);
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

	// select position of secondary within step
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	// check if in tabulated energy range
	if (en < tabFracEnergy.front() or (en > tabFracEnergy.back()))
		return;

	// Use Kafexhiu et al. 2014 approach to compute gamma-ray fluxes.
	// Use Kelner et al. 2006 to scale lepton spectra wrt gamma rays.
	int mult = neutralPionMultiplicity(candidate);
	double eps = 1e-6;
	double fmin = 1e-6;
	double fmax = 1;
	int i = 0;
	double xtot = 0;
	double dEtot_gamma = 0;
	double dEtot_charged = 0;
	while(i < mult and xtot < 1) {
		double x_pi0 = interpolate2d(random.randUniform(fmin, fmax), en, tabFracProb, tabFracEnergy, tabFrac);
		double en_pi0 = x_pi0 * en;
		double x_gamma = 0.5;
		if (havePhotons) {
			candidate->addSecondary(22, x_gamma * en_pi0, pos);
			candidate->addSecondary(22, x_gamma * en_pi0, pos);
		}
		fmax -= x_pi0;
		xtot += x_pi0;
		dEtot_gamma += xtot * en_pi0;
		i++;
	} 

	// // Now compute the secondary leptons.
	// // Use the relation wrt to photons.
	// fmax = 1;
	// i = 0;
	// xtot = 0;
	// mult = 

	// while(i < mult and xtot < 1) {
	// 	double x_pi0 = interpolate2d(random.randUniform(fmin, fmax), en, tabFracProb, tabFracEnergy, tabFrac);
	// 	double en_pi0 = x_pi0 * en;
	// 	double x_gamma = 0.5;
	// 	double x_numu2 = muonNuSpec2.energyFraction(fmin, fmax) * en_pi0;
	// 	double x_e = electronSpec.energyFraction(fmin, fmax) * en_pi0;
	// 	double x_nue = electronNuSpec.energyFraction(fmin, fmax) * en_pi0;
	// 	// double x_numu1 = muonNuSpec1.energyFraction(fmin, fmax) * en_pi0;
	// 	double x_numu1 = 1 - x_e - x_nue - x_numu2;
	// 	if (havePhotons) {
	// 		candidate->addSecondary(22, x_gamma * en_pi0, pos);
	// 		candidate->addSecondary(22, x_gamma * en_pi0, pos);
	// 	}
	// 	if (haveElectrons) {
	// 		candidate->addSecondary(11, x_e * en_pi0, pos);
	// 	}
	// 	if (haveNeutrinos) {
	// 		candidate->addSecondary( 12, x_nue * en_pi0, pos);
	// 		candidate->addSecondary(-14, x_numu1 * en_pi0, pos);
	// 		candidate->addSecondary( 14, x_numu2 * en_pi0, pos);
	// 	}
	// 	fmax -= x_pi0;
	// 	xtot = x_gamma * 2 + x_numu1 + x_numu2 + x_e + x_nue;
	// 	dEtot += xtot * en_pi0;
	// 	i++;
	// } 

	// // \pi^+ -> \mu^+ + \nu_\mu
	candidate->current.setEnergy(en - dEtot_gamma);
	
}

int ProtonProtonInteraction::neutralPionMultiplicity(Candidate *candidate) const {
	// Parametrisation from Kafexhiu et al. 2014 for the pi0 multiplicity
	// Parameters are for GEANT4 (will be changed in future releases).
	double en = candidate->current.getEnergy() - mass_proton * c_squared;
	double a1 = 0.728;
	double a2 = 0.596;
	double a3 = 0.491;
	double a4 = 0.2503;
	double a5 = 0.117;
	double csi = (en - 3e9 * eV) / (mass_proton * c_squared);
	double res = a1 * pow(csi, a4) * (1 + exp(-a2 * pow(csi, a5))) * (1 - exp(-a3 * pow(csi, 0.25)));
	return (int) std::floor(res);
}

int ProtonProtonInteraction::chargedPionMultiplicity(Candidate *candidate) const {
	// Parametrisation from Kelner et al. 2006 (Eq. 4).
	// Parameters are for QGSJet.
	double en = candidate->current.getEnergy();
	double x = en / mass_proton;


	double L = log(en / (1e12 * eV));
	double Bpi = 5.58 + 0.78 * L + 0.10 * L * L;
	double r = 3.1 / pow(Bpi, 1.5);
	double alpha = 0.89 / (sqrt(Bpi) * (1 - exp(-0.33 * Bpi)));
	double a = 1 - pow(x, alpha);
	double b = pow(1 + r * pow(x, alpha), 3);
	// return (int) std::floor(Bpi * pow(a / b, 4));
	double pf = 4 * alpha * Bpi * pow(x, alpha - 1);
	double f1 = pow((1 - pow(x, alpha)) / pow(1 + r * pow(x, alpha), 3), 4);
	double f2 = 1 / (1 - pow(x, alpha)) + 3 * r / (1 + r * pow(x, alpha));
	double f3 = sqrt(1 - mChargedPion / (x * en));
	std::cout << pf * f1 *f2 * f3 << std::endl;
	return (int) std::floor(pf * f1 *f2 * f3 * x);
}

PionSpectrum ProtonProtonInteraction::pionDistribution(double en) const {
	return pionSpec;
}

LeptonSpectrum::LeptonSpectrum(int id) {
	prob.clear();
	ratio.clear();
	frac.clear();
	switch(id) {
		case -11:
			positronDistribution();
		case 12:
			electronNeutrinoDistribution();
		case 14:
			muonNeutrinoDistribution();
		case -14:
			muonAntiNeutrinoDistribution();
		default:
			throw std::invalid_argument("Unknown id for lepton.");
	}
	setLepton(id);
}

void LeptonSpectrum::setLepton(int id) {
	leptonId = id;
}

void LeptonSpectrum::setFraction(std::vector<double> fraction) {
	frac = fraction;
}

void LeptonSpectrum::setRatio(std::vector<double> ratio) {
	ratio = ratio;
}

void LeptonSpectrum::setProbability(std::vector<double> probability) {
	prob = probability;
}

void LeptonSpectrum::computeMultiplicity(PionSpectrum ps) {

}

double LeptonSpectrum::energyFraction(double pmin, double pmax) const {
	Random &random = Random::instance();
	return interpolate(random.randUniform(pmin, pmax), prob, frac);
}

double LeptonSpectrum::ratioToPhoton(double x) const {
	return interpolate(x, frac, ratio);
}

void LeptonSpectrum::muonNeutrinoDistribution() {
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
			ratio.push_back(2.367);
		} else {
			prob.push_back(res);
			ratio.push_back(0.);
		}
	}
}

void LeptonSpectrum::positronDistribution() { 
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
			ratio.push_back(h1 + h2);
		} else {
			res += g;
			prob.push_back(res);
			ratio.push_back(g);
		}
	}
}

void LeptonSpectrum::electronNeutrinoDistribution() {
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
			ratio.push_back(h1 + h2);
		} else {
			res += g;
			prob.push_back(res);
			ratio.push_back(g);
		}
	}
}

void LeptonSpectrum::muonAntiNeutrinoDistribution() {
	// Initiates the table computed following Eq. 36 from Kelner et al. 2006.
	positronDistribution();
}

PionSpectrum::PionSpectrum() {
}

PionSpectrum::PionSpectrum(std::string filename) {
	initSpectrum(filename);
}

void PionSpectrum::initSpectrum(std::string filename) {
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

double PionSpectrum::energyFraction(double en) const {
	Random &random = Random::instance();
	return interpolate2d(random.rand(), en, prob, energy, frac);
}

