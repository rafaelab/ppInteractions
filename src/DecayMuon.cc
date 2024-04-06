#include "ppInteractions/DecayMuon.h"

using namespace crpropa;

DecayMuon::DecayMuon(bool neutrinos, bool electrons, double thinning, double limit) {
	setHaveNeutrinos(neutrinos);
	setHaveElectrons(electrons);
	setLimit(limit);
	setThinning(thinning);
	initSpectra();
}

void DecayMuon::initSpectra() {

	// clear previously loaded tables
	probabilities.clear();
	electronNeutrinoFraction.clear();
	electronFraction.clear();
	muonNeutrinoFraction.clear();

	const int nx = 1000; // number of points to sample energy fraction
	const int np = 10000; // number of points to sample probability array

	// Create an array with energy fractions
	const double xmin = 0.; // minimum energy fraction
	double dx = (1. - xmin) / nx;

	std::vector<double> fraction(nx);
	for (int i = 0; i < nx; i++) {
		fraction[i] = xmin + i;
	}

	// Create an array with probabilities for sampling
	const double pmin = 0; // minimum energy fraction
	const double pmax = 1; // maximum energy fraction
	double dp = (pmax - pmin) / np;
	for (int i = 0; i < np; i++) {
		probabilities.push_back(pmin + i * dp);
	}

	// For a fixed incident proton energy, get distribution of energy fractions.
	std::vector<double> fElectron; // from muon decay
	std::vector<double> fElectronNeutrino; // from muon decay
	std::vector<double> fMuonNeutrino; // from muon decay
	std::vector<double> logFElectron; // from muon decay
	std::vector<double> logFElectronNeutrino; // from muon decay
	std::vector<double> logFMuonNeutrino; // from muon decay

	for (size_t j = 0; j < nx; j++) {
		// Equations from Gaisser's book, 2016, tab. 6.2
		double x = fraction[j];

		// electrons/muon neutrinos from muon decay
		double g0 = 5. / 3. - 3. * pow_integer<2>(x) + 4. / 3. * pow_integer<3>(x);
		double g1 = 1. / 3. - 3. * pow_integer<2>(x) + 8. / 3. * pow_integer<3>(x);
		fElectron.push_back(g0 + g1);
		fMuonNeutrino.push_back(g0 + g1);

		// electron neutrinos from muon decay
		double h0 = 2. - 6. * pow_integer<2>(x) + 4. * pow_integer<3>(x);
		double h1 = -2. + 12. * x - 18. * pow_integer<2>(x) + 8. * pow_integer<3>(x);
		fElectronNeutrino.push_back(h0 + h1);
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

	for (size_t j = 0; j < nx; j++) {
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
	for (size_t j = 0; j < probabilities.size(); j++) {
		double y1 = interpolate(probabilities[j], logFMuonNeutrino, fraction);
		double y2 = interpolate(probabilities[j], logFElectron, fraction);
		double y3 = interpolate(probabilities[j], logFElectronNeutrino, fraction);
		muonNeutrinoFraction.push_back(y1);
		electronFraction.push_back(y2);
		electronNeutrinoFraction.push_back(y3);
	}
}

void DecayMuon::setHaveElectrons(bool b) {
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

double DecayMuon::lossLength(const double& lf) const {
	return c_light * tauMuon * lf;
}

double DecayMuon::energyFractionElectron(Random& random, const double& xmin, const double& xmax) const {
	return interpolate(random.randUniform(xmin, xmax), probabilities, electronFraction);
}

double DecayMuon::energyFractionElectronNeutrino(Random& random, const double& xmin, const double& xmax) const {
	return interpolate(random.randUniform(xmin, xmax),probabilities, electronNeutrinoFraction);
}

double DecayMuon::energyFractionMuonNeutrino(Random& random, const double& xmin, const double& xmax) const {
	return interpolate(random.randUniform(xmin, xmax), probabilities, muonNeutrinoFraction);
}

void DecayMuon::performInteraction(Candidate* candidate) const {

	double E = candidate->current.getEnergy();    
	int id = candidate->current.getId();
	double sign = (id > 0) ? 1 : ((id < 0) ? -1 : 0);

	Random& random = Random::instance();
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
	// // It doesn't conserve energy, but it holds statistically.
	double fe = energyFractionElectron(random, 1e-10, 1);
	double fnue = energyFractionElectronNeutrino(random, 1e-10, 1);
	double fnumu = energyFractionMuonNeutrino(random, 1e-10, 1);
	double ftot = fe + fnue + fnumu;

	if (haveElectrons) {
		if (random.rand() < pow(fe, thinning)) {
			double w = 1. / pow(fe, thinning);
			candidate->addSecondary(sign * 11, E * fe, pos, w);
		}
	}
	if (haveNeutrinos) {
		if (random.rand() < pow(fnue, thinning)) {
			double w = 1. / pow(fnue, thinning);
			candidate->addSecondary(- sign * 12, E * fnue, pos, w);
		}
		if (random.rand() < pow(fnumu, thinning)) {
			double w = 1. / pow(fnumu, thinning);
			candidate->addSecondary(sign * 14, E * fnumu, pos, w);	
		}
	}

}
