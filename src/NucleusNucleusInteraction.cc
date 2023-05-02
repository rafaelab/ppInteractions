#include "ppInteractions/NucleusNucleusInteraction.h"

using namespace crpropa;

NucleusNucleusInteraction::NucleusNucleusInteraction(double normMatterField, double thinning, double limit) : Module() {
	setFieldNorm(normMatterField);
	setLimit(limit);
	setIsDensityConstant(true);
	setThinning(thinning);
	initMesonSpectra();
	setDescription("NucleusNucleusInteraction");
}

NucleusNucleusInteraction::NucleusNucleusInteraction(ref_ptr<Grid1f> grid, double normMatterField, double thinning, double limit) : Module() {
	setFieldNorm(normMatterField);
	setLimit(limit);
	setThinning(thinning);
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


void NucleusNucleusInteraction::setFieldNorm(double x) {
	if (x == 0) 
		x = std::numeric_limits<double>::min();
	normMatterField = x;
}

void NucleusNucleusInteraction::setLimit(double l) {
	limit = l;
}

void NucleusNucleusInteraction::setThinning(double t) {
	thinning = t;
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

double NucleusNucleusInteraction::crossSection(const double& energy) const {
	// Parametrisation from:
	//   Kafexhiu et al. PRD 90 (2014) 123014.
	// Note: only works for en >> Ethr
	// Values are given in mb, hence the 1e-31 factor
	double ethr = 2.797e8 * eV;
	double x = (energy - mass_proton * c_squared) / ethr;
	double res = (30.7 - 0.96 * log(x) + 0.18 * log(x) * log(x)) * pow(1. - pow(1. / x, 1.9), 3);
	return res * 1e-31; 
}

double NucleusNucleusInteraction::energyFractionChargedPion(const double& energy, const double& xmin, const double& xmax) const {
	Random &random = Random::instance();
	return interpolate2d(log10(energy), random.randUniform(xmin, xmax), logIncidentEnergy, probabilities, chargedPionFraction);
}

double NucleusNucleusInteraction::energyFractionNeutralPion(const double& energy, const double& xmin, const double& xmax) const {
	Random &random = Random::instance();
	return interpolate2d(log10(energy), random.randUniform(xmin, xmax), logIncidentEnergy, probabilities, neutralPionFraction);
}

double NucleusNucleusInteraction::energyFractionEtaMeson(const double& energy, const double& xmin, const double& xmax) const {
	Random &random = Random::instance();
	return interpolate2d(log10(energy), random.randUniform(xmin, xmax), logIncidentEnergy, probabilities, etaMesonFraction);
}

double NucleusNucleusInteraction::lossLength(const int& id, const double& energy) const {
	double rate = crossSection(energy) * normMatterField;
	return 1. / rate;
}

double NucleusNucleusInteraction::lossLength(const int& id, const double& energy, const Vector3d& position) const {
	double rate = crossSection(energy) * densityGrid->interpolate(position);
	return 1. / rate;
}

void NucleusNucleusInteraction::process(Candidate* candidate) const {

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

void NucleusNucleusInteraction::performInteraction(Candidate* candidate) const {
	
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