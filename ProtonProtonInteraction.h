#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Common.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>

#ifdef _OPENMP
	#include "omp.h"
#endif

// particle masses in eV
const double mChargedPion = 1.39570 * 1e8 * crpropa::eV / crpropa::c_squared;
const double mNeutralPion = 1.34977 * 1e8 * crpropa::eV / crpropa::c_squared;
const double mMuon = 1.0566 * 1e8 * crpropa::eV / crpropa::c_squared;
/**
 @class PionSpectrum
 @brief Simple class for the pion spectrum.
 This spectrum is for muons due to pi-0 decay.
 Follows: Kafexhiu et al. Phys. Rev. D 90 (2014) 123014.
 */
class PionSpectrum {
public:
	std::vector<double> energy;
	std::vector<double> frac;
	std::vector<double> ratio;
	std::vector<double> prob;
	PionSpectrum();
	void initSpectrum();
	double energyFraction(double en) const;
	double computeSlopeInInterval(double xmin, double xmax) const;
	int neutralPionMultiplicity(double en) const;
	int chargedPionMultiplicity(double en) const;
};


/**
 @class LeptonSpectrum
 @brief Simple protype class to store information about the spectrum of secondary particles.
 All spectra are written in terms of the gamma-ray spectrum.
 Follows:  Kelner et al. Phys. Rev. D 74 (2006) 034018.
		   Kafexhiu et al. Phys. Rev. D 90 (2014) 123014.
 */
class LeptonSpectrum {
protected:
	static const size_t nsamples = 1000;
public:
	int leptonId;
	std::vector<double> frac;
	std::vector<double> ratio;
	std::vector<double> prob;
	LeptonSpectrum();
	LeptonSpectrum(int id);
	void init();
	void setLepton(int id);
	void setFraction(std::vector<double> fraction);
	void setRatio(std::vector<double> ratioToPhoton);
	void setProbability(std::vector<double> probability);
	void muonNeutrinoDistribution();
	void muonAntiNeutrinoDistribution();
	void electronNeutrinoDistribution();
	void positronDistribution();
	int computeMultiplicity(PionSpectrum *ps, int mult) const;
	double energyFraction(double pmin = 0, double pmax = 1) const;
	double ratioToPhoton(double x) const;
};

/**
 @class ProtonProtonInteraction
 @brief Handles pp interactions following the parametrisations from:
 
 Kelner et al. Phys. Rev. D 74 (2006) 034018.
 Kafexhiu et al. Phys. Rev. D 90 (2014) 123014.

 Note that this is a hybrid approach: Kafexhiu is used only to model the gamma-ray
 spectrum due to the pi-0 decay, whereas Kelner+ is used to model everything else.
 */
class ProtonProtonInteraction : public crpropa::Module
{
protected:
	std::vector<double> tabEnergy;
	std::vector<double> tabRate;
	std::vector<double> tabFrac;		
	std::vector<double> tabFracProb;
	std::vector<double> tabFracEnergy;

	bool haveElectrons;
	bool havePhotons;
	bool haveNeutrinos;
	double limit;
	double normBaryonField;

	LeptonSpectrum *positronSpec;
	LeptonSpectrum *electronNuSpec;
	LeptonSpectrum *muonAntiNuSpec;
	LeptonSpectrum *muonNuSpec;
	PionSpectrum *pionSpec;

public:
	ProtonProtonInteraction(double normBaryonField = 1., bool photons = false, bool electrons = false, bool neutrinos = false, double limit = 0.1);
	void process(crpropa::Candidate *candidate) const;
	void performInteraction(crpropa::Candidate *candidate) const;
	void initSpectra();
	void setLimit(double limit);
	void setHaveElectrons(bool electrons);
	void setHaveNeutrinos(bool neutrinos);
	void setHavePhotons(bool photons);
	void setFieldNorm(double normBaryonField);
	double crossSection(double en) const;
	double lossLength(int id, double energy) const;
};

