#include <sstream>
#include <fstream>
#include <cstring>

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
const double mChargedPion = 1.39570 * 1e8;
const double mNeutralPion = 1.34977 * 1e8;
const double mMuon = 1.0566 * 1e8;


class ElectronSpectrum {
public:
    std::vector<double> frac;
    std::vector<double> prob;
	void initSpectrum();
	double energyFraction();
};

class ElectronNeutrinoSpectrum {
public:
    std::vector<double> frac;
    std::vector<double> prob;
	void initSpectrum();
	double energyFraction();
};

class MuonNeutrinoSpectrum1 {
public:
    std::vector<double> frac;
    std::vector<double> prob;
	void initSpectrum();
	double energyFraction();
};

class MuonNeutrinoSpectrum2 {
public:
    std::vector<double> frac;
    std::vector<double> prob;
	void initSpectrum();
	double energyFraction();
};

class GammaSpectrum {
public:
    std::vector<double> frac;
    std::vector<double> prob;
    std::vector<double> energy;
	void initSpectrum(std::string);
	double energyFraction(double energy);
};

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

	ElectronSpectrum *electronSpec;
	ElectronNeutrinoSpectrum *electronNuSpec;
	MuonNeutrinoSpectrum1 *muonNuSpec1;
	MuonNeutrinoSpectrum2 *muonNuSpec2;
	GammaSpectrum *gammaSpec;
	std::string gammaSpecFile;

public:
	/// The parent's constructor need to be called on initialization!
	ProtonProtonInteraction(std::string fieldName, std::string dataDir, double normBaryonField = 1., bool photons = false, bool electrons = false, bool neutrinos = false, double limit = 0.1);
	void process(crpropa::Candidate *candidate) const;
	void performInteraction(crpropa::Candidate *candidate) const;
	void initSpectra();
	void initRate(std::string filename);
	void initFraction(std::string filename);
	void setLimit(double limit);
	void setHaveElectrons(bool electrons);
	void setHaveNeutrinos(bool neutrinos);
	void setHavePhotons(bool photons);
	void setFieldNorm(double normBaryonField);
	void neutralPionChannel(crpropa::Candidate *candidate) const;
	void chargedPionChannel(crpropa::Candidate *candidate) const;
	void etaMesonChannel(crpropa::Candidate *candidate) const;
};

