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

/**
 @class ElectronSpectrum
 @brief Simple class for the electron spectrum due to pi-+ decay.
 Follows:  Kelner et al. Phys. Rev. D 74 (2006) 034018.
 */
class ElectronSpectrum {
public:
    std::vector<double> frac;
    std::vector<double> prob;
	void initSpectrum();
	double energyFraction();
};

/**
 @class ElectronNeutrinoSpectrum
 @brief Simple class for the electron neutrino spectrum due to pi-+ decay.
 Follows:  Kelner et al. Phys. Rev. D 74 (2006) 034018.
 */
class ElectronNeutrinoSpectrum {
public:
    std::vector<double> frac;
    std::vector<double> prob;
	void initSpectrum();
	double energyFraction();
};

/**
 @class MuonNeutrinoSpectrum1
 @brief Simple class for the muon spectrum due to pi-+ decay.
 This spectrum is for muons due to pi-+ decay.
 Follows:  Kelner et al. Phys. Rev. D 74 (2006) 034018.
 */
class MuonNeutrinoSpectrum1 {
public:
    std::vector<double> frac;
    std::vector<double> prob;
	void initSpectrum();
	double energyFraction();
};

/**
 @class MuonNeutrinoSpectrum2
 @brief Simple class for the muon spectrum due to pi-+ decay.
 This spectrum is for muon anti-neutrinos due to muon decay.
 Follows:  Kelner et al. Phys. Rev. D 74 (2006) 034018.
 */
class MuonNeutrinoSpectrum2 {
public:
    std::vector<double> frac;
    std::vector<double> prob;
	void initSpectrum();
	double energyFraction();
};

/**
 @class GammaSpectrum
 @brief Simple class for the gamma-ray spectrum.
 This spectrum is for muons due to pi-0 decay.
 Follows: Kafexhiu et al. Phys. Rev. D 90 (2014) 123014.
 */
class GammaSpectrum {
public:
    std::vector<double> frac;
    std::vector<double> prob;
    std::vector<double> energy;
	void initSpectrum(std::string);
	double energyFraction(double energy);
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

	ElectronSpectrum *electronSpec;
	ElectronNeutrinoSpectrum *electronNuSpec;
	MuonNeutrinoSpectrum1 *muonNuSpec1;
	MuonNeutrinoSpectrum2 *muonNuSpec2;
	GammaSpectrum *gammaSpec;
	std::string gammaSpecFile;

public:
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
	void posChargedPionChannel(crpropa::Candidate *candidate) const;
	void negChargedPionChannel(crpropa::Candidate *candidate) const;
	void etaMesonChannel(crpropa::Candidate *candidate) const;
	double pionMultiplicity(crpropa::Candidate *candidate);
};

