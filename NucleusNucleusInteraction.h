#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <functional>

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Common.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Grid.h>
#include <crpropa/GridTools.h>

#ifdef _OPENMP
	#include "omp.h"
#endif

// particle masses in eV
const double mChargedPion = 1.39570 * 1e8 * crpropa::eV / crpropa::c_squared;
const double mNeutralPion = 1.34977 * 1e8 * crpropa::eV / crpropa::c_squared;
const double mEtaMeson = 5.47862 * 1e8 * crpropa::eV / crpropa::c_squared;
const double mMuon = 1.0566 * 1e8 * crpropa::eV / crpropa::c_squared;
const double tauChargedPion = 2.6033e-8;
const double tauNeutralPion = 8.4e-17;
const double tauEtaMeson = 5e-19; 
const double tauMuon = 2.1969811e-6;

/**
 @class NucleusNucleusInteraction
 @brief Handles pp interactions following the parametrisations from:
 
 Kelner et al. Phys. Rev. D 74 (2006) 034018.
 Kafexhiu et al. Phys. Rev. D 90 (2014) 123014.

 Note that this is a hybrid approach: Kafexhiu is used only to model the gamma-ray
 spectrum due to the pi-0 decay, whereas Kelner+ is used to model everything else.
 */
class NucleusNucleusInteraction : public crpropa::Module
{
public:
	std::vector<double> logFraction;
	std::vector<double> logIncidentEnergy;
	std::vector<double> logProbabilities;
	std::vector<double> neutralPionFraction;
	std::vector<double> chargedPionFraction;
	std::vector<double> etaMesonFraction;
	std::vector<double> probabilities;

	crpropa::ref_ptr<crpropa::ScalarGrid> densityGrid;

	bool haveElectrons;
	bool havePhotons;
	bool haveNeutrinos;
	bool havePions;
	bool haveMuons;
	bool doDecayChargedPion;
	bool doDecayNeutralPion;
	bool doDecayMuon;
	double limit;
	double thinning;
	double normMatterField;
	bool isDensityConstant;

	NucleusNucleusInteraction(double normMatterField = 1., bool photons = false, bool electrons = false, bool neutrinos = false, double thinning = 0, double limit = 0.1);
	NucleusNucleusInteraction(crpropa::ref_ptr<crpropa::ScalarGrid> densityGrid, double normMatterField = 1., bool photons = false, bool electrons = false, bool neutrinos = false, double limit = 0.1);
	void setLimit(double limit);
	void setThinning(double thinning);
	void setHaveElectrons(bool electrons);
	void setHaveNeutrinos(bool neutrinos);
	void setHavePhotons(bool photons);
	void setFieldNorm(double normMatterField);
	void setIsDensityConstant(bool densityConstant);
	void initMesonSpectra();
	void decayAll(bool b);
	void decayNeutralPion(bool b);
	void decayChargedPion(bool b);
	void decayMuon(bool b);
	double crossSection(double energy) const;
	double energyFractionNeutralPion(double energy, double xmin = 0, double xmax = 1) const;
	double energyFractionChargedPion(double energy, double xmin = 0, double xmax = 1) const;
	double energyFractionEtaMeson(double energy, double xmin = 0, double xmax = 1) const;
	double lossLength(int id, double E) const; 
	double lossLength(int id, double E, crpropa::Vector3d position) const; 
	void process(crpropa::Candidate *candidate) const;
	void performInteraction(crpropa::Candidate *candidate) const;
};



class DecayChargedPion
{
public:

	bool haveMuons;
	bool haveNeutrinos;
	double limit;
	double thinning;

	DecayChargedPion(bool muons = true, bool neutrinos = true, double thinning = 0, double limit = 0.05);
	void setLimit(double limit);
	void setThinning(double thinning);
	void setHaveMuons(bool muons);
	void setHaveNeutrinos(bool neutrinos);
	double lossLength(double lorentzFactor) const; 
	double energyFractionMuon() const; 
	void performInteraction(crpropa::Candidate *candidate) const;
};

class DecayNeutralPion {
public:

	bool havePhotons;
	double limit;
	double thinning;

	DecayNeutralPion(bool photons = true, double thinning = 0, double limit = 0.05);
	void setLimit(double limit);
	void setThinning(double thinning);
	void setHavePhotons(bool photons);
	double lossLength(double lorentzFactor) const; 
	double energyFractionPhoton() const; 
	void performInteraction(crpropa::Candidate *candidate) const;
};

class DecayEtaMeson {
public:

	bool havePhotons;
	double limit;
	double thinning;

	DecayEtaMeson(bool photons = true, double thinning = 0, double limit = 0.05);
	void setLimit(double limit);
	void setThinning(double thinning);
	void setHavePhotons(bool photons);
	double lossLength(double lorentzFactor) const; 
	double energyFractionPhoton() const; 
	void performInteraction(crpropa::Candidate *candidate) const;
};


class DecayMuon {
public:

	std::vector<double> logFraction;
	std::vector<double> probabilities;
	std::vector<double> electronFraction;
	std::vector<double> electronNeutrinoFraction;
	std::vector<double> muonNeutrinoFraction;

	bool haveElectrons;
	bool haveNeutrinos;
	double limit;
	double thinning;

	DecayMuon(bool neutrinos = true, bool electrons = true, double thinning = 0, double limit = 0.05);
	void initSpectra();
	void setLimit(double limit);
	void setThinning(double thinning);
	void setHaveElectrons(bool electrons);
	void setHaveNeutrinos(bool neutrinos);
	double lossLength(double lorentzFactor) const; 
	double energyFractionElectron(double xmin = 0, double xmax = 1) const; 
	double energyFractionElectronNeutrino(double xmin = 0, double xmax = 1) const; 
	double energyFractionMuonNeutrino(double xmin = 0, double xmax = 1) const; 
	void performInteraction(crpropa::Candidate *candidate) const;
};


/**
 @class ParticleDecay
 @brief Handles the decay of pions, eta mesons, or muons.
 Functions can be easily adapted for any other particle.
 */
class ParticleDecay : public crpropa::Module {
public:
	std::vector<double> logInitialEnergy;
	std::vector<double> logFraction;
	std::vector<double> probabilities;
	std::vector<double> md_electronFraction;
	std::vector<double> md_muonNeutrinoFraction;
	std::vector<double> md_electronNeutrinoFraction;
	std::vector<double> pd_muonNeutrinoFraction;
	std::vector<double> pd_muonFraction;

	DecayChargedPion *chargedPionDecay;
	DecayNeutralPion *neutralPionDecay;
	DecayEtaMeson *etaMesonDecay;
	DecayMuon *muonDecay;

	std::vector<int> decayParticles;
	bool haveElectrons;
	bool havePhotons;
	bool haveNeutrinos;
	bool haveMuons;
	double limit;
	double thinning;

	ParticleDecay(bool photons = true, bool neutrinos = true, bool electrons = true, bool muons = true, const std::vector<int>& pList = std::vector<int>(), double thinning = 0, double limit = 0.1);
	void setLimit(double limit);
	void setThinning(double thinning);
	void setHaveElectrons(bool electrons);
	void setHaveMuons(bool muons);
	void setHaveNeutrinos(bool neutrinos);
	void setHavePhotons(bool photons);
	void setDecayParticles(std::vector<int> v);
	double lossLength(int id, double lorentzFactor) const; 
	void process(crpropa::Candidate *candidate) const;
	void performInteraction(crpropa::Candidate *candidate) const;
};