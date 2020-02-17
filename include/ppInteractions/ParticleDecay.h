#ifndef PARTICLEDECAY_H
#define PARTICLEDECAY_H

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Common.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Grid.h>
#include <crpropa/GridTools.h>
#include <crpropa/Candidate.h>

#include "ppInteractions/DecayMuon.h"
#include "ppInteractions/DecayChargedPion.h"
#include "ppInteractions/DecayNeutralPion.h"
#include "ppInteractions/DecayEtaMeson.h"

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

	ParticleDecay(bool photons = true, bool neutrinos = true, bool electrons = true, bool muons = true, double thinning = 0, const std::vector<int>& pList = std::vector<int>(), double limit = 0.1);
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

#endif