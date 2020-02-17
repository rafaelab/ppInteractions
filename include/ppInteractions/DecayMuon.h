#ifndef DECAYMUON_H
#define DECAYMUON_H

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Common.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Grid.h>
#include <crpropa/GridTools.h>
#include <crpropa/Candidate.h>
#include <crpropa/Vector3.h>

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


#endif