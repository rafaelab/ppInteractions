#ifndef DECAYCHARGEDPION_H
#define DECAYCHARGEDPION_H

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Common.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Grid.h>
#include <crpropa/GridTools.h>
#include <crpropa/Candidate.h>


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


#endif