#ifndef DECAYCHARGEDPION_H
#define DECAYCHARGEDPION_H

#include <crpropa/Units.h>
#include <crpropa/Common.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Candidate.h>
#include <crpropa/Vector3.h>

#include "ppInteractions/Constants.h"


class DecayChargedPion {
	protected:
		bool haveMuons;
		bool haveNeutrinos;
		double limit;
		double thinning;

	public:
		DecayChargedPion(bool muons = true, bool neutrinos = true, double thinning = 0, double limit = 0.05);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setHaveMuons(bool muons);
		void setHaveNeutrinos(bool neutrinos);
		double lossLength(const double& lorentzFactor) const; 
		double energyFractionMuon() const; 
		void performInteraction(crpropa::Candidate* candidate) const;
};


#endif