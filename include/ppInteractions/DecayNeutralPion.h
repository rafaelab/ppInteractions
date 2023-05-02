#ifndef DECAYNEUTRALPION_H
#define DECAYNEUTRALPION_H

#include <crpropa/Units.h>
#include <crpropa/Common.h>
#include <crpropa/ParticleID.h>
#include <crpropa/ParticleMass.h>
#include <crpropa/Random.h>
#include <crpropa/Candidate.h>
#include <crpropa/Vector3.h>

#include "ppInteractions/Constants.h"

class DecayNeutralPion {
	protected:
		bool havePhotons;
		double limit;
		double thinning;

	public:
		DecayNeutralPion(bool photons = true, double thinning = 0, double limit = 0.05);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setHavePhotons(bool photons);
		double lossLength(const double& lorentzFactor) const; 
		double energyFractionPhoton() const; 
		void performInteraction(crpropa::Candidate* candidate) const;
};


#endif