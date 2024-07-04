#ifndef PPINTERACTIONS_DECAYMUON_H
#define PPINTERACTIONS_DECAYMUON_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

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

#include "ppInteractions/Constants.h"

class DecayMuon {
	protected:
		std::vector<double> probabilities;
		std::vector<double> electronFraction;
		std::vector<double> electronNeutrinoFraction;
		std::vector<double> muonNeutrinoFraction;
		std::string interactionTag;
		bool haveElectrons;
		bool haveNeutrinos;
		double limit;
		double thinning;
		
	public:
		DecayMuon(bool neutrinos = true, bool electrons = true, double thinning = 0, double limit = 0.05);
		void initSpectra();
		void setLimit(double limit);
		void setThinning(double thinning);
		void setHaveElectrons(bool electrons);
		void setHaveNeutrinos(bool neutrinos);
		void setInteractionTag(std::string tag);
		std::string getInteractionTag() const;
		double lossLength(const double& lorentzFactor) const; 
		double energyFractionElectron(crpropa::Random& random, const double& xmin = 0, const double& xmax = 1) const; 
		double energyFractionElectronNeutrino(crpropa::Random& random, const double& xmin = 0, const double& xmax = 1) const; 
		double energyFractionMuonNeutrino(crpropa::Random& random, const double& xmin = 0, const double& xmax = 1) const; 
		void performInteraction(crpropa::Candidate* candidate) const;
};


#endif