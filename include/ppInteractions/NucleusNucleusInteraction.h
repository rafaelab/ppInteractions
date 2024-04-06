#ifndef PPINTERACTIONS_NUCLEUSNUCLEUSINTERACTIONS_H
#define PPINTERACTIONS_NUCLEUSNUCLEUSINTERACTIONS_H

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
#include <crpropa/Candidate.h>
#include <crpropa/Vector3.h>

#include "ppInteractions/Constants.h"
#include "ppInteractions/ParticleDecay.h"


/**
 @class NucleusNucleusInteraction
 @brief Handles pp interactions following the parametrisations from:
 
 Kelner et al. Phys. Rev. D 74 (2006) 034018.
 Kafexhiu et al. Phys. Rev. D 90 (2014) 123014.

 Note that this is a hybrid approach: Kafexhiu is used only to model the gamma-ray
 spectrum due to the pi-0 decay, whereas Kelner+ is used to model everything else.
 */
class NucleusNucleusInteraction : public crpropa::Module {
		std::vector<double> logFraction;
		std::vector<double> logIncidentEnergy;
		std::vector<double> logProbabilities;
		std::vector<double> neutralPionFraction;
		std::vector<double> chargedPionFraction;
		std::vector<double> etaMesonFraction;
		std::vector<double> probabilities;
		crpropa::ref_ptr<crpropa::Grid1f> densityGrid;
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

	public:
		NucleusNucleusInteraction(double normMatterField = 1., double thinning = 0, double limit = 0.1);
		NucleusNucleusInteraction(crpropa::ref_ptr<crpropa::Grid1f> densityGrid, double normMatterField = 1., double thinning = 0, double limit = 0.1);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setFieldNorm(double normMatterField);
		void setIsDensityConstant(bool densityConstant);
		void setHavePhotons(bool photons);
		void setHaveElectrons(bool electorns);
		void setHaveNeutrinos(bool neutrinos);
		void initMesonSpectra();
		double crossSection(const double& energy) const;
		double energyFractionNeutralPion(const double& energy, const double& xmin = 0, const double& xmax = 1) const;
		double energyFractionChargedPion(const double& energy, const double& xmin = 0, const double& xmax = 1) const;
		double energyFractionEtaMeson(const double& energy, const double& xmin = 0, const double& xmax = 1) const;
		double lossLength(const int& id, const double& E) const; 
		double lossLength(const int& id, const double& E, const crpropa::Vector3d& position) const; 
		void process(crpropa::Candidate* candidate) const;
		void performInteraction(crpropa::Candidate* candidate) const;
};



#endif



