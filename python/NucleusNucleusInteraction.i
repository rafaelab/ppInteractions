%module(directors="1", threads="1", allprotected="1") NucleusNucleusInteraction

%include "exception.i"

%{
#include "CRPropa.h"
#include "ppInteractions/include/ppInteractions/NucleusNucleusInteraction.h"
#include "ppInteractions/include/ppInteractions/ParticleDecay.h"
#include "ppInteractions/include/ppInteractions/DecayMuon.h"
#include "ppInteractions/include/ppInteractions/DecayEtaMeson.h"
#include "ppInteractions/include/ppInteractions/DecayChargedPion.h"
#include "ppInteractions/include/ppInteractions/DecayNeutralPion.h"
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/* include plugin parts to generate wrappers for */
%include "include/ppInteractions/NucleusNucleusInteraction.h"
%include "ppInteractions/include/ppInteractions/ParticleDecay.h"
%include "ppInteractions/include/ppInteractions/DecayMuon.h"
%include "ppInteractions/include/ppInteractions/DecayEtaMeson.h"
%include "ppInteractions/include/ppInteractions/DecayChargedPion.h"
%include "ppInteractions/include/ppInteractions/DecayNeutralPion.h"



