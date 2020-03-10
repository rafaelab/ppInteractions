%module(directors="1", threads="1", allprotected="1") NucleusNucleusInteraction

%include "std_string.i"
%include "std_iostream.i"
%include "std_container.i"
%include "exception.i"

%{
#include "CRPropa.h"
#include "ppInteractions/NucleusNucleusInteraction.h"
#include "ppInteractions/Constants.h"
#include "ppInteractions/ParticleDecay.h"
#include "ppInteractions/DecayMuon.h"
#include "ppInteractions/DecayEtaMeson.h"
#include "ppInteractions/DecayChargedPion.h"
#include "ppInteractions/DecayNeutralPion.h"
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/* include plugin parts to generate wrappers for */
%include "ppInteractions/NucleusNucleusInteraction.h"
%include "ppInteractions/Constants.h"
%include "ppInteractions/ParticleDecay.h"
%include "ppInteractions/DecayMuon.h"
%include "ppInteractions/DecayEtaMeson.h"
%include "ppInteractions/DecayChargedPion.h"
%include "ppInteractions/DecayNeutralPion.h"



