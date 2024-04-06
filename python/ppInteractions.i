%module(directors = "1", threads = "1", allprotected = "1") ppInteractions


%include "attribute.i"
%include "cpointer.i"
%include "exception.i"
%include "std_array.i"
%include "std_container.i"
%include "std_iostream.i"
%include "std_map.i"
%include "std_pair.i"
%include "std_set.i"
%include "std_string.i"
%include "std_unordered_map.i"
%include "std_vector.i"
%include "stl.i"
%include "typemaps.i"


%{
	#include "CRPropa.h"
%}

/* Import CRPropa and BSMPropa in wrapper */
%import (module = "crpropa") "crpropa.i"

%{
	#include "ppInteractions/Constants.h"
	#include "ppInteractions/NucleusNucleusInteraction.h"
	#include "ppInteractions/ParticleDecay.h"
	#include "ppInteractions/DecayMuon.h"
	#include "ppInteractions/DecayEtaMeson.h"
	#include "ppInteractions/DecayChargedPion.h"
	#include "ppInteractions/DecayNeutralPion.h"
%}


/* Include plugin parts to generate wrappers for */
%include "ppInteractions/Constants.h"
%include "ppInteractions/NucleusNucleusInteraction.h"
%include "ppInteractions/ParticleDecay.h"
%include "ppInteractions/DecayMuon.h"
%include "ppInteractions/DecayEtaMeson.h"
%include "ppInteractions/DecayChargedPion.h"
%include "ppInteractions/DecayNeutralPion.h"


/* Ignore list */
%ignore operator<<;
%ignore operator>>;
%ignore *::operator=;


/* Hide warnings */
#pragma SWIG nowarn=302,312,325,361,389,401,508,509





